! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine te0465(option, nomte)
!
    use HHO_type
    use HHO_basis_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_Neumann_module
    use HHO_init_module
    use HHO_eval_module
    use HHO_rhs_module
    use HHO_utils_module
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/readVector.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
!---------------------------------------------------------------------------------------------------
!
!  HHO METHODS
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN THERMIQUE
!          SUR UNE CELLULE POUR HHO
!          (LE CHARGEMENT PEUT ETRE DONNE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'CHAR_THER_SOUR_R'
!                    'CHAR_THER_SOUR_F'
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!---------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: maxpara = 4
    real(kind=8) :: valpar(maxpara)
    character(len=8) :: nompar(maxpara)
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_Quadrature) :: hhoQuadCell
    type(HHO_basis_cell) :: hhoBasisCell
    real(kind=8), dimension(MSIZE_CELL_SCAL) :: rhs_T, temp_T_curr
    real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: rhs
    real(kind=8) :: VoluValuesQP(MSIZE_QP_CELL)
    real(kind=8) :: theta, time_curr, tg
    integer(kind=8) :: fbs, nbpara, npg, faces_dofs, cbs, total_dofs
    integer(kind=8) :: j_time, j_sour, ipg, iret
!
!
! -- Get number of Gauss points
!
    call elrefe_info(fami='RIGI', npg=npg)
!
! -- Retrieve HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCell)
    call hhoBasisCell%initialize(hhoCell)
!
! ---- number of dofs
!
    call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    faces_dofs = total_dofs-cbs
!
    ASSERT(hhoQuadCell%nbQuadPoints <= MSIZE_QP_CELL)
!
    VoluValuesQP = 0.d0
    nompar(:) = 'XXXXXXXX'
    valpar(:) = 0.d0
!
    call jevech('PINSTR', 'L', j_time)
    time_curr = zr(j_time)
    theta = zr(j_time+2)
    ASSERT(theta < -0.5)
!
! ---- Which option ?
!
    if (option .eq. 'CHAR_THER_SOUR_R') then
!
! ----- Get real value COEF_H
!
        call jevech('PSOURCR', 'L', j_sour)
        VoluValuesQP(1:npg) = zr(j_sour-1+1:j_sour-1+npg)
!
    else if (option .eq. 'CHAR_THER_SOUR_F') then
        call jevech('PSOURCF', 'L', j_sour)
!
! ---- Get Function Parameters
!
        if (hhoCell%ndim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (hhoCell%ndim == 2) then
            nbpara = 3
            nompar(1:2) = (/'X', 'Y'/)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ---- Time +
!
        nompar(nbpara) = 'INST'
        valpar(nbpara) = time_curr
!
! ----- Evaluate the analytical function at T+
!
        call hhoFuncFScalEvalQp(hhoQuadCell, zk8(j_sour), nbpara, nompar, valpar, &
                                hhoCell%ndim, VoluValuesQP)
!
    else if (option .eq. 'CHAR_THER_SOURNL') then
        call jevech('PSOURNL', 'L', j_sour)
        if (zk8(j_sour) (1:7) .eq. '&FOZERO') goto 999
!
        temp_T_curr = 0.d0
        call readVector('PTEMPER', cbs, temp_T_curr, total_dofs-cbs)
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCell%nbQuadPoints
!
! --------- Evaluate temperature
!
            tg = hhoEvalScalCell2( &
                 hhoBasisCell, hhoData%cell_degree(), hhoQuadCell%points(1:3, ipg), temp_T_curr)
!
! --------- Evaluate source
!
            call fointe('FM', zk8(j_sour), 1, ['TEMP'], [tg], &
                        VoluValuesQP(ipg), iret)
        end do
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
! ---- compute surface load
!
    call hhoMakeRhsCellScal(hhoCell, hhoQuadCell, VoluValuesQP, hhoData%cell_degree(), rhs_T)
!
! ---- save result
!
999 continue
!
    rhs = 0.0
    rhs(faces_dofs+1:total_dofs) = rhs_T(1:cbs)
    call writeVector('PVECTTR', total_dofs, rhs)
!
!
end subroutine
