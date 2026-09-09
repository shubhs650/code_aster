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
subroutine te0476(option, nomte)
!
    use HHO_type
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_Neumann_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_eval_module
    use HHO_utils_module
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/tefrep.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/writeVector.h"
!
    character(len=16), intent(in) :: option, nomte
!
!---------------------------------------------------------------------------------------------------
!
!  HHO METHODS
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          CORRESPONDANT A UN CHARGEMENT VOLUMIQUE POUR HHO
!          (LE CHARGEMENT PEUT ETRE DONNE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'CHAR_MECA_FF3D3D'
!                    'CHAR_MECA_FR3D3D'
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
    real(kind=8) :: rhs_forces(MSIZE_CELL_VEC), rhs(MSIZE_TDOFS_VEC), VoluValQP(3, MSIZE_QP_CELL)
    integer(kind=8) :: fbs, total_dofs, cbs, nbpara, idim
    integer(kind=8) :: j_time, j_forc
!
! -- Retrieve HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData, hhoQuad=hhoQuadCell)
!
    ASSERT(hhoQuadCell%nbQuadPoints <= MSIZE_QP_CELL)
!
    VoluValQP = 0.d0
    nompar(:) = 'XXXXXXXX'
    valpar(:) = 0.d0
!
! ---- Which option ?
!
    if (option .eq. 'CHAR_MECA_FF3D3D' .or. option .eq. 'CHAR_MECA_FF2D2D') then
!
! ---- Get Function Parameters
!
        if (hhocell%ndim == 3) then
            ASSERT(option .eq. 'CHAR_MECA_FF3D3D')
            call jevech('PFF3D3D', 'L', j_forc)
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (hhocell%ndim == 2) then
            ASSERT(option .eq. 'CHAR_MECA_FF2D2D')
            call jevech('PFF2D2D', 'L', j_forc)
            nbpara = 3
            nompar(1:2) = (/'X', 'Y'/)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ---- Time
!
        call jevech('PINSTR', 'L', j_time)
        nompar(nbpara) = 'INST'
        valpar(nbpara) = zr(j_time)
!
! ----- Evaluate the analytical function (FX,FY,FZ)
!
        do idim = 1, hhocell%ndim
            call hhoFuncFScalEvalQp(hhoQuadCell, zk8(j_forc-1+idim), nbpara, nompar, valpar, &
                                    hhocell%ndim, VoluValQP(idim, 1:MSIZE_QP_CELL))
        end do
!
    else if (option .eq. 'CHAR_MECA_FR2D2D' .or. option .eq. 'CHAR_MECA_FR3D3D') then
!
! ---- Get Forces
!
        if (hhocell%ndim == 3) then
            ASSERT(option .eq. 'CHAR_MECA_FR3D3D')
            call tefrep(option, 'PFR3D3D', j_forc)
        else if (hhocell%ndim == 2) then
            ASSERT(option .eq. 'CHAR_MECA_FR2D2D')
            call tefrep(option, 'PFR2D2D', j_forc)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ---- Compute the load at the quadrature points
!
        call hhoFuncRVecEvalCellQp(hhoCell, hhoQuadCell, zr(j_forc), VoluValQP)
!
    else
!
        ASSERT(ASTER_FALSE)
    end if
!
! ---- compute volumic load
!
    call hhoMecaVoluForces(hhoCell, hhoData, hhoQuadCell, VoluValQP, rhs_forces)
!
! ---- number of dofs
!
    call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    rhs = 0.d0
    rhs(total_dofs-cbs+1:total_dofs) = rhs_forces(1:cbs)
!
! ---- save result
!
    call writeVector('PVECTUR', total_dofs, rhs)
!
end subroutine
