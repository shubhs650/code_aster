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
subroutine te0461(option, nomte)
!
    use HHO_type
    use HHO_basis_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_Neumann_module
    use HHO_init_module, only: hhoInfoInitFace
    use HHO_eval_module
    use HHO_utils_module
!
    implicit none
!
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
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
!          SUR UNE FACE POUR HHO
!          (LE CHARGEMENT PEUT ETRE DONNE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'CHAR_THER_ECHA_R'
!                    'CHAR_THER_ECHA_F'
!                    'CHAR_THER_FLUN_R'
!                    'CHAR_THER_FLUN_F'
!                    'CHAR_THER_FLUNL'
!                    'CHAR_THER_RAYO_R'
!                    'CHAR_THER_RAYO_F'
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
    type(HHO_Face) :: hhoFace
    type(HHO_Quadrature) :: hhoQuadFace
    type(HHO_basis_face) :: hhoBasisFace
    real(kind=8), dimension(MSIZE_FACE_SCAL) :: rhs, temp_F_curr
    real(kind=8) :: CoeffQP_curr(MSIZE_QP_FACE)
    real(kind=8) :: ParaQP_curr(MSIZE_QP_FACE)
    real(kind=8) :: NeumValuesQP(MSIZE_QP_FACE)
    real(kind=8) :: time_curr, theta, temp_eval_curr, tz0
    real(kind=8) :: sigmEner(MSIZE_QP_FACE), epsil(MSIZE_QP_FACE), rbid
    integer(kind=8) :: fbs, celldim, ipg, nbpara, npg
    integer(kind=8) :: j_time, j_coefh, j_para
!
!
! -- Get number of Gauss points
!
    call elrefe_info(fami='RIGI', npg=npg)
!
! -- Retrieve HHO informations
!
    call hhoInfoInitFace(hhoFace, hhoData, npg, hhoQuadFace)
    call hhoBasisFace%initialize(hhoFace)
!
! ---- number of dofs
!
    call hhoTherFaceDofs(hhoFace, hhoData, fbs)
!
    ASSERT(hhoQuadFace%nbQuadPoints <= MSIZE_QP_FACE)
!
    celldim = hhoFace%ndim+1
    CoeffQP_curr = 0.d0
    ParaQP_curr = 0.d0
    NeumValuesQP = 0.d0
    nompar(:) = 'XXXXXXXX'
    valpar(:) = 0.d0
!
    call jevech('PINSTR', 'L', j_time)
    time_curr = zr(j_time)
    theta = zr(j_time+2)
    ASSERT(theta < -0.5)
!
!
! ---- Which option ?
!
    if (option .eq. 'CHAR_THER_ECHA_R') then
!
! ----- Get real value COEF_H
!
        call jevech('PCOEFHR', 'L', j_coefh)
        CoeffQP_curr = zr(j_coefh)
!
! ----- Get real value Text
!
        call jevech('PT_EXTR', 'L', j_para)
        ParaQP_curr = zr(j_para)
!
    else if (option .eq. 'CHAR_THER_ECHA_F') then
        call jevech('PCOEFHF', 'L', j_coefh)
        call jevech('PT_EXTF', 'L', j_para)
!
! ---- Get Function Parameters
!
        if (celldim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (celldim == 2) then
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
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_coefh), nbpara, nompar, valpar, &
                                celldim, CoeffQP_curr)
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para), nbpara, nompar, valpar, &
                                celldim, ParaQP_curr)
!
    else if (option .eq. 'CHAR_THER_RAYO_R') then
!
! ----- Get real value (sigmEner, epsil, temp_inf)
!
        call jevech('PRAYONR', 'L', j_para)
        CoeffQP_curr = zr(j_para)*zr(j_para+1)
        ParaQP_curr = zr(j_para+2)
!
    else if (option .eq. 'CHAR_THER_RAYO_F') then
        call jevech('PRAYONF', 'L', j_para)
!
! ---- Get Function Parameters (sigmEner, epsil, temp_inf)
!
        if (celldim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (celldim == 2) then
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
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para), nbpara, nompar, valpar, &
                                celldim, sigmEner)
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para+1), nbpara, nompar, valpar, &
                                celldim, epsil)
        CoeffQP_curr = sigmEner*epsil
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para+2), nbpara, nompar, valpar, &
                                celldim, ParaQP_curr)
!
    else if (option .eq. 'CHAR_THER_FLUN_R') then
!
! ----- Get real value FLUXN
!
        call jevech('PFLUXNR', 'L', j_para)
        ParaQP_curr = zr(j_para)
!
    else if (option .eq. 'CHAR_THER_FLUN_F') then
        call jevech('PFLUXNF', 'L', j_para)
!
! ---- Get Function Parameters
!
        if (celldim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (celldim == 2) then
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
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para), nbpara, nompar, valpar, &
                                celldim, ParaQP_curr)
!
    else if (option .eq. 'CHAR_THER_FLUNL') then
        call jevech('PFLUXNL', 'L', j_para)
!
    else
!
        ASSERT(ASTER_FALSE)
    end if
!
    if (option(1:15) .eq. 'CHAR_THER_ECHA_') then
!
        call readVector('PTEMPER', fbs, temp_F_curr)
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
            temp_eval_curr = hhoEvalScalFace( &
                             hhoBasisFace, hhoData%face_degree(), hhoQuadFace%points(1:3, ipg), &
                             temp_F_curr)
!
            NeumValuesQP(ipg) = CoeffQP_curr(ipg)*(ParaQP_curr(ipg)-temp_eval_curr)
        end do
    else if (option(1:15) .eq. 'CHAR_THER_RAYO_') then
!
        tz0 = r8t0()
!
        call readVector('PTEMPER', fbs, temp_F_curr)
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
            temp_eval_curr = hhoEvalScalFace( &
                             hhoBasisFace, hhoData%face_degree(), hhoQuadFace%points(1:3, ipg), &
                             temp_F_curr)
!
            NeumValuesQP(ipg) = CoeffQP_curr(ipg)*((ParaQP_curr(ipg)+tz0)**4-(temp_eval_curr+tz0)&
                                                  &**4)
        end do
!
    else if (option(1:15) .eq. 'CHAR_THER_FLUN_') then
        NeumValuesQP = ParaQP_curr
    else if (option .eq. 'CHAR_THER_FLUNL') then
        call readVector('PTEMPER', fbs, temp_F_curr)
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
            temp_eval_curr = hhoEvalScalFace( &
                             hhoBasisFace, hhoData%face_degree(), hhoQuadFace%points(1:3, ipg), &
                             temp_F_curr)
            call foderi(zk8(j_para), temp_eval_curr, NeumValuesQP(ipg), rbid)
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
! ---- compute surface load
!
    call hhoTherNeumForces(hhoFace, hhoData, hhoQuadFace, NeumValuesQP, rhs)
!
! ---- save result
!
    call writeVector('PVECTTR', fbs, rhs)
!
end subroutine
