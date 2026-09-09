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
subroutine te0075(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_rhs_module
    use FE_eval_module
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/writeVector.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES CHARGEMENTS ELEMENTAIRES

!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
    type(FE_Skin) :: FESkin
    type(FE_Quadrature) :: FEQuad
    type(FE_Basis) :: FEBasis
!
    integer(kind=8), parameter :: nbres = 4
    character(len=8) :: nompar(nbres)
    real(kind=8) :: valpar(nbres), theta, time_curr, time_prev
    real(kind=8) :: rhs(MAX_BS_CG), valQP(MAX_QP)
    real(kind=8) :: valQPC(MAX_BS_CG), valQPP(MAX_QP)
    real(kind=8) :: para1, para2, para3, tz0, tpg, normal(3)
    integer(kind=8) :: kp, itemps, ipara, icode, ipara2
    real(kind=8), pointer :: tempi(:) => null()
!
    call FESkin%init()
    call FEQuad%initFace(FESkin, "RIGI")
    call FEBasis%initFace(FESkin)
!
    theta = 1.d0
    tz0 = r8t0()
!
    if (option .ne. "CHAR_THER_FLUN_R") then
        call jevech('PINSTR', 'L', itemps)
!
        theta = zr(itemps+2)
        time_curr = zr(itemps)
        time_prev = zr(itemps)-zr(itemps+1)
    end if
!
    nompar(1:3) = ['X', 'Y', 'Z']
    nompar(4) = 'INST'
!
    valQPC = 0.d0
    valQPP = 0.d0
!
    if (option == "CHAR_THER_FLUN_F") then
        call jevech('PFLUXNF', 'L', ipara)
!
        do kp = 1, FEQuad%nbQuadPoints
            valpar(1:3) = FEQuad%points(1:3, kp)
            valpar(4) = time_curr
            call fointe('FM', zk8(ipara), 4, nompar, valpar, valQPC(kp), icode)
            if (theta > -0.5d0) then
                valpar(4) = time_prev
                call fointe('FM', zk8(ipara), 4, nompar, valpar, valQPP(kp), icode)
            end if
        end do
    elseif (option == "CHAR_THER_FLUN_R") then
        call jevech('PFLUXNR', 'L', ipara)
        valQPC = zr(ipara)
    else if (option == "CHAR_THER_FLUX_F") then
!
        call jevech('PFLUXVF', 'L', ipara)
        theta = -1.d0
!
        do kp = 1, FEQuad%nbQuadPoints
!
            valpar(1:3) = FEQuad%points(1:3, kp)
            valpar(4) = time_curr
            ! FLUX_X
            call fointe('FM', zk8(ipara), 4, nompar, valpar, para1, icode)
            ! FLUX_Y
            call fointe('FM', zk8(ipara+1), 4, nompar, valpar, para2, icode)
            ! FLUX_Z
            if (FESkin%ndim+1 == 3) then
                call fointe('FM', zk8(ipara+2), 4, nompar, valpar, para3, icode)
            else
                para3 = 0.d0
            end if
            normal = FESkin%normal(FEQuad%points_param(1:2, kp))
            !
            ! FLUX.NORMAL
            valQPC(kp) = normal(1)*para1+normal(2)*para2+normal(3)*para3
        end do
    else if (option == "CHAR_THER_RAYO_F") then
!
        call jevech('PRAYONF', 'L', ipara)
        call jevech('PTEMPER', 'L', vr=tempi)
!
        do kp = 1, FEQuad%nbQuadPoints
            tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuad%points_param(1:3, kp))
!
            valpar(1:3) = FEQuad%points(1:3, kp)
            valpar(4) = time_curr
            ! SIGM
            call fointe('FM', zk8(ipara), 4, nompar, valpar, para1, icode)
            ! EPS
            call fointe('FM', zk8(ipara+1), 4, nompar, valpar, para2, icode)
            ! TPF
            call fointe('FM', zk8(ipara+2), 4, nompar, valpar, para3, icode)
            !
            if (theta > -0.5d0) then
                valpar(4) = time_prev
                ! SIGM
                call fointe('FM', zk8(ipara), 4, nompar, valpar, para1, icode)
                ! EPS
                call fointe('FM', zk8(ipara+1), 4, nompar, valpar, para2, icode)
                ! TPF
                call fointe('FM', zk8(ipara+2), 4, nompar, valpar, para3, icode)
                !
                valQPC(kp) = para1*para2*((para3+tz0)**4)
                valQPP(kp) = para1*para2*((para3+tz0)**4-(tpg+tz0)**4)
            else
                ! SIGM * EPS * ((TPF+T0)^4-(T+T0)^4)
                valQPC(kp) = para1*para2*((para3+tz0)**4-(tpg+tz0)**4)
            end if
        end do
    else if (option == "CHAR_THER_RAYO_R") then
!
        call jevech('PRAYONR', 'L', ipara)
        call jevech('PTEMPER', 'L', vr=tempi)
!
        do kp = 1, FEQuad%nbQuadPoints
            tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuad%points_param(1:3, kp))
            ! SIGM
            para1 = zr(ipara)
            ! EPS
            para2 = zr(ipara+1)
            ! TPF
            para3 = zr(ipara+2)
            !
            if (theta > -0.5d0) then
                valQPC(kp) = para1*para2*((para3+tz0)**4)
                valQPP(kp) = para1*para2*((para3+tz0)**4-(tpg+tz0)**4)
            else
                ! SIGM * EPS * ((TPF+T0)^4-(T+T0)^4)
                valQPC(kp) = para1*para2*((para3+tz0)**4-(tpg+tz0)**4)
            end if
        end do
    else if (option == "CHAR_THER_ECHA_F") then
        !
        call jevech('PT_EXTF', 'L', ipara)
        call jevech('PCOEFHF', 'L', ipara2)
        call jevech('PTEMPER', 'L', vr=tempi)
!
        do kp = 1, FEQuad%nbQuadPoints
            tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuad%points_param(1:3, kp))
!
            valpar(1:3) = FEQuad%points(1:3, kp)
            valpar(4) = time_curr
            ! TEXT
            call fointe('FM', zk8(ipara), 4, nompar, valpar, para1, icode)
            ! COEFH
            call fointe('FM', zk8(ipara2), 4, nompar, valpar, para2, icode)
            !
            if (theta > -0.5d0) then
                valQPC(kp) = para2*para1
                valpar(4) = time_prev
                ! TEXT
                call fointe('FM', zk8(ipara), 4, nompar, valpar, para1, icode)
                ! COEFH
                call fointe('FM', zk8(ipara2), 4, nompar, valpar, para2, icode)
                !
                valQPP(kp) = para2*(para1-tpg)
            else
                ! COEF_H*(TEXT-T)
                valQPC(kp) = para2*(para1-tpg)
            end if
        end do
    else if (option == "CHAR_THER_ECHA_R") then
!
        call jevech('PT_EXTR', 'L', ipara)
        call jevech('PCOEFHR', 'L', ipara2)
        call jevech('PTEMPER', 'L', vr=tempi)
!
        do kp = 1, FEQuad%nbQuadPoints
            tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuad%points_param(1:3, kp))
            ! TEXT
            para1 = zr(ipara)
            ! COEFH
            para2 = zr(ipara2)
            !
            if (theta > -0.5d0) then
                valQPC(kp) = para2*para1
                valQPP(kp) = para2*(para1-tpg)
            else
                ! COEF_H*(TEXT-T)
                valQPC(kp) = para2*(para1-tpg)
            end if
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
    do kp = 1, FEQuad%nbQuadPoints
        if (theta > -0.5d0) then
            valQP(kp) = theta*valQPC(kp)+(1.0d0-theta)*valQPP(kp)
        else
            valQP(kp) = valQPC(kp)
        end if
    end do
!
    call FeMakeRhsScal(FEQuad, FEBasis, valQP, rhs)
    call writeVector("PVECTTR", FEBasis%size, rhs)
!
end subroutine
