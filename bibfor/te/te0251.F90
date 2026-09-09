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
subroutine te0251(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_mass_module
    use FE_eval_module
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/foderi.h"
#include "asterfort/jevech.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES TANGENTES ELEMENTAIRES
!                          EN THERMIQUE
!
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
    real(kind=8) :: valpar(nbres), time_curr
    real(kind=8) :: mass(MAX_BS_CG, MAX_BS_CG), valQP(MAX_QP)
    real(kind=8) :: para1, para2, tz0, tpg, rbid
    integer(kind=8) :: kp, itemps, ipara, icode, ipara2
    real(kind=8), pointer :: tempi(:) => null()
!
    call FESkin%init()
    call FEQuad%initFace(FESkin, "RIGI")
    call FEBasis%initFace(FESkin)
!
    tz0 = r8t0()
!
    call jevech('PINSTR', 'L', itemps)
    time_curr = zr(itemps)
!
    nompar(1:3) = ['X', 'Y', 'Z']
    nompar(4) = 'INST'
!
    valQP = 0.d0
!
    if (option == "MTAN_THER_FLUXNL") then
        call jevech('PFLUXNL', 'L', ipara)
        call jevech('PTEMPEI', 'L', vr=tempi)
!
        do kp = 1, FEQuad%nbQuadPoints
            tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuad%points_param(1:3, kp))
            call foderi(zk8(ipara), tpg, rbid, para1)
            ! - d alpha/dT
            valQP(kp) = -para1
        end do
    else if (option == "MTAN_THER_RAYO_F") then
!
        call jevech('PRAYONF', 'L', ipara)
        call jevech('PTEMPEI', 'L', vr=tempi)
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
            !
            ! 4*SIGM * EPS * (T+T0)^3
            valQP(kp) = 4.d0*para1*para2*(tpg+tz0)**3
        end do
    else if (option == "MTAN_THER_RAYO_R") then
!
        call jevech('PRAYONR', 'L', ipara)
        call jevech('PTEMPEI', 'L', vr=tempi)
!
        ! SIGM
        para1 = zr(ipara)
        ! EPS
        para2 = zr(ipara+1)
!
        do kp = 1, FEQuad%nbQuadPoints
            tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuad%points_param(1:3, kp))
            ! 4*SIGM * EPS * (T+T0)^3
            valQP(kp) = 4.d0*para1*para2*(tpg+tz0)**3
        end do
    else if (option == "RIGI_THER_ECHA_F") then
!
        call jevech('PCOEFHF', 'L', ipara2)
!
        do kp = 1, FEQuad%nbQuadPoints
            valpar(1:3) = FEQuad%points(1:3, kp)
            valpar(4) = time_curr
            ! COEFH
            call fointe('FM', zk8(ipara2), 4, nompar, valpar, valQP(kp), icode)
        end do
    else if (option == "RIGI_THER_ECHA_R") then
!
        call jevech('PCOEFHR', 'L', ipara2)
        ! COEFH
        valQP(1:FEQuad%nbQuadPoints) = zr(ipara2)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call FEMassMatScal(FEQuad, FEBasis, mass, valQP)
    call writeMatrix("PMATTTR", FEBasis%size, FEBasis%size, ASTER_TRUE, mass)
!
end subroutine
