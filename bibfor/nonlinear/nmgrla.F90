! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1504
!
subroutine nmgrla(FECell, FEBasis, FEQuad, option, typmod, &
                  imate, &
                  ndim, nno, npg, lgpg, &
                  compor, carcri, mult_comp, &
                  instam, instap, &
                  dispPrev, dispIncr, &
                  angmas, sigmPrev, sigmCurr, &
                  vim, vip, &
                  matsym, matuu, vectu, &
                  codret)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_eval_module
    use FE_mechanics_module
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/codere.h"
#include "asterfort/lcdetf.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmgrtg.h"
#include "asterfort/pk2sig.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
#include "FE_module.h"
!
    type(FE_Cell), intent(in) :: FECell
    type(FE_Quadrature), intent(in) :: FEQuad
    type(FE_basis), intent(in) :: FEBasis
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: typmod(*)
    integer, intent(in) :: imate
    integer, intent(in) :: ndim, nno, npg, lgpg
    character(len=16), intent(in) :: compor(*)
    real(kind=8), intent(in) :: carcri(*), angmas(*)
    character(len=16), intent(in) :: mult_comp
    real(kind=8), intent(in) :: instam, instap
    real(kind=8), intent(inout) :: dispPrev(ndim*nno), dispIncr(ndim*nno)
    real(kind=8), intent(inout) :: sigmPrev(2*ndim, npg), sigmCurr(2*ndim, npg)
    real(kind=8), intent(inout) :: vim(lgpg, npg), vip(lgpg, npg)
    aster_logical, intent(in) :: matsym
    real(kind=8), intent(inout) :: matuu(*)
    real(kind=8), intent(inout) :: vectu(ndim*nno)
    integer, intent(inout) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: C_PLAN, D_PLAN, AXIS
!
! Options: RIGI_MECA_TANG, RAPH_MECA and FULL_MECA - Large displacements/rotations (GREEN_LAGRANGE)
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  imate            : coded material address (JEVEUX)
! In  ndim             : dimension (2 ou 3)
! In  nno              : number of nodes
! In  npg              : number of Gauss integration point
! In  lgpg             : total length of vector for internal state variable
! In  ipoids           : Gauss point weight address (JEVEUX)
! In  ivf              : shape functions address (JEVEUX)
! In  idfde            : derivative of shape functions address (JEVEUX)
! In  carcri           : parameters for behaviour
! In  compor           : behaviour
! In  mult_comp        : multi-comportment (DEFI_COMPOR for PMF)
! In  instam           : time at beginning of time step
! In  instap           : time at end of time step
! IO  dispPrev        : displacements at beginning of time step
! IO  dispIncr        : displacements from beginning of time step
! IO  sigmPrev             : stresses at beginning of time step
! IO  sigmCurr             : stresses at end of time step
! IO  vim              : internal state variables at beginning of time step
! IO  vip              : internal state variables at end of time step
! In  matsym           : .true. if symmetric matrix
! IO  matuu            : matrix
! IO  vectu            : vector (internal forces)
! Out codret           : code for error
!
! --------------------------------------------------------------------------------------------------
!
!
    aster_logical :: lVect, lMatr, lSigm, lMatrPred, lPred
    integer :: kpg, ipoids, ivf, idfde
    integer :: cod(MAX_QP)
    real(kind=8) :: dsidep(6, 6), coorpg(3), BGSEval(3, MAX_BS)
    real(kind=8) :: fPrev(3, 3), fCurr(3, 3), fIncr(3, 3), gPrev(3, 3), gCurr(3, 3)
    real(kind=8) :: epsgPrev(6), epsgIncr(6), epsgCurr(6)
    real(kind=8) :: detfPrev, detfCurr
    real(kind=8) :: dispCurr(ndim*nno)
    real(kind=8) :: sigmPost(6), sigmPrep(6)
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    type(Behaviour_Integ) :: BEHinteg
!
! --------------------------------------------------------------------------------------------------
!
    cod = 0
    lSigm = L_SIGM(option)
    lVect = L_VECT(option)
    lMatr = L_MATR(option)
    lPred = L_PRED(option)
    lMatrPred = L_MATR_PRED(option)
    dispCurr = 0.d0
!
! - Initialisation of behaviour datastructure
!
    call behaviourInit(BEHinteg)
!
! - Prepare external state variables
!
    call elrefe_info(fami=FEQuad%fami, jpoids=ipoids, jvf=ivf, jdfde=idfde)
    call behaviourPrepESVAElem(carcri, typmod, &
                               nno, npg, ndim, &
                               ipoids, ivf, idfde, &
                               FECell%coorno(1:ndim, 1:nno), BEHinteg, &
                               dispPrev, dispIncr)
!
! - Update displacements
!
    dispCurr(:) = dispPrev(:)+dispIncr(:)
!
! - Loop on Gauss points
!
    do kpg = 1, FEQuad%nbQuadPoints
        coorpg = FEQuad%points_param(1:3, kpg)
        BGSEval = FEBasis%grad(coorpg, FEQuad%jacob(1:3, 1:3, kpg))

! ----- Kinematic - Previous strains
        gPrev = FEEvalGradMat(FEBasis, dispPrev, coorpg, BGSEval)
        fPrev = matG2F(gPrev)
        epsgPrev = matG2E(gPrev)

! ----- Kinematic - Current strains
        gCurr = FEEvalGradMat(FEBasis, dispCurr, coorpg, BGSEval)
        fCurr = matG2F(gCurr)
        epsgCurr = matG2E(gCurr)

! ----- Stresses: convert Cauchy to PK2
        call lcdetf(ndim, fPrev, detfPrev)
        call pk2sig(ndim, fPrev, detfPrev, sigmPrep, sigmPrev(1, kpg), -1)
        sigmPrep(4:2*ndim) = sigmPrep(4:2*ndim)*rac2

! ----- Compute behaviour
        sigmPost = 0
! ----- Check if the behavior law is MFRONT
        if (compor(RELA_NAME) == 'MFRONT') then
! --------- Compute the increment of f for MFRONT
            fIncr = fCurr - fPrev

            call nmcomp(BEHinteg, &
                        FEQuad%fami, kpg, 1, ndim, typmod, &
                        imate, compor, carcri, instam, instap, &
                        9, fPrev, fIncr, 6, sigmPrep, &
                        vim(1, kpg), option, angmas, &
                        sigmPost, vip(1, kpg), 36, dsidep, &
                        cod(kpg), mult_comp)
        else
! --------- Original behavior
            epsgIncr = epsgCurr - epsgPrev

            call nmcomp(BEHinteg, &
                        FEQuad%fami, kpg, 1, ndim, typmod, &
                        imate, compor, carcri, instam, instap, &
                        6, epsgPrev, epsgIncr, 6, sigmPrep, &
                        vim(1, kpg), option, angmas, &
                        sigmPost, vip(1, kpg), 36, dsidep, &
                        cod(kpg), mult_comp)
        endif
    
        if (cod(kpg) .eq. 1) goto 999
!        write (6,*) 'option = ',option
!        write (6,*) 'epsm   = ',epsgPrev
!        write (6,*) 'deps   = ',epsgIncr
!        write (6,*) 'sigm   = ',sigmPrep
!        if (lSigm) write(6,*) 'sigp   = ',sigmPost
!        if (lVect) write(6,*) 'vip    = ',vip(:,kpg)
!        if (lMatr) write(6,*) 'dsdiep = ',dsidep

! ----- Compute internal forces vector and rigidity matrix
        call nmgrtg(FEBasis, coorpg, FEQuad%weights(kpg), BGSEval, &
                    lVect, lMatr, lMatrPred, &
                    fPrev, fCurr, dsidep, sigmPrep, &
                    sigmPost, matsym, matuu, vectu)

! ----- Stresses: convert PK2 to Cauchy
        if (option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL') then
            call lcdetf(ndim, fCurr, detfCurr)
            call pk2sig(ndim, fCurr, detfCurr, sigmPost, sigmCurr(1, kpg), 1)
        end if
    end do
!
999 continue
!
! - Return code summary
!
    call codere(cod, npg, codret)
!
end subroutine
