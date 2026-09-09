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
! aslint: disable=W1504
!
subroutine nmplxd(FECell, FEBasis, FEQuad, &
                  nno, npg, ndim, &
                  typmod, option, &
                  compor, carcri, multComp, &
                  BEHInteg, &
                  instam, instap, &
                  dispPrev, dispIncr, &
                  lgpg, sigmPrev, vim, &
                  sigmCurr, vip, &
                  matsym, matuu, vectu, &
                  codret)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_eval_module
    use FE_mechanics_module
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/crirup.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/nmcomp.h"
#include "FE_module.h"
!
    type(FE_Cell), intent(in) :: FECell
    type(FE_Quadrature), intent(in) :: FEQuad
    type(FE_basis), intent(in) :: FEBasis
    integer(kind=8), intent(in) :: nno, npg, ndim
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: compor(COMPOR_SIZE), multComp
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    real(kind=8), intent(in) :: instam, instap
    real(kind=8), intent(inout) :: dispPrev(ndim, nno), dispIncr(ndim, nno)
    integer(kind=8), intent(in) :: lgpg
    real(kind=8), intent(inout) :: sigmPrev(2*ndim, npg), vim(lgpg, npg)
    real(kind=8), intent(inout) :: sigmCurr(2*ndim, npg), vip(lgpg, npg)
    aster_logical, intent(in) :: matsym
    real(kind=8), intent(inout) :: matuu(*), vectu(ndim, nno)
    integer(kind=8), intent(inout) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D / D_PLAN / AXIS
!
! Options: RIGI_MECA_TANG, RAPH_MECA and FULL_MECA - Hypoelasticity (PETIT/PETIT_REAC)
!
! --------------------------------------------------------------------------------------------------
!
! IN  NNO     : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  POIDSG  : POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME
! IN  DFDE    : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  DFDK    : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  IMATE   : MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  CARCRI  : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT PRECEDENT
! IN  INSTAP  : INSTANT DE CALCUL
! IN  DEPLM   : DEPLACEMENT A L'INSTANT PRECEDENT
! IN  DEPLP   : INCREMENT DE DEPLACEMENT
! IN  ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT DFDI    : DERIVEE DES FONCTIONS DE FORME  AU DERNIER PT DE GAUSS
! OUT DEF     : PRODUIT DER. FCT. FORME PAR F   AU DERNIER PT DE GAUSS
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
! OUT MATUU   : MATRICE DE RIGIDITE PROFIL (RIGI_MECA_TANG ET FULL_MECA)
! OUT VECTU   : FORCES NODALES (RAPH_MECA ET FULL_MECA)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    aster_logical :: lVect, lMatr, lSigm
    integer(kind=8) :: kpg, i_tens, ipoids, ivf, idfde
    integer(kind=8) :: cod(MAX_QP)
    real(kind=8) :: BGSEval(3, MAX_BS_CG)
    real(kind=8) :: def(6, MAX_BS_CG, 3)
    real(kind=8) :: coorpg(3)
    real(kind=8) :: eps(6), deps(6)
    real(kind=8) :: dsidep(6, 6), sigmPost(6), sigmPrep(6)
!
! --------------------------------------------------------------------------------------------------
!
    cod = 0

! - Finite element parameters
    ASSERT(BEHInteg%materPara%schemePara%fami .eq. FEQuad%fami)
    call elrefe_info(fami=FEQuad%fami, jpoids=ipoids, jvf=ivf, jdfde=idfde)

! - Prepare external state variables (geometry)
    call behaviourPrepESVAGeom(nno, npg, ndim, &
                               ipoids, ivf, idfde, &
                               FECell%coorno(1:ndim, 1:nno), BEHInteg, &
                               dispPrev, dispIncr)

! - Quantities to compute
    lSigm = L_SIGM(option)
    lVect = L_VECT(option)
    lMatr = L_MATR(option)

! - Loop on Gauss points
    do kpg = 1, npg
        coorpg = FEQuad%points_param(1:3, kpg)
        BGSEval = FEBasis%grad(coorpg, FEQuad%jacob(1:3, 1:3, kpg))

! ----- Kinematic - Previous strains
        eps = FEEvalGradSymMat(FEBasis, dispPrev, coorpg, BGSEval)

! ----- Kinematic - Increment of strains
        deps = FEEvalGradSymMat(FEBasis, dispIncr, coorpg, BGSEval)

! ----- Kinematic - Product [B]
        call FEMatB(FEBasis, coorpg, BGSEval, def)

! ----- Prepare stresses
        do i_tens = 1, 3
            sigmPrep(i_tens) = sigmPrev(i_tens, kpg)
        end do
        do i_tens = 4, 2*ndim
            sigmPrep(i_tens) = sigmPrev(i_tens, kpg)*rac2
        end do

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- Integrator
        sigmPost = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    instam, instap, &
                    compor, carcri, multComp, &
                    6, eps, deps, &
                    6, sigmPrep, &
                    vim(1, kpg), &
                    sigmPost, vip(1, kpg), &
                    36, dsidep, &
                    cod(kpg))
        if (cod(kpg) .eq. 1) then
            goto 999
        end if

! ----- Rigidity matrix
        if (lMatr) then
            call FEStiffJacoVectSymAdd(FEBasis, def, FEQuad%weights(kpg), dsidep, matsym, matuu)
        end if

! ----- Internal forces
        if (lVect) then
            call FEStiffResiVectSymAdd(FEBasis, def, FEQuad%weights(kpg), sigmPost, vectu)
        end if

! ----- Cauchy stresses
        if (lSigm .or. option .eq. 'RIGI_MECA_IMPLEX') then
            do i_tens = 1, 3
                sigmCurr(i_tens, kpg) = sigmPost(i_tens)
            end do
            do i_tens = 4, 2*ndim
                sigmCurr(i_tens, kpg) = sigmPost(i_tens)/rac2
            end do
        end if
    end do

! - For POST_ITER='CRIT_RUPT'
    if (carcri(13) .gt. 0.d0) then
        call crirup(BEHInteg%materPara, &
                    ndim, npg, lgpg, &
                    option, sigmCurr, vip, vim, &
                    instam, instap)
    end if
!
999 continue

! - Return code summary
    call codere(cod, FEQuad%nbQuadPoints, codret)
!
end subroutine
