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
subroutine nmdlog(FECell, FEBasis, FEQuad, &
                  nno, npg, ndim, &
                  typmod, option, &
                  compor, carcri, multComp, &
                  BEHInteg, &
                  instam, instap, &
                  dispPrev, dispIncr, &
                  lgpg, sigmPrev, vim, &
                  sigmCurr, vip, &
                  matsym, fint, matuu, &
                  codret)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
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
#include "asterfort/elrefe_info.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmgrtg.h"
#include "asterfort/poslog.h"
#include "asterfort/prelog.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "MeshTypes_type.h"
#include "FE_module.h"
!
    type(FE_Cell), intent(in) :: FECell
    type(FE_Quadrature), intent(in) :: FEQuad
    type(FE_basis), intent(in) :: FEBasis
    integer(kind=8), intent(in) :: nno, npg, ndim, lgpg
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: compor(COMPOR_SIZE), multComp
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    real(kind=8) :: dispPrev(*), dispIncr(*), sigmPrev(2*ndim, npg)
    real(kind=8) :: vim(lgpg, npg), sigmCurr(2*ndim, npg), vip(lgpg, npg)
    real(kind=8) :: matuu(*), fint(ndim*nno)
    integer(kind=8) :: codret
    real(kind=8) :: instam, instap
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: D_PLAN, C_PLAN, AXIS
!           AXIS_SI, C_PLAN_SI (QUAD8), D_PLAN_SI (QUAD8)
!           3D
!           3D_SI (HEXA20)
!
! Options: RIGI_MECA_TANG, RAPH_MECA and FULL_MECA - GDEF_LOG
!
! --------------------------------------------------------------------------------------------------
!
! IN  FAMI    : FAMILLE DE POINTS DE GAUSS
! IN  OPTION  : OPTION DE CALCUL
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO     : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : PTR. POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME
! IN  IDFF    : PTR. DERIVEE DES FONCTIONS DE FORME ELEMENT DE REF.
! IN  GEOMI   : COORDONNEES DES NOEUDS (CONFIGURATION INITIALE)
! MEM DFF     : ESPACE MEMOIRE POUR LA DERIVEE DES FONCTIONS DE FORME
!               DIM :(NNO,3) EN 3D, (NNO,4) EN AXI, (NNO,2) EN D_PLAN
! IN  COMPOR  : COMPORTEMENT
! IN  MATE    : MATERIAU CODE
! IN  LGPG    : DIMENSION DU VECTEUR DES VAR. INTERNES POUR 1 PT GAUSS
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTM   : VALEUR DE L'INSTANT T-
! IN  INSTP   : VALEUR DE L'INSTANT T+
! IN  MATSYM  : .TRUE. SI MATRICE SYMETRIQUE
! IN  DEPLM   : DEPLACEMENT EN T-
! IN  DEPLD   : INCREMENT DE DEPLACEMENT ENTRE T- ET T+
! IN  SIGM    : CONTRAINTES DE CAUCHY EN T-
! IN  VIM     : VARIABLES INTERNES EN T-
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA_*)
! OUT FINT    : FORCES INTERIEURES (RAPH_MECA ET FULL_MECA_*)
! OUT MATUU   : MATR. DE RIGIDITE NON SYM. (RIGI_MECA_* ET FULL_MECA_*)
! OUT CODRET  : CODE RETOUR DE L'INTEGRATION DE LA LDC
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    aster_logical :: cplan
    aster_logical :: matsym, lintbo
    aster_logical :: lVect, lMatr, lSigm, lMatrPred, lCorr, lVari
    integer(kind=8) :: kpg, nddl, cod(MAX_QP), ivf
    integer(kind=8) :: iw, idff, iret
    real(kind=8) :: dtde(6, 6)
    real(kind=8) :: epslPrev(6)
    real(kind=8) :: fPrev(3, 3), fCurr(3, 3), dispCurr(3*MT_NNOMAX3D)
    real(kind=8) :: tlogPrev(6), tlogCurr(6), epslIncr(6)
    real(kind=8) :: gn(3, 3), lamb(3), logl(3)
    real(kind=8) :: gPrev(3, 3), gCurr(3, 3), coorpg(3), BGSEval(3, MAX_BS_CG)
    real(kind=8) :: dsidep(6, 6), pk2Curr(6), pk2Prev(6)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    cplan = typmod(1) .eq. 'C_PLAN'
    lSigm = L_SIGM(option)
    lVari = L_VARI(option)
    lVect = L_VECT(option)
    lMatr = L_MATR(option)
    lMatrPred = L_MATR_PRED(option)
    lCorr = L_CORR(option)
    nddl = ndim*nno
    ASSERT(nno .le. 27)
    ASSERT(compor(PLANESTRESS) .ne. 'DEBORST')

! - Finite element parameters
    ASSERT(BEHInteg%materPara%schemePara%fami .eq. FEQuad%fami)
    call elrefe_info(fami=FEQuad%fami, jpoids=iw, jvf=ivf, jdfde=idff)

! - Prepare external state variables (geometry)
    call behaviourPrepESVAGeom(nno, npg, ndim, &
                               iw, ivf, idff, &
                               FECell%coorno(1:ndim, 1:nno), BEHInteg, &
                               dispPrev, dispIncr)

! - Update configuration
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, dispPrev, b_incx, dispCurr, b_incy)
    if (lCorr) then
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, dispIncr, b_incx, dispCurr, &
                   b_incy)
    end if

! - Loop on Gauss points
    lintbo = ASTER_FALSE
    cod = 0
    do kpg = 1, npg
        coorpg = FEQuad%points_param(1:3, kpg)
        BGSEval = FEBasis%grad(coorpg, FEQuad%jacob(1:3, 1:3, kpg))

! ----- Kinematic - Previous strains
        gPrev = FEEvalGradMat(FEBasis, dispPrev, coorpg, BGSEval)
        fPrev = matG2F(gPrev)

! ----- Kinematic - Current strains
        gCurr = FEEvalGradMat(FEBasis, dispCurr, coorpg, BGSEval)
        fCurr = matG2F(gCurr)

! ----- Pre-treatment of kinematic quantities
        call prelog(ndim, lgpg, vim(1, kpg), gn, lamb, &
                    logl, fPrev, fCurr, epslPrev, epslIncr, &
                    tlogPrev, lCorr, cod(kpg))
        if (cod(kpg) .ne. 0) then
            goto 999
        end if

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- Integrator
        cod(kpg) = 0
        dtde = 0.d0
        tlogCurr = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    instam, instap, &
                    compor, carcri, multComp, &
                    6, epslPrev, epslIncr, &
                    6, tlogPrev, &
                    vim(1, kpg), &
                    tlogCurr, vip(1, kpg), &
                    36, dtde, cod(kpg))
        if (cod(kpg) .eq. 1) then
            goto 999
        end if
        if (cod(kpg) .eq. 4) then
            lintbo = .true.
        end if

! ----- Post-treatment of sthenic quantities
        call poslog(BEHInteg, &
                    lCorr, lMatr, lSigm, lVari, &
                    tlogPrev, tlogCurr, fPrev, &
                    lgpg, vip(1, kpg), ndim, &
                    fCurr, dtde, sigmPrev(1, kpg), cplan, &
                    instap, gn, &
                    lamb, logl, sigmCurr(1, kpg), dsidep, pk2Prev, &
                    pk2Curr, iret)
        if (iret .eq. 1) then
            cod(kpg) = 1
            goto 999
        end if

! ----- Compute internal forces and matrix
        call nmgrtg(FEBasis, coorpg, FEQuad%weights(kpg), BGSEval, lVect, &
                    lMatr, lMatrPred, fPrev, fCurr, dsidep, &
                    pk2Prev, pk2Curr, matsym, matuu, fint)
    end do
    if (lintbo) then
        cod(1) = 4
    end if
!
999 continue

! - Return code summary
    call codere(cod, npg, codret)
!
end subroutine
