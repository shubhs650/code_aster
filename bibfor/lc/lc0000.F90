! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1306,W1501,W1504,C1505
!
subroutine lc0000(BEHinteg, &
                  fami, kpg, ksp, ndim, typmod, &
                  l_epsi_varc, imate, materi, compor, mult_comp, &
                  carcri, instam, instap, neps, epsm_tot, &
                  deps_tot, nsig, sigm_all, vim, option, &
                  angmas, numlc, sigp, vip, &
                  ndsde, dsidep, icomp, nvi_all, codret)
!
    use calcul_module, only: calcul_status
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lc0001.h"
#include "asterfort/lc0002.h"
#include "asterfort/lc0003.h"
#include "asterfort/lc0004.h"
#include "asterfort/lc0007.h"
#include "asterfort/lc0008.h"
#include "asterfort/lc0009.h"
#include "asterfort/lc0015.h"
#include "asterfort/lc0016.h"
#include "asterfort/lc0017.h"
#include "asterfort/lc0018.h"
#include "asterfort/lc0019.h"
#include "asterfort/lc0021.h"
#include "asterfort/lc0022.h"
#include "asterfort/lc0023.h"
#include "asterfort/lc0024.h"
#include "asterfort/lc0025.h"
#include "asterfort/lc0026.h"
#include "asterfort/lc0028.h"
#include "asterfort/lc0029.h"
#include "asterfort/lc0030.h"
#include "asterfort/lc0031.h"
#include "asterfort/lc0032.h"
#include "asterfort/lc0033.h"
#include "asterfort/lc0034.h"
#include "asterfort/lc0035.h"
#include "asterfort/lc0036.h"
#include "asterfort/lc0040.h"
#include "asterfort/lc0042.h"
#include "asterfort/lc0050.h"
#include "asterfort/lc0054.h"
#include "asterfort/lc0055.h"
#include "asterfort/lc0058.h"
#include "asterfort/lc0059.h"
#include "asterfort/lc0060.h"
#include "asterfort/lc0062.h"
#include "asterfort/lc0075.h"
#include "asterfort/lc0076.h"
#include "asterfort/lc0077.h"
#include "asterfort/lc0078.h"
#include "asterfort/lc0079.h"
#include "asterfort/lc0120.h"
#include "asterfort/lc0137.h"
#include "asterfort/lc0145.h"
#include "asterfort/lc0152.h"
#include "asterfort/lc0165.h"
#include "asterfort/lc0166.h"
#include "asterfort/lc0167.h"
#include "asterfort/lc0168.h"
#include "asterfort/lc0169.h"
#include "asterfort/lc1002.h"
#include "asterfort/lc1015.h"
#include "asterfort/lc1037.h"
#include "asterfort/lc1137.h"
#include "asterfort/lc2001.h"
#include "asterfort/lc2002.h"
#include "asterfort/lc2036.h"
#include "asterfort/lc3053.h"
#include "asterfort/lc4047.h"
#include "asterfort/lc6036.h"
#include "asterfort/lc6046.h"
#include "asterfort/lc6057.h"
#include "asterfort/lc6075.h"
#include "asterfort/lc6076.h"
#include "asterfort/lc7010.h"
#include "asterfort/lc7011.h"
#include "asterfort/lc7013.h"
#include "asterfort/lc7040.h"
#include "asterfort/lc7041.h"
#include "asterfort/lc7043.h"
#include "asterfort/lc7045.h"
#include "asterfort/lc7046.h"
#include "asterfort/lc7047.h"
#include "asterfort/lc7048.h"
#include "asterfort/lc7049.h"
#include "asterfort/lc7051.h"
#include "asterfort/lc7056.h"
#include "asterfort/lc7058.h"
#include "asterfort/lc7077.h"
#include "asterfort/lc8028.h"
#include "asterfort/lc8029.h"
#include "asterfort/lc8057.h"
#include "asterfort/lc8146.h"
#include "asterfort/lc8331.h"
#include "asterfort/lcvisc.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHinteg
    integer :: imate, ndim, nvi_all, kpg, ksp
    aster_logical, intent(in) :: l_epsi_varc
    integer :: neps, nsig, ndsde
    real(kind=8) :: carcri(*), angmas(3)
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm_tot(neps), deps_tot(neps)
    real(kind=8) :: sigm_all(nsig), sigp(nsig)
    real(kind=8) :: vim(nvi_all), vip(nvi_all)
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    character(len=16) :: compor(*), option
    character(len=8), intent(in) :: materi
    character(len=16), intent(in) :: mult_comp
    character(len=8) :: typmod(*)
    character(len=*) :: fami
    integer :: icomp
    integer :: numlc
    integer :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     INTEGRATION DES LOIS DE COMPORTEMENT NON LINEAIRE POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES OU GRANDES DEFORMATIONS
!
! --------------------------------------------------------------------------------------------------
!
! IO  BEHinteg         : parameters for integration of behaviour
! IN  FAMI,KPG,KSP  : FAMILLE ET NUMERO DU (SOUS)POINT DE GAUSS
!     NDIM    : DIMENSION DE L'ESPACE
!               3 : 3D , 2 : D_PLAN ,AXIS OU  C_PLAN
!     TYPMOD(2): MODELISATION ex: 1:3D, 2:INCO
!     IMATE   : ADRESSE DU MATERIAU CODE
!     COMPOR  : COMPORTEMENT :  (1) = TYPE DE RELATION COMPORTEMENT
!                               (2) = NB VARIABLES INTERNES / PG
!                               (3) = HYPOTHESE SUR LES DEFORMATIONS
!                               (4) etc... (voir grandeur COMPOR)
!     CRIT    : CRITERES DE CONVERGENCE LOCAUX (voir grandeur CARCRI)
!     INSTAM  : INSTANT DU CALCUL PRECEDENT
!     INSTAP  : INSTANT DU CALCUL
!     NEPS    : NOMBRE DE CMP DE EPSM ET DEPS (SUIVANT MODELISATION)
!     EPSM    : DEFORMATIONS A L'INSTANT DU CALCUL PRECEDENT
!     DEPS    : INCREMENT DE DEFORMATION TOTALE :
!                DEPS(T) = DEPS(MECANIQUE(T)) + DEPS(DILATATION(T))
!     NSIG    : NOMBRE DE CMP DE SIGM ET SIGP (SUIVANT MODELISATION)
!     SIGM_ALL: CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
!     VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
!     OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
!     ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM),
!               + UN REEL QUI VAUT 0 SI NAUTIQUIES OU 2 SI EULER
!               + LES 3 ANGLES D'EULER
!     NUMLC   : NUMERO DE LOI DE COMPORTEMENT ISSUE DU CATALOGUE DE LC
!     ICOMP   : COMPTEUR DE REDECOUPAGE PRODUIT PAR REDECE
!     NVI_ALL : NOMBRE DE VARIABLES INTERNES DU POINT D'INTEGRATION
!
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! VAR VIP     : VARIABLES INTERNES
!                IN  : ESTIMATION (ITERATION PRECEDENTE OU LAG. AUGM.)
!                OUT : EN T+
!     NDSDE   : DIMENSION DE DSIDEP
!     DSIDEP  : OPERATEUR TANGENT DSIG/DEPS OU DSIG/DF
!     CODRET  : CODE RETOUR LOI DE COMPORMENT :
!               CODRET=0 : TOUT VA BIEN
!               CODRET=1 : ECHEC DANS L'INTEGRATION DE LA LOI
!               CODRET=2 : LOI OK MAIS UNE CONDITION INTERDIT LA CONVERGENCE GLOBALE
!
! PRECISIONS :
! -----------
!  LES TENSEURS ET MATRICES SONT RANGES DANS L'ORDRE :
!         XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ
!
! -SI DEFORMATION = SIMO_MIEHE
!   EPSM(3,3)    GRADIENT DE LA TRANSFORMATION EN T-
!   DEPS(3,3)    GRADIENT DE LA TRANSFORMATION DE T- A T+
!
!  OUTPUT SI RESI (RAPH_MECA, FULL_MECA_*)
!   VIP      VARIABLES INTERNES EN T+
!   SIGP(6)  CONTRAINTE DE KIRCHHOFF EN T+ RANGES DANS L'ORDRE
!         XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ
!
!  OUTPUT SI RIGI (RIGI_MECA_*, FULL_MECA_*)
!   DSIDEP(6,3,3) MATRICE TANGENTE D(TAU)/D(FD) * (FD)T
!                 (AVEC LES RACINES DE 2)
!
! -SINON (DEFORMATION = PETIT OU PETIT_REAC OU GDEF_...)
!   EPSM(6), DEPS(6)  SONT LES DEFORMATIONS (LINEARISEES OU GREEN OU ..)
!
! ----------------------------------------------------------------------
!
!    ATTENTION  VIM    VARIABLES INTERNES A T MODIFIEES SI REDECOUPAGE
!       ----------------------------------------------------------------

!     VARIABLES LOCALES POUR LE REDECOUPAGE DU PAS DE TEMPS
!             TD      INSTANT T
!             TF      INSTANT T+DT
!             TEMD    TEMPERATURE A T
!             TEMF    TEMPERATURE A T+DT
!             DEPS    INCREMENT DE DEFORMATION TOTALE
!             VD      VARIABLES INTERNES A T    + INDICATEUR ETAT T
!             DSIDEPLO MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!             ICOMP           COMPTEUR POUR LE REDECOUPAGE DU PAS DE
!                                  TEMPS
!             RETURN1 EN CAS DE NON CONVERGENCE LOCALE
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter             :: rac2 = sqrt(2.d0)
    real(kind=8), dimension(6), parameter:: r2 = [1.d0, 1.d0, 1.d0, rac2, rac2, rac2]
    integer, parameter                 :: nvi_regu_visc = 8
    integer, parameter                 :: nvi_gdef_log = 6
    character(len=16) :: defo_ldc, defo_comp, regu_visc, postIncr
    aster_logical :: l_pred, l_czm, l_defo_meca, l_large, lVari, lSigm, lMatr, l_grad_vari
    aster_logical :: l_regu_visc, l_gdef_log, lImplex, lAnnealing
    integer:: nvi, idx_regu_visc, ndimsi
    real(kind=8):: sigm(nsig), epsm(neps), deps(neps)
! --------------------------------------------------------------------------------------------------
    integer :: ndt, ndi
    common/tdim/ndt, ndi
!
! --------------------------------------------------------------------------------------------------
    ASSERT(neps*nsig .eq. ndsde .or. (ndsde .eq. 36 .and. neps .le. 9 .and. nsig .le. 6))

    ndt = 2*ndim
    ndi = ndim

! - Get parameters from COMPOR map
    read (compor(DEFO_LDC), '(A16)') defo_ldc
    read (compor(DEFO), '(A16)') defo_comp
    read (compor(REGUVISC), '(A16)') regu_visc
    read (compor(POSTINCR), '(A16)') postIncr
    l_large = (defo_comp .eq. 'SIMO_MIEHE' .or. defo_comp .eq. 'GROT_GDEP' .or. (compor(RELA_NAME) .eq. 'MFRONT' .and. defo_comp .eq. 'GREEN_LAGRANGE'))
    l_gdef_log = defo_comp .eq. 'GDEF_LOG'
    l_defo_meca = defo_ldc .eq. 'MECANIQUE'
    l_regu_visc = regu_visc .eq. 'REGU_VISC_ELAS'
    lAnnealing = postIncr .eq. "REST_ECRO"

! - Get parameters from modelization
    l_czm = typmod(2) .eq. 'ELEMJOIN'
    l_grad_vari = typmod(2) .eq. 'GRADVARI'

! - Get parameters from option
    lImplex = option .eq. "RIGI_MECA_IMPLEX" .or. option .eq. "RAPH_MECA_IMPLEX"
    l_pred = L_PRED(option)
    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)

! - Initializations
    codret = 0
    if (lSigm) sigp = 0.d0
    if (lMatr) dsidep = 0.d0

! - Prepare parameters at Gauss point
    call behaviourPrepESVAGauss(carcri, defo_ldc, imate, &
                                fami, kpg, ksp, &
                                neps, instap, BEHinteg)

! - Prepare input strains for the behaviour law
    epsm = epsm_tot
    deps = deps_tot
    call behaviourPrepStrain(l_czm, l_large, l_defo_meca, &
                             l_grad_vari, imate, fami, kpg, ksp, &
                             neps, BEHinteg%esva, epsm, deps)

! - Prepare external state variables for external solvers (UMAT/MFRONT)
    if (BEHinteg%l_mfront .or. BEHinteg%l_umat) then
        call behaviourPrepESVAExte(compor, fami, kpg, ksp, BEHinteg)
    end if

! - Prepare index of behaviour law
    if (typmod(2) .eq. 'GDVARINO') then
        numlc = numlc+3000
    end if
    if (typmod(2) .eq. 'GRADSIGM') then
        numlc = numlc+4000
    end if
    if (typmod(2) .eq. 'GRADVARI') then
        numlc = numlc+6000
    end if
    if (typmod(2) .eq. 'EJ_HYME' .or. &
        typmod(2) .eq. 'ELEMJOIN' .or. typmod(2) .eq. 'INTERFAC') then
        numlc = numlc+7000
    end if
    if (lImplex) then
        numlc = numlc+2000
    end if

! - How many internal variables really for the constitutive law
    nvi = nvi_all
    if (l_gdef_log) then
        nvi = nvi-nvi_gdef_log
    end if
    if (l_regu_visc) then
        nvi = nvi-nvi_regu_visc
        idx_regu_visc = nvi+1
    end if
    ASSERT(nvi .ge. 1)

! - What is the stress at t- for the constitutive law
    sigm(1:nsig) = sigm_all(1:nsig)
    if (l_regu_visc) then
        ASSERT(nsig .ge. 2*ndim)
        sigm(1:2*ndim) = sigm(1:2*ndim)-vim(idx_regu_visc:idx_regu_visc-1+2*ndim)*r2(1:2*ndim)
    end if

! - Copy internal states variables
    if (lVari .and. lAnnealing) then
        vip(nvi_all) = vim(nvi_all)
    end if

! --------------------------------------------------------------------------------------------------
    select case (numlc)
    case (1)
!     ELAS
        call lc0001(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (2)
!     VMIS_ISOT_XXX, VISC_ISOT_XXX
        call lc0002(fami, kpg, ksp, ndim, imate, l_epsi_varc, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, sigp, vip, typmod, ndsde, &
                    dsidep, codret)
    case (3)
!     VMIS_CINE_LINE, VMIS_ECMI_XXXX
        call lc0003(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (4)
!     VMIS_CINX_CHAB/MEMO VISC_CINX_CHAB/MEMO,
        call lc0004(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (7)
!     ENDO_ORTH_BETON
        call lc0007(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (8)
!     MAZARS
        call lc0008(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (9)
!     BETON_REGLE_PR
        call lc0009(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (15)
! ----- KIT_META
        call lc0015(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)

    case (16)
!     DRUCK_PRAGER
        call lc0016(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (17)
!     NORTON_HOFF
        call lc0017(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (18)
!     VISC_TAHERI
        call lc0018(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (19)
!     ELAS_HYPER
        call lc0019(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (21)
        call lc0021(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (22)
        call lc0022(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (23)
        call lc0023(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (24)
        call lc0024(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (25)
!     KIT_DDI : NE PAS UTILISER COMME EXEMPLE
        call lc0025(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, &
                    epsm, deps, sigm, vim, option, &
                    sigp, vip, typmod, icomp, &
                    nvi, numlc, dsidep, codret)
    case (26)
        call lc0026(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, materi, &
                    nvi, dsidep, codret)
    case (28)
        call lc0028(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (29)
        call lc0029(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (30)
        call lc0030(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, &
                    typmod, icomp, nvi, dsidep, &
                    codret)
    case (31)
        call lc0031(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, sigm, vim, option, &
                    angmas, sigp, vip, typmod, &
                    icomp, nvi, dsidep, codret)
    case (32)
        call lc0032(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, sigm, vim, option, &
                    angmas, sigp, vip, typmod, icomp, nvi, &
                    dsidep, codret)
    case (33)
        call lc0033(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, &
                    typmod, icomp, nvi, dsidep, &
                    codret)
    case (34)
        call lc0034(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (35)
        call lc0035(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (36)
!     ENDO_ISOT_BETON
        call lc0036(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, angmas, sigp, vip, &
                    typmod, icomp, nvi, ndsde, &
                    dsidep, codret)
    case (40)
!       DRUCKER_PRAGER_NA
        call lc0040(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, angmas, sigp, vip, &
                    typmod, icomp, ndsde, &
                    dsidep, codret)
    case (42)
        call lc0042(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (50)
!     UMAT
        call lc0050(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    imate, compor, carcri, instam, instap, &
                    neps, epsm, deps, nsig, sigm, &
                    nvi, vim, option, angmas, &
                    sigp, vip, dsidep, codret)
    case (54)
        call lc0054(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (55)
        call lc0055(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (58)
!     MFRONT
        call lc0058(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    imate, compor, carcri, instam, instap, &
                    neps, epsm, deps, nsig, sigm, &
                    nvi, vim, option, angmas, &
                    sigp, vip, ndsde, dsidep, codret)
    case (59)
        call lc0059(BEHinteg, &
                    fami, kpg, ksp, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, &
                    typmod, icomp, dsidep, codret)
    case (60)
        call lc0060(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (62)
        call lc0062(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (75)
        call lc0075(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (76)
        call lc0076(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (77)
        call lc0077(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (78)
        call lc0078(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (79)
        call lc0079(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (120)
!     BETON_DOUBLE_DP
        call lc0120(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, l_epsi_varc, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)

    case (137)
!     MONOCRISTAL, POLYCRISTAL
        call lc0137(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, mult_comp, carcri, instam, instap, neps, &
                    epsm, deps, sigm, vim, option, &
                    angmas, sigp, vip, &
                    typmod, icomp, nvi, &
                    dsidep, codret)

    case (145)
!       BETON_RAG : nouvelle
        call lc0145(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)

    case (152)
!     CABLE_GAINE
        call lc0152(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (165)
!     FLUA_PORO_BETON
        call lc0165(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
    case (166)
!     ENDO_PORO_BETON
        call lc0166(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
    case (167)
!     FLUA_ENDO_PORO
        call lc0167(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
    case (168)
!     RGI_BETON
        call lc0168(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
    case (169)
!     RGI_BETON_BA
        call lc0169(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With SIMO_MIEHE
! --------------------------------------------------------------------------------------------------
!
    case (1002)
!     VMIS_ISOT_XXX, VISC_ISOT_XXX
        call lc1002(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, sigp, vip, typmod, ndsde, &
                    dsidep, codret)

    case (1015)
! ----- KIT_META
        call lc1015(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)

    case (1037)
!     ROUSSELIER
        call lc1037(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)

    case (1137)
!     MONOCRISTAL, POLYCRISTAL
        call lc1137(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, mult_comp, carcri, instam, instap, neps, &
                    epsm, deps, sigm, vim, option, &
                    angmas, sigp, vip, &
                    typmod, icomp, nvi, &
                    dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With IMPLEX
! --------------------------------------------------------------------------------------------------
!
    case (2001)
!     ELAS
        call lc2001(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    neps, deps, nsig, sigm, option, &
                    angmas, sigp, vip, typmod, ndsde, &
                    dsidep, codret)

    case (2002)
!     VMIS_ISOT_XXX, VISC_ISOT_XXX
        call lc2002(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, sigp, vip, typmod, ndsde, &
                    dsidep, codret)

    case (2036)
!     ENDO_ISOT_BETON
        call lc2036(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, angmas, sigp, vip, &
                    typmod, icomp, nvi, ndsde, &
                    dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With GDVARINO
! --------------------------------------------------------------------------------------------------
!
    case (3053)
!     ENDO_CARRE
        call lc3053(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)

!
! --------------------------------------------------------------------------------------------------
! - With GRADSIGM
! --------------------------------------------------------------------------------------------------
!
    case (4047)
!     ENDO_HETEROGENE
        call lc4047(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, &
                    icomp, nvi, dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With GRADVARI
! --------------------------------------------------------------------------------------------------
!
    case (6036)
!     ENDO_ISOT_BETON
        call lc6036(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, angmas, sigp, vip, &
                    typmod, icomp, nvi, ndsde, &
                    dsidep, codret)
    case (6046)
!     ENDO_SCALAIRE
        call lc6046(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, angmas, sigp, vip, &
                    typmod, icomp, ndsde, &
                    dsidep, codret)
!
    case (6057)
!     ENDO_FISS_EXP
        call lc6057(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
!
    case (6075)
!     GTN
        call lc6075(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
!
    case (6076)
!     VMIS_ISOT_NL
        call lc6076(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With INTERFAC/EJ_HYME/ELEMJOIN
! --------------------------------------------------------------------------------------------------
!
    case (7010)
!     CZM_EXP_REG
        call lc7010(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (7011)
!     CZM_LIN_REG
        call lc7011(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (7013)
!     JOINT_BA
        call lc7013(fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (7040)
        call lc7040(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (7041)
        call lc7041(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (7043)
        call lc7043(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (7045)
        call lc7045(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (7046)
        call lc7046(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)

    case (7047)
!     JOINT_MECA_ENDO
        call lc7047(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    instam, instap, epsm, &
                    deps, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (7048)
        call lc7048(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    nvi, dsidep, codret)
    case (7049)
        call lc7049(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (7051)
        call lc7051(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (7056)
        call lc7056(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
    case (7058)
!     MFRONT
        call lc7058(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    imate, compor, carcri, instam, instap, &
                    neps, epsm, deps, nsig, sigm, &
                    nvi, vim, option, angmas, &
                    sigp, vip, dsidep, codret)
    case (7077)
        call lc7077(BEHinteg, fami, kpg, ksp, ndim, imate, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, angmas, &
                    sigp, vip, typmod, icomp, &
                    ndsde, dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - For KIT_DDI
! --------------------------------------------------------------------------------------------------
!
    case (8028)
        call lc8028(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, mult_comp, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, angmas, sigp, nvi, vip, &
                    typmod, icomp, ndsde, dsidep, codret)
!
    case (8029)
        call lc8029(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, mult_comp, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, angmas, sigp, nvi, vip, &
                    typmod, icomp, ndsde, dsidep, codret)
!
    case (8057)
        call lc8057(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, mult_comp, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, angmas, sigp, nvi, vip, &
                    typmod, icomp, ndsde, dsidep, codret)
!
    case (8146)
        call lc8146(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, mult_comp, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, angmas, sigp, nvi, vip, &
                    typmod, icomp, ndsde, dsidep, codret)
!
    case (8331)
        call lc8331(BEHinteg, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, mult_comp, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, vim, &
                    option, angmas, sigp, nvi, vip, &
                    typmod, icomp, ndsde, dsidep, codret)
    case default
        call utmess('F', 'COMPOR1_43', si=numlc)
    end select
! --------------------------------------------------------------------------------------------------

! - Prediction
    if (l_Pred .and. lSigm .and. .not. l_defo_meca) then
        sigp = sigm
    end if

! - Viscous regularisation
    if (l_regu_visc .and. codret .ne. 1) then
        ndimsi = 2*ndim
        ASSERT(.not. l_large)
        ASSERT(typmod(2) .eq. ' ' .or. typmod(2) .eq. 'GRADVARI')
        ASSERT(neps .ge. ndimsi)
        ASSERT(nsig .ge. ndimsi)
        call lcvisc(fami, kpg, ksp, ndim, imate, lSigm, lMatr, lVari, &
                    instam, instap, deps(1:ndimsi), &
                    vim(idx_regu_visc:idx_regu_visc+nvi_regu_visc-1), &
                    sigp(1:ndimsi), &
                    vip(idx_regu_visc:idx_regu_visc+nvi_regu_visc-1), &
                    dsidep(1:ndimsi, 1:ndimsi))
    end if
!
    ASSERT(codret .ge. 0 .and. codret .le. 2)
end subroutine
