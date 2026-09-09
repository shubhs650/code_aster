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
! aslint: disable=W0413
!
subroutine nmvprk(BEHInteg, &
                  option, typmod, ndim, &
                  compor, carcri, &
                  instam, instap, &
                  neps, epsdt, depst, sigd, &
                  nvi, vind, sigf, &
                  vinf, dsde, iret, multComp_)
!
    use Behaviour_type
    use MaterialPara_type
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterfort/calsig.h"
#include "asterfort/gerpas.h"
#include "asterfort/lcdpeq.h"
#include "asterfort/lcmate.h"
#include "asterfort/lcopli.h"
#include "asterfort/lcrkin.h"
#include "asterfort/lcrksg.h"
#include "asterfort/lcsmelas.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(in) :: instam, instap
    integer(kind=8), intent(in) :: nvi
    character(len=16), optional, intent(in) :: multComp_
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour - The RUNGE_KUTTA environment (prefer MFront !)
!
! Main subroutine
!
! --------------------------------------------------------------------------------------------------
!
!     INTEGRATION DE LOIS DE COMPORTEMENT ELASTO-VISCOPLASTIQUE
!     PAR UNE METHODE DE RUNGE KUTTA
!             AVEC    .N VARIABLES INTERNES
!
!     INTEGRATION DES CONTRAINTES           = SIG(T+DT)
!     INTEGRATION DES VARIABLES INTERNES    = VIN(T+DT)
!
!     CETTE METHODE NE FOURNIT PAS DE MATRICE TANGENTE POUR LA
!     RESOLUTION GLOBALE.
!     ON FOURNIT DONC LA MATRICE DSDE OBTENUE EN ELASTICITE
!     ================================================================
!
!     ARGUMENTS
!
!     IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!           KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!           NDIM   DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!           TYPMOD TYPE DE MODELISATION
!           IMAT   ADRESSE DU MATERIAU CODE
!           COMP   COMPORTEMENT DE L ELEMENT
!                  COMP(1) = RELATION DE COMPORTEMENT (CHABOCHE...)
!                  COMP(2) = NB DE VARIABLES INTERNES
!                  COMP(3) = TYPE DE DEFORMATION (PETIT,JAUMANN...)
!           OPT    OPTION DE CALCUL A FAIRE
!                          'RIGI_MECA_TANG'> DSDE(T)
!           CRIT   CRITERES  LOCAUX
!                  CRIT(1) = NOMBRE D ITERATIONS MAXI A CONVERGENCE
!                            (ITER_INTE_MAXI == ITECREL)
!                  CRIT(3) = CRITERE DE PRECISION POUR L INTEGRATION
!                            PAR LA METHODE DE RUNGE KUTTA
!                  CRIT(6) = TYPE D INTEGRATION LOCAL POUR LA LOI DE
!                            COMPORTEMENT (ALGO_INTE)
!           TIMED   INSTANT T
!           TIMEF   INSTANT T+DT
!           EPSDT   DEFORMATION TOTALE A T
!           DEPST   INCREMENT DE DEFORMATION TOTALE
!           SIGD    CONTRAINTE A T
!           VIND    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!           ANGMAS  3 ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
!     OUT   SIGF    CONTRAINTE A T+DT
!           VINF    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!           DSDE    MATRICE DE COMPORTEMENT TANGENT ELASTIQUE
!     ----------------------------------------------------------------
!     INFO  MATERD        (*,1) = CARACTERISTIQUES ELASTIQUES A T
!                         (*,2) = CARACTERISTIQUES PLASTIQUES A T
!           MATERF        (*,1) = CARACTERISTIQUES ELASTIQUES A T+DT
!                         (*,2) = CARACTERISTIQUES PLASTIQUES A T+DT
!           MATCST        'OUI' SI MATERIAU CST ENTRE T ET T+DT
!                         'NAP' SI LE PARAMETRE K_D EST UNE NAPPE
!                         'NON' SINON
!           NDT            NB DE COMPOSANTES TOTALES DES TENSEURS
!                                  = 6  3D
!                                  = 4  AXIS  C_PLAN  D_PLAN
!           NDI            NB DE COMPOSANTES DIRECTES DES TENSEURS
!           NVI            NB DE VARIABLES INTERNES
!           NR             NB EQUATIONS SYSTEME INTEGRE A RESOUDRE
!     ----------------------------------------------------------------
!     ROUTINE LC....UTILITAIRES POUR INTEGRATION LOI DE COMPORTEMENT
!     ----------------------------------------------------------------
!     ORDRE DES TENSEURS  3D XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ
!                         DP XX YY ZZ SQRT(2)*XY
!                         AX RR ZZ TT SQRT(2)*RZ
!     ----------------------------------------------------------------
!     ATTENTION
!     SI OPT = 'RIGI_MECA_TANG' NE PAS TOUCHER AUX VARIABLES SIGF,VINF
!     QUI N ONT PAS DE PLACE MEMOIRE ALLOUEE
!
!     SIG EPS DEPS  ONT DEJA LEURS COMPOSANTES DE CISAILLEMENT
!     MULTIPLIES PAR RACINE DE 2 > PRISE EN COMPTE DES DOUBLES
!     PRODUITS TENSORIELS ET CONSERVATION DE LA SYMETRIE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: rungeKutta = 1, nmat = 6000
    integer(kind=8), parameter ::  nsg = 30, nfs = 5, nhsr = 5
    integer(kind=8) :: ndim, ndt, ndi, nr, i, nbphas, itmax
    integer(kind=8) :: ioptio, idnr, neps
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
!     POUR LCMATE (MONOCRISTAL) DIMENSIONS MAX
!        NSG=NOMBRE DE SYSTEMES DE GLISSEMENT MAXIMUM
!        NFS=NOMBRE DE FAMILLES DE SYSTEMES DE GLISSEMENT MAXIMUM
    integer(kind=8), allocatable :: nbcomm(:, :)
    real(kind=8), allocatable :: materd(:, :), materf(:, :)
    character(len=24), allocatable :: cpmono(:)
    integer(kind=8), allocatable :: numhsr(:)
    real(kind=8), allocatable :: cothe(:), dcothe(:)
    real(kind=8), allocatable :: coeff(:), dcoeff(:), coel(:)
    integer(kind=8) :: iret
    real(kind=8) :: epsdt(neps), depst(neps)
    real(kind=8) :: rbid
    real(kind=8) :: toler, ymfs, vind(*), vinf(*)
    real(kind=8) :: sigd(6), sigf(6), dsde(6, *)
    real(kind=8) :: pgl(3, 3), epsd(9)
    real(kind=8) :: dtime, x
!     POUR POLYCRISTAL, 5 MATRICE HSR MAXI. POUR MONOCRISTAL, 1 MAXI
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg, nhsr), detot(9)
    character(len=3) :: matcst
    character(len=8) :: typmod1, typma
    character(len=11) :: meting
    character(len=16) ::  relaComp, defoComp, multComp
    blas_int :: b_incx, b_incy, b_n
    common/tdim/ndt, ndi
    common/opti/ioptio, idnr
    common/meti/meting
    common/polycr/irr, decirr, nbsyst, decal, gdef
    type(Material_Para) :: materPara
    character(len=8) :: fami
    integer(kind=8) :: jvMaterCode, kpg, ksp
!
! --------------------------------------------------------------------------------------------------
!
    materPara = BEHInteg%materPara
    jvMaterCode = materPara%jvMaterCode
    fami = materPara%schemePara%fami
    kpg = materPara%schemePara%kpg
    ksp = materPara%schemePara%ksp

    itmax = int(carcri(1))
    toler = carcri(3)
    meting = 'RUNGE_KUTTA'
    typmod1 = typmod(1)
    relaComp = compor(RELA_NAME)
    defoComp = compor(DEFO)
    multComp = ' '
    if (present(multComp_)) then
        multComp = multComp_
    end if
    gdef = 0
    if (defoComp .eq. 'SIMO_MIEHE') gdef = 1

! - Allocate big objects
    allocate (nbcomm(nmat, 3))
    allocate (materd(nmat, 2))
    allocate (materf(nmat, 2))
    allocate (cpmono(5*nmat+1))
    allocate (numhsr(nmat))
    allocate (cothe(nmat))
    allocate (dcothe(nmat))
    allocate (coeff(nmat))
    allocate (dcoeff(nmat))
    allocate (coel(nmat))
!
!     YMFS EST UTILISE LORS DU CALCUL D ERREUR COMME MINIMUM DE
!     CHAQUE COMPOSANTE DE VINT. L IDEAL SERAIT DE RENTRER CE
!     PARAMETRE EN DONNEE POUR CHAQUE VARIABLE INTERNE
!
    ymfs = 0.0001d0
!
! --  RECUPERATION COEF(TEMP(T))) LOI ELASTO-PLASTIQUE A T ET/OU T+DT
!                    NB DE CMP DIRECTES/CISAILLEMENT + NB VAR. INTERNES
!
    call lcmate(materPara, &
                carcri, relaComp, typmod1, &
                nmat, rbid, rbid, rbid, rungeKutta, &
                typma, hsr, materd, materf, matcst, &
                nbcomm, cpmono, pgl, 0, &
                toler, ndt, ndi, nr, &
                nvi, vind, nfs, nsg, toutms, &
                nhsr, numhsr, sigd, multComp)
!
    if (option(1:9) .eq. 'RIGI_MECA') goto 900
!
    b_n = to_blas_int(neps)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, depst, b_incx, detot, b_incy)
    b_n = to_blas_int(neps)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsdt, b_incx, epsd, b_incy)
!
    dtime = instap-instam
!
! --  INITIALISATION DES VARIABLES INTERNES A T
!
    do i = 1, nmat
        cothe(i) = materd(i, 1)
        dcothe(i) = -cothe(i)+materf(i, 1)
    end do
!
    do i = 1, nmat
        coeff(i) = materd(i, 2)
        dcoeff(i) = -coeff(i)+materf(i, 2)
    end do
!
!     INITIALISATIONS PARTICULIERES POUR CERTAINES LOIS
!
    call lcrkin(ndim, option, relaComp, materf, nbcomm, &
                cpmono, nmat, typmod1, nvi, sigd, &
                sigf, vind, vinf, nbphas, iret)
    if (iret .eq. 9) then
!        ENDOMMAGEMENT MAXI AU POINT DE GAUSS
        iret = 0
        goto 999
    end if
!
    call gerpas(materPara, &
                relaComp, typmod1, &
                matcst, nbcomm, cpmono, nbphas, &
                nvi, nmat, vinf, dtime, itmax, &
                toler, ymfs, cothe, coeff, dcothe, &
                dcoeff, coel, pgl, neps, &
                epsd, detot, x, nfs, nsg, &
                nhsr, numhsr, hsr, iret)
    if (iret .ne. 0) then
        goto 999
    end if
!
! --  CALCUL DES CONTRAINTES
!
    if ((relaComp .eq. 'MONOCRISTAL') .and. (gdef .eq. 1)) then
        call lcrksg(relaComp, nvi, vinf, epsd, detot, &
                    nmat, coel, sigf)
    else
        call calsig(fami, kpg, ksp, vinf, typmod1, &
                    relaComp, vinf, x, dtime, epsd, &
                    detot, nmat, coel, sigf)
    end if
!
    call lcdpeq(vind, vinf, relaComp, nbcomm, cpmono, &
                nmat, nvi, sigf, detot, epsd, &
                materf, pgl)
!
900 continue
!
    if (option(1:10) .eq. 'RIGI_MECA_' .and. gdef .eq. 1 .and. relaComp .eq. 'MONOCRISTAL') then
        call lcsmelas(epsdt, depst, dsde, nmat=nmat, materd_=materd)
        iret = 0
        goto 999
    end if
!
!     OPERATEUR TANGENT = ELASTIQUE OU SECANT (ENDOMMAGEMENT)
    if (materf(nmat, 1) .eq. 0) then
        call lcopli('ISOTROPE', typmod1, materf(1, 1), dsde)
    else if (materf(nmat, 1) .eq. 1) then
        call lcopli('ORTHOTRO', typmod1, materf(1, 1), dsde)
    end if
!
999 continue

! - De-Allocate big objects
    deallocate (nbcomm)
    deallocate (materd)
    deallocate (materf)
    deallocate (cpmono)
    deallocate (numhsr)
    deallocate (cothe)
    deallocate (dcothe)
    deallocate (coeff)
    deallocate (dcoeff)
    deallocate (coel)
!
end subroutine
