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
subroutine nmhuj(BEHInteg, &
                 fami, kpg, ksp, typmod, imat, &
                 carcri, epsd, &
                 deps, sigd, vind, option, sigf, &
                 vinf, dsde, iret)
!
    use Behaviour_type
    use MaterialPara_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/get_varc.h"
#include "asterfort/hujcdc.h"
#include "asterfort/hujcic.h"
#include "asterfort/hujcrd.h"
#include "asterfort/hujcri.h"
#include "asterfort/hujdp.h"
#include "asterfort/hujmat.h"
#include "asterfort/hujori.h"
#include "asterfort/hujpre.h"
#include "asterfort/hujprj.h"
#include "asterfort/hujres.h"
#include "asterfort/hujtel.h"
#include "asterfort/hujtid.h"
#include "asterfort/MaterialPara_type.h"
#include "asterfort/mgauss.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
!  INTEGRATION DE LA LOI DE COMPORTEMENT ELASTO PLASTIQUE DE HUJEUX
!  AVEC    . 50 VARIABLES INTERNES
!          . 4 FONCTIONS SEUIL ELASTIQUE DEDOUBLEES AVEC CYCLIQUE
!
!  INTEGRATION DES CONTRAINTES           = SIG(T+DT)
!  INTEGRATION DES VARIABLES INTERNES    = VIN(T+DT)
!  ET CALCUL DU JACOBIEN ASSOCIE         = DS/DE(T+DT) OU DS/DE(T)
!
! --------------------------------------------------------------------------------------------------
!
!  IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!          KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!          TYPMOD  TYPE DE MODELISATION
!          IMAT    ADRESSE DU MATERIAU CODE
!          CRIT    CRITERES  LOCAUX
!                  CRIT(1) = NOMBRE D ITERATIONS MAXI A CONVERGENCE
!                            (ITER_INTE_MAXI == ITECREL)
!                  CRIT(2) = TYPE DE JACOBIEN A T+DT
!                            (TYPE_MATR_COMP == MACOMP)
!                            0 = EN VITESSE     > SYMETRIQUE
!                            1 = EN INCREMENTAL > NON-SYMETRIQUE
!                  CRIT(3) = VALEUR DE LA TOLERANCE DE CONVERGENCE
!                            (RESI_INTE == RESCREL)
!                  CRIT(5) = NOMBRE D'INCREMENTS POUR LE
!                            REDECOUPAGE LOCAL DU PAS DE TEMPS
!                            (RESI_INTE_PAS == ITEDEC )
!                            0 = PAS DE REDECOUPAGE
!                            N = NOMBRE DE PALIERS
!          ANGMAS  LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM),
!                    + UN REEL QUI VAUT 0 SI NAUTIQUIES OU 2 SI EULER
!                    + LES 3 ANGLES D'EULER
!          EPSD    DEFORMATION TOTALE A T
!          DEPS    INCREMENT DE DEFORMATION TOTALE
!          SIGD    CONTRAINTE A T
!          VIND    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!          OPT     OPTION DE CALCUL A FAIRE
!                          'RIGI_MECA_TANG'> DSDE(T)
!                          'FULL_MECA'     > DSDE(T+DT) , SIG(T+DT)
!                          'RAPH_MECA'     > SIG(T+DT)
!  OUT     SIGF    CONTRAINTE A T+DT
!          VINF    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!          DSDE    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!          IRET    CODE RETOUR DE  L'INTEGRATION DE LA LOI CJS
!                         IRET=0 => PAS DE PROBLEME
!                         IRET=1 => ECHEC
!  ----------------------------------------------------------------
!  INFO    MATERD        (*,1) = CARACTERISTIQUES ELASTIQUES A T
!                        (*,2) = CARACTERISTIQUES PLASTIQUES A T
!          MATERF        (*,1) = CARACTERISTIQUES ELASTIQUES A T+DT
!                        (*,2) = CARACTERISTIQUES PLASTIQUES A T+DT
!          NDT             NB DE COMPOSANTE TOTALES DES TENSEURS
!                                  = 6  3D
!                                  = 4  AXIS  C_PLAN  D_PLAN
!                                  = 1  1D
!          NDI             NB DE COMPOSANTE DIRECTES DES TENSEURS
!          NVI             NB DE VARIABLES INTERNES
!  ----------------------------------------------------------------
!  ROUTINE LC....UTILITAIRES POUR INTEGRATION LOI DE COMPORTEMENT
!  ----------------------------------------------------------------
!  ORDRE DES TENSEURS      3D      XX YY ZZ XY XZ YZ
!                          DP      XX YY ZZ XY
!                          AX      RR ZZ TT RZ
!                          1D      XX YY ZZ
!  ----------------------------------------------------------------
!  ATTENTION
!  SI OPT = 'RIGI_MECA_TANG' NE PAS TOUCHER AUX VARIABLES SIGF,VINF
!  QUI N ONT PAS DE PLACE MEMOIRE ALLOUEE
!
!  SIG EPS DEPS  ONT DEJA LEURS COMPOSANTES DE CISAILLEMENT
!  MULTIPLIES PAR RACINE DE 2 > PRISE EN COMPTE DES DOUBLES
!  PRODUITS TENSORIELS ET CONSERVATION DE LA SYMETRIE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)      :: imat, ndt, ndi, nvi, iret, iret1, kpg, ksp
    integer(kind=8)      :: i, inc, incmax, ndtt, limsup
    real(kind=8) :: carcri(CARCRI_SIZE), vind(50), vinf(50), vind0(50), variTmp(50)
    real(kind=8) :: epsd(6), deps(6), deps0(6)
    real(kind=8) :: sigd(6), sigf(6), dsde(6, 6), seuil
    real(kind=8) :: piso, depsr(6), depsq(6), tin(3)
    real(kind=8) :: d, q, m, phi, b
    real(kind=8) :: pc0, sigd0(6), hill, dsig(6)
    character(len=7)  :: etatd, etatf
    character(len=8)  :: typmod1, typmod(2)
    character(len=16) :: option
    character(len=*)  :: fami
    real(kind=8) :: depsth(6), alpha(3), tempm, tempf, tref
    real(kind=8) :: det, bid16(6), bid66(6, 6)
    real(kind=8) :: materf(22, 2)
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0, deux = 2.d0, trois = 3.d0
    real(kind=8), parameter :: degr = 0.0174532925199d0, tole = 0.1d0
    real(kind=8) :: neps, nsig, ptrac, rtrac
    real(kind=8) :: crit, dpiso
    aster_logical:: debug, conv, reorie, tract, lVari
    common/tdim/ndt, ndi
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = BEHInteg%materPara
!
    iret = 0
! --- DEBUG = ASTER_TRUE : MODE AFFICHAGE ENRICHI
    debug = ASTER_FALSE
    tract = ASTER_FALSE
    conv = ASTER_TRUE
    reorie = ASTER_FALSE

! - Flag to modify internal state variable
    lVari = L_VARI(option)
!
    variTmp = 0.d0
    deps0 = 0.d0
    seuil = 0.d0
    piso = 0.d0
    depsr = 0.d0
    depsq = 0.d0
    tin = 0.d0
    d = 0.d0
    q = 0.d0
    m = 0.d0
    phi = 0.d0
    b = 0.d0
    pc0 = 0.d0
    sigd0 = 0.d0
    hill = 0.d0
    dsig = 0.d0
    depsth = 0.d0
    alpha = 0.d0
    tempm = 0.d0
    tempf = 0.d0
    tref = 0.d0
    det = 0.d0
    bid16 = 0.d0
    bid66 = 0.d0
    materf = 0.d0
    neps = 0.d0
    nsig = 0.d0
    ptrac = 0.d0
    rtrac = 0.d0
    crit = 0.d0
    dpiso = 0.d0
    dsde = 0.d0
!
    if (debug) then
        write (6, *)
        write (6, '(A)') '!!!!(@_@)!!!!                        !!!!(@_@)!!!!'
        write (6, '(A)') '!!!!(@_@)!!!!    MODE DEBUG ACTIF    !!!!(@_@)!!!!'
        write (6, '(A)') '!!!!(@_@)!!!!                        !!!!(@_@)!!!!'
    end if
!
    typmod1 = typmod(1)
!
! - Get temperatures
!
    call get_varc(fami, kpg, ksp, 'T', tempm, tempf, tref)
!
! ---> RECUPERATION COEF DE LA LOI HUJEUX
!      (INDEPENDANTS DE LA TEMPERATURE)
!      NB DE CMP DIRECTES/CISAILLEMENT
!      NB VARIABLES INTERNES
    call hujmat(fami, kpg, ksp, typmod1, imat, &
                tempf, materf, ndt, ndi, nvi)
!
    ptrac = materf(21, 2)
    rtrac = abs(1.d-6*materf(8, 2))
!
! --- REORIENTATION DES PLANS DE GLISSEMENT SUR LES AXES DU
!     REPERE LOCAL DONNE PAR LES ANGLES NAUTIQUES (ANGMAS)
    if (materPara%lcsPara%lcsType .eq. MATER_LCS_NONE) then
        call utmess('F', 'ALGORITH8_20')
    end if
!
    reorie = (materPara%lcsPara%lcsAngle(1) .ne. zero) .or. &
             (materPara%lcsPara%lcsAngle(2) .ne. zero) .or. &
             (materPara%lcsPara%lcsAngle(3) .ne. zero)
!
    call hujori('LOCAL', 1, reorie, materPara%lcsPara%lcsAngle, sigd, bid66)
    call hujori('LOCAL', 1, reorie, materPara%lcsPara%lcsAngle, epsd, bid66)
    call hujori('LOCAL', 1, reorie, materPara%lcsPara%lcsAngle, deps, bid66)
!
! --- ON TRAVAILLE TOUJOURS AVEC UN TENSEUR CONTRAINTES
!     DEFINI EN 3D
!
    ndtt = 6
    if (ndt .lt. 6) then
        ndtt = 4
        ndt = 6
    end if
!
!     CALCUL DE DEPSTH ET EPSDTH
!     --------------------------
! ---> COEF DE DILATATION LE MEME A TPLUS ET TMOINS
    if (materf(17, 1) .eq. un) then
!
        if ((isnan(tempm) .or. isnan(tref)) .and. &
            materf(3, 1) .ne. zero) then
            call utmess('F', 'CALCULEL_15')
        end if
!
        alpha(1) = materf(3, 1)
        alpha(2) = materf(3, 1)
        alpha(3) = materf(3, 1)
!
    else if (materf(17, 1) .eq. deux) then
!
        alpha(1) = materf(10, 1)
        alpha(2) = materf(11, 1)
        alpha(3) = materf(12, 1)
        if ((isnan(tempm) .or. isnan(tref)) .and. &
            (alpha(1) .ne. zero .or. alpha(2) .ne. zero .or. &
             alpha(3) .ne. zero)) then
            call utmess('F', 'CALCULEL_15')
        end if
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (isnan(tempm) .or. isnan(tempf) .or. isnan(tref)) then
        do i = 1, ndi
            depsth(i) = deps(i)
        end do
    else
        do i = 1, ndi
            depsth(i) = deps(i)- &
                        alpha(i)*(tempf-tref)+alpha(i)*(tempm-tref)
        end do
    end if
!
    do i = ndi+1, ndt
        depsth(i) = deps(i)
    end do
!
    if (ndtt .lt. 6) then
        do i = ndtt+1, 6
            depsth(i) = zero
            sigd(i) = zero
        end do
    end if
!
! ---> INITIALISATION SEUIL DEVIATOIRE SI NUL
!
    do i = 1, ndi
        if (vind(i) .eq. zero) then
!
            if (materf(13, 2) .eq. zero) then
                vind(i) = 1.d-3
            else
                vind(i) = materf(13, 2)
            end if
!
            call hujcrd(i, materf, sigd, vind, seuil, iret)
            if (iret .ne. 0) then
                goto 999
            end if
!
! --- SI LE SEUIL EST DESEQUILIBRE A L'ETAT INITIAL
!     ON EQUILIBRE LE SEUIL EN CALCULANT LA VALEUR DE R
!     APPROPRIEE
!
            if (seuil .gt. zero) then
                call hujprj(i, sigd, tin, piso, q)
                piso = piso-ptrac
                b = materf(4, 2)
                phi = materf(5, 2)
                m = sin(degr*phi)
                pc0 = materf(7, 2)
                vind(i) = -q/(m*piso*(un-b*log(piso/pc0)))
                vind(23+i) = un
            end if
!
        end if
    end do
!
! ---> INITIALISATION SEUIL ISOTROPE SI NUL
    if (vind(4) .eq. zero) then
        if (materf(14, 2) .eq. zero) then
            vind(4) = 1.d-3
        else
            vind(4) = materf(14, 2)
        end if
!
        call hujcri(materf, sigd, vind, seuil)
!
! --- SI LE SEUIL EST DESEQUILIBRE A L'ETAT INITIAL
!     ON EQUILIBRE LE SEUIL EN CALCULANT LA VALEUR DE R
!     APPROPRIEE
!
        if (seuil .gt. zero) then
            piso = trace(3, sigd)/trois
            d = materf(3, 2)
            pc0 = materf(7, 2)
            vind(4) = piso/(d*pc0)
            if (vind(4) .gt. 1.d0) then
                call utmess('F', 'COMPOR1_83')
            end if
            vind(27) = un
        end if
!
    end if
!
! ---> INITIALISATION SEUIL CYCLIQUE SI NUL
    do i = 1, ndi
        if (vind(4+i) .eq. zero) then
            if (materf(18, 2) .eq. zero) then
                vind(4+i) = 1.d-3
            else
                vind(4+i) = materf(18, 2)
            end if
        end if
    end do
!
    if (vind(8) .eq. zero) then
        if (materf(19, 2) .eq. zero) then
            vind(8) = 1.d-3
        else
            vind(8) = materf(19, 2)
        end if
    end if
!
!ONTROLE DES INDICATEURS DE PLASTICITE
    do i = 1, 4
        if (abs(vind(27+i)-un) .lt. r8prem()) vind(23+i) = -un
    end do
!
    if (option(1:9) .ne. 'RIGI_MECA') then
        variTmp(1:50) = vind(1:50)
    end if
!
! ---> ETAT ELASTIQUE OU PLASTIQUE A T
    if (((vind(24) .eq. zero) .or. &
         (vind(24) .eq. -un .and. vind(28) .eq. zero)) .and. &
        ((vind(25) .eq. zero) .or. &
         (vind(25) .eq. -un .and. vind(29) .eq. zero)) .and. &
        ((vind(26) .eq. zero) .or. &
         (vind(26) .eq. -un .and. vind(30) .eq. zero)) .and. &
        ((vind(27) .eq. zero) .or. &
         (vind(27) .eq. -un .and. vind(31) .eq. zero))) then
        etatd = 'ELASTIC'
    else
        etatd = 'PLASTIC'
    end if
!
! -------------------------------------------------------------
! OPTIONS 'FULL_MECA' ET 'RAPH_MECA' = CALCUL DE SIG(T+DT)
! -------------------------------------------------------------
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
!
        if (debug) write (6, *) ' * DEPS =', (depsth(i), i=1, 3)
!
        do i = 1, 3
            call hujprj(i, sigd, tin, piso, q)
            if (abs(piso+deux*rtrac-ptrac) .lt. r8prem()) &
                tract = ASTER_TRUE
        end do
!
! INTEGRATION ELASTIQUE SUR DT
        do i = 1, ndt
            depsq(i) = zero
        end do
!
! -----------------------------------------------
! ---> INCREMENT TOTAL DE DEFORMATION A APPLIQUER
! -----------------------------------------------
! ENREGISTREMENT DE L'ETAT DE CONTRAINTES A T
        sigd0(1:ndt) = sigd(1:ndt)
!
! ENREGISTREMENT DE L'INCREMENT TOTAL DEPS0
        deps0(1:ndt) = depsth(1:ndt)
!
! INITIALISATION DES DEFORMATIONS RESTANTES
        depsq(1:ndt) = depsth(1:ndt)
        vind0(1:nvi) = vind(1:nvi)
!
! INITIALISATION DU COMPTEUR D'ITERATIONS LOCALES
!        vind(35) = zero
!
! -----------------------------------------------------
! ---> PREDICTION VIA TENSEUR ELASTIQUE DES CONTRAINTES
! -----------------------------------------------------
        inc = 0
        incmax = 1
        limsup = int(max(20.d0, abs(carcri(1))))
!
100     continue
!
        inc = inc+1
        depsr(1:ndt) = depsq(1:ndt)
        call hujpre(fami, kpg, ksp, etatd, typmod1, &
                    imat, materf, depsr, sigd, &
                    sigf, vind0, iret)
!
        if (iret .eq. 1) then
            if (debug) &
                write (6, '(A)') &
                '!!!@_@!!! NMHUJ :: ARRET DANS HUJPRE !!!@_@!!!'
            goto 999
        end if
!
! ----------------------------------------------------
! ---> CONTROLE DE L EVOLUTION DE LA PRESSION ISOTROPE
! ----------------------------------------------------
        iret1 = 0
        call hujdp(typmod1, depsr, sigd, sigf, materf, &
                   vind, incmax, iret1)

        if (iret1 .eq. 1) then
            if (debug) &
                write (6, '(A)') &
                '!!!@_@!!! NMHUJ :: ARRET DANS HUJDP :: PAS DE RESUBDIVISION !!!@_@!!!'
!           goto 999
        end if
!
! --- ON LIMITE LE REDECOUPAGE LOCAL A MAX(20,ITER_INTE_MAXI)
        if (incmax .ge. limsup) then
            incmax = limsup
        else if (incmax .le. 1) then
            incmax = 1
        end if
!
        if (inc .eq. 1 .and. incmax .gt. 1) then
            do i = 1, ndt
                depsq(i) = deps0(i)/incmax
                depsr(i) = deps0(i)/incmax
            end do
            call hujpre(fami, kpg, ksp, etatd, typmod1, &
                        imat, materf, depsr, sigd, &
                        sigf, vind0, iret)
        end if
!
! ---------------------------------------------
! CALCUL DE L'ETAT DE CONTRAINTES CORRESPONDANT
! ---------------------------------------------
        if (debug) write (6, *) &
            '!!!@_@!!! NMHUJ -- VINF =', (variTmp(i), i=24, 31), ' !!!@_@!!!'
!
        call hujres(fami, kpg, ksp, typmod1, carcri, &
                    materf, imat, nvi, depsr, sigd, &
                    vind, sigf, variTmp, iret, etatf)
        if (iret .eq. 1) goto 999
!
! -------------------------------------------
! - CONTROLE DES DEFORMATIONS DEJA APPLIQUEES
! -------------------------------------------
        if (inc .lt. incmax) then
            conv = ASTER_FALSE
        else
            conv = ASTER_TRUE
        end if
!
        if (.not. conv) then
            sigd(1:ndt) = sigf(1:ndt)
            vind(1:nvi) = variTmp(1:nvi)
            goto 100
        end if
!
! --- CALCUL DU CRITERE DE HILL: DSIG*DEPS
        hill = zero
        nsig = zero
        neps = zero
        do i = 1, ndt
            dsig(i) = sigf(i)-sigd0(i)
            hill = hill+dsig(i)*deps0(i)
            nsig = nsig+dsig(i)**2
            neps = neps+deps0(i)**2
        end do
!
! --- NORMALISATION DU CRITERE : VARIE ENTRE -1 ET 1
        if (neps .gt. r8prem() .and. nsig .gt. r8prem()) then
            variTmp(32) = hill/sqrt(neps*nsig)
        else
            variTmp(32) = zero
        end if
!
    end if
!af 07/05/07 fin <IF RAPH_MECA et FULL_MECA>
!
!       ----------------------------------------------------------------
!       OPTIONS 'FULL_MECA' ET 'RIGI_MECA_TANG' = CALCUL DE DSDE
!       ----------------------------------------------------------------
!       CALCUL ELASTIQUE ET EVALUATION DE DSDE A (T)
!       POUR 'RIGI_MECA_TANG' ET POUR 'FULL_MECA'
!       ----------------------------------------------------------------
    if (option .eq. 'RIGI_MECA_TANG') then
!
        dsde(:, :) = zero
!
! REMARQUE: CALCUL DE DSDE A T AVEC MATERF CAR PARAMETRES HUJEUX
! --------  INDEPENDANTS DE LA TEMPERATURE
!
! ---> CALCUL MATRICE DE RIGIDITE ELASTIQUE
        if (etatd .eq. 'ELASTIC') then
            call hujtel(typmod1, materf, sigd, dsde)
        end if
!
! ---> CALCUL MATRICE TANGENTE DU PROBLEME CONTINU
        if (etatd .eq. 'PLASTIC') then
            call hujtid(fami, kpg, ksp, typmod1, imat, &
                        sigd, vind, dsde, iret)
            if (iret .eq. 1) goto 999
        end if
!
        call hujori('GLOBA', 2, reorie, materPara%lcsPara%lcsAngle, bid16, dsde)
!
    else if (option .eq. 'FULL_MECA') then
!
        dsde(:, :) = zero
!
! ---> CALCUL MATRICE DE RIGIDITE ELASTIQUE
        if (etatf .eq. 'ELASTIC') then
            call hujtel(typmod1, materf, sigf, dsde)
        end if
!
! ---> CALCUL MATRICE TANGENTE DU PROBLEME CONTINU
        if (etatf .eq. 'PLASTIC') then
            call hujtid(fami, kpg, ksp, typmod1, imat, &
                        sigf, variTmp, dsde, iret)
            if (iret .eq. 1) goto 999
        end if
!
    else if (option .eq. 'FULL_MECA_ELAS') then
!
        dsde(:, :) = zero
        call hujtel(typmod1, materf, sigf, dsde)
!
    else if (option .eq. 'RIGI_MECA_ELAS') then
!
        dsde(:, :) = zero
        call hujtel(typmod1, materf, sigd, dsde)
        call hujori('GLOBA', 2, reorie, materPara%lcsPara%lcsAngle, bid16, dsde)
!
    end if
! fin <IF RIGI_MECA_TANG>
!
! ---> CALCUL DETERMINANT DE LA MATRICE TANGENTE + INDICATEUR
! --- RELIE AUX MECANISMES ACTIFS
    if (option(1:9) .ne. 'RIGI_MECA') then
!
        call hujori('GLOBA', 2, reorie, materPara%lcsPara%lcsAngle, bid16, &
                    dsde)
!
        if (option .eq. 'FULL_MECA') then
            call mgauss('NCSD', dsde, sigd, 6, 6, &
                        1, det, iret)
            if (iret .eq. 1) then
                variTmp(33) = un
                iret = 0
            else
                variTmp(33) = det
            end if
        end if
!
    end if
! --- ON RENVOIE LA VALEUR ADEQUATE DE NDT
!     POUR MODELISATION D_PLAN
    if (ndtt .eq. 4) ndt = 4
!
    if (option .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') &
        call hujori('GLOBA', 1, reorie, materPara%lcsPara%lcsAngle, sigf, bid66)
!
999 continue
!
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA' &
        .or. option(1:14) .eq. 'RIGI_MECA_ELAS') then
        if (iret .eq. 1) then
            if (.not. tract) then
                dsde(:, :) = zero
                call hujtid(fami, kpg, ksp, typmod1, imat, &
                            sigd, vind, dsde, iret1)
                if (iret1 .eq. 1) then
                    dsde(:, :) = zero
                    call hujtel(typmod1, materf, sigd, dsde)
                end if
! debut ---new dvp 23/01/2019---
                dsig(1:ndt) = matmul(dsde(1:ndt, 1:ndt), deps0(1:ndt))
!
! on limite la variation de dsig a 2% par rapport a piso
                piso = trace(3, sigd0)
                dpiso = trace(3, dsig)
                if (piso .lt. -100. .and. abs(dpiso/piso) .le. tole) then
                    det = 1.
                elseif (piso .lt. -100.) then
                    det = tole*abs(piso/dpiso)
                else
                    det = 0.
                end if
                dsig(1:ndt) = det*dsig(1:ndt)

                sigf = sigd0+dsig
!
! y-a-t-il traction?
                conv = ASTER_TRUE
                do i = 1, 3
                    call hujprj(i, sigf, tin, piso, q)
                    if (abs(piso+deux*rtrac-ptrac) .lt. r8prem()) &
                        conv = ASTER_FALSE
                end do
!
                if (.not. conv) then
                    do i = 1, 3
                        sigf(i) = -deux*rtrac+ptrac
                        sigf(i+3) = zero
                    end do
                    dsde(:, :) = zero
                    call hujtel(typmod1, materf, sigd, dsde)
                end if
!
                variTmp(1:50) = vind0(1:50)
! fin   ---new dvp 23/01/2019---
                if (debug) then
                    write (6, '(A)') ' ----------- FIN NMHUJ -----------------'
                    write (6, *) ' * DEPS =', (deps0(i), i=1, ndt)
                    write (6, *) ' * SIGD =', (sigd0(i), i=1, ndt)
                    write (6, *) ' * VIND =', (vind0(i), i=1, 50)
                    write (6, *) ' * SIGF =', (sigf(i), i=1, ndt)
                    write (6, *) ' * VINF =', (variTmp(i), i=1, 50)
                    write (6, '(A)') ' ----------------------------------------'
                end if
            else
!
                dsde(:, :) = zero
                call hujtel(typmod1, materf, sigd, dsde)
                do i = 1, 3
                    sigf(i) = -deux*rtrac+ptrac
                    sigf(i+3) = zero
                end do
                variTmp(1:50) = vind0(1:50)
                iret = 0
            end if
        end if
        sigd(1:ndt) = sigd0(1:ndt)
        deps(1:ndt) = deps0(1:ndt)
        vind(1:50) = vind0(1:50)
! debut ---new dvp 23/01/2019---
! sauvegarde d'une estimation de l'erreur cumulee
! L'erreur est mesuree sur F=(sigm, vari_interne)=0
!
        crit = zero
!
        if (conv) then
            det = abs(trace(3, deps0))
            do i = 1, 3
!
! On normalise les seuils de la meme facon que dans hujjid
! i.e. par le module d'Young materf(1,1)/Pcr0
! pour assurer la coherence du controle avec RESI_INTE
                call hujcrd(i, materf, sigf, variTmp, seuil, iret)
                if (iret .ne. 0) then
                    goto 999
                end if
                seuil = seuil*det

                if (seuil .gt. zero) then
                    bid16(i) = un
                else
                    bid16(i) = zero
                end if
                crit = max(seuil, crit)
!
            end do
!
! On normalise les seuils de la meme facon que dans hujjid
! i.e. par le module d'Young materf(1,1)/Pcr0
! pour assurer la coherence du controle avec RESI_INTE
!
            call hujcri(materf, sigf, variTmp, seuil)
            seuil = seuil/materf(1, 1)*abs(materf(7, 2))

            if (seuil .gt. zero) then
                bid16(4) = un
            else
                bid16(4) = zero
            end if
            crit = max(seuil, crit)
!
            do i = 5, 8
!
! On normalise les seuils de la meme facon que dans hujjid
! i.e. par le module d'Young materf(1,1)/Pcr0
! pour assurer la coherence du controle avec RESI_INTE
!
                if (i .lt. 8 .and. bid16(i-4) .eq. zero) then

                    call hujcdc(i-4, materf, sigf, variTmp, seuil)
                    seuil = seuil*det

                elseif (bid16(4) .eq. zero) then

                    call hujcic(materf, sigf, variTmp, seuil)
                    seuil = seuil/materf(1, 1)*abs(materf(7, 2))
                end if
                crit = max(seuil, crit)
!
            end do
        end if
!
        variTmp(34) = crit
        variTmp(35) = zero
! fin   ---new dvp 23/01/2019---
    end if

! - Copy internal state variables
    if (lVari) then
        vinf = variTmp
    end if
!
end subroutine
