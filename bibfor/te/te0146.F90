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
subroutine te0146(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/clcplq.h"
#include "asterfort/jevech.h"
#include "asterfort/sandwich.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!  BUT: CALCUL DE L'OPTION FERRAILLAGE POUR LES ELEMENTS DE COQUE
!
! CALCUL DES DENSITES DE FERRAILLAGE DANS LE BETON ARME
!              (METHODE DE CAPRA ET MAURY)
!
! --------------------------------------------------------------------------------------------------
!
! PARAMETRES D'ECHANGE ENTRE CODE_ASTER ET CLCPLQ
! (POINT D'ENTREE DU CALCUL DE FERRAILLAGE PAR CAPRA ET MAURY)
!
!   PARAMETRES D'ENTREE (FOURNIS PAR CODE_ASTER)
!
!     TYPCMB     TYPE DE COMBINAISON :
!                   0 = ELU, 1 = ELS, 2 = ELS QP
!     TYPCO      TYPE DE CODIFICATION :
!                   1 = BAEL91, 2 = EUROCODE 2
!     METH2D     CHOIX DE METHODE POUR LE CALCUL DU FERRAILLAGE DES
!                   ELEMENTS 2D
!                   1 = CapraMaury, 2 = Sandwich/Multicouches
!     THITER     ANGLE D'ITERATION DANS LE CAS DU CALCUL DU FERRAILLAGE DES PLAQUES
!                   (POUR LA METHODE SANDWICH, IL S'AGIT DE L'INCLINAISON DES BIELLES)
!                   (POUR LA METHODE CAPRA-MAURY, IL S'AGIT DE L'ORIENTATION DES FACETTES)
!     COND109    PRISE EN COMPTE OU PAS DE LA CLAUSE §109 POUR LE CALCUL
!                   DE RESISTANCE DANS LA METHODE SANDWICH
!                   1 = OUI, 2 = NON
!     TYPSTRU    TYPE DE STRUCTURE :
!                   0 = 2D, 1 = 1D
!     FERRSYME   FERRAILLAGE SYMETRIQUE?
!                   0 = NON, 1 = OUI
!     SLSYME     SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!     FERRCOMP   FERRAILLAGE DE COMPRESSION ADMIS?
!                   0 = NON, 1 = OUI
!     EPUCISA    IMPACT DE L'EFFORT TRANCHANT ET DE LA TORSION SUR LE
!                   FERRAILLAGE LONGITUDINAL?
!                   0 = NON, 1 = OUI
!     FERRMIN    PRISE EN COMPTE D'UN FERRAILLAGE MINIMUN
!                   0 = NON, 1 = OUI, 2 = CODE
!     FERRMINFISS   PRISE EN COMPTE DU Asmin Fiss ( 1 = EC2, 2 = RCC-CW)
!     EFFECHEL      PRISE EN COMPTE DE LEFFET DECHELLE ( 1 = OUI)
!     RHOLMIN    RATIO DE FERRAILLAGE LONGI MINI (A RENSEIGNER SI FERMIN='OUI')
!     RHOTMIN    RATIO DE FERRAILLAGE TRNSV MINI (A RENSEIGNER SI FERMIN='OUI')
!     COMPRESS   VALORISATION DE LA COMPRESSION POUR LES ACIERS TRANSVERSAUX
!                   0 = COMPRESSION NON PRISE EN COMPTE
!                   1 = COMPRESSION PRISE EN COMPTE
!     CEQUI      COEFFICIENT D'EQUIVALENCE ACIER/BETON
!     ENROBI     ENROBAGE DES ARMATURES INFERIEURES (2D)
!     ENROBS     ENROBAGE DES ARMATURES SUPERIEURES (2D)
!     ENROBYI    ENROBAGE DES ARMATURES INFERIEURES SUIVANT L'AXE Y (1D)
!     ENROBYS    ENROBAGE DES ARMATURES SUPERIEURES SUIVANT L'AXE Y (1D)
!     ENROBZI    ENROBAGE DES ARMATURES INFERIEURES SUIVANT L'AXE Z (1D)
!     ENROBZS    ENROBAGE DES ARMATURES SUPERIEURES SUIVANT L'AXE Z (1D)
!     SIGS       CONTRAINTE ULTIME DES ACIERS À L'ELS
!     SGICI      CONTRAINTE ULTIME DU BÉTON COMPRIME EN FIBRE INFERIEURE À L'ELS (2D)
!     SGICS      CONTRAINTE ULTIME DU BÉTON COMPRIME EN FIBRE SUPERIEURE À L'ELS (2D)
!     SGICYI     CONTRAINTE ULTIME DU BÉTON COMPRIME
!                   EN FIBRE INFERIEURE SUIVANT L'AXE Y À L'ELS (1D)
!     SGICYS     CONTRAINTE ULTIME DU BÉTON COMPRIME
!                   EN FIBRE SUPERIEURE SUIVANT L'AXE Y À L'ELS (1D)
!     SGICZI     CONTRAINTE ULTIME DU BÉTON COMPRIME
!                   EN FIBRE INFERIEURE SUIVANT L'AXE Z À L'ELS (1D)
!     SGICZS     CONTRAINTE ULTIME DU BÉTON COMPRIME
!                   EN FIBRE SUPERIEURE SUIVANT L'AXE Z À L'ELS (1D)
!     ALPHACC    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DU BETON EN COMPRESSION
!     GAMMAS     COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DES ACIERS
!     GAMMAC     COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DU BETON
!     FACIER     LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!     EYS        MODULE D'YOUNG DE L'ACIER
!     TYPDIAG    TYPE DE DIAGRAMME UTILISÉ POUR L'ACIER
!                   TYPDIAG = 1 ("B1" ==> PALIER INCLINÉ)
!                   TYPDIAG = 2 ("B2" ==> PALIER HORIZONTAL)
!     FBETON     RESISTANCE EN COMPRESSION DU BETON (CONTRAINTE)
!     CLACIER    CLASSE DE DUCTILITE DES ACIERS (POUR L'EC2) :
!                   CLACIER = 0 ACIER PEU DUCTILE (CLASSE A)
!                   CLACIER = 1 ACIER MOYENNEMENT DUCTILE (CLASSE B)
!                   CLACIER = 2 ACIER FORTEMENT DUCTILE (CLASSE C)
!     UC         UNITE DES CONTRAINTES :
!                   0 = CONTRAINTES EN Pa
!                   1 = CONTRAINTES EN MPa
!     UM         UNITE DES DIMENSIONS :
!                   0 = DIMENSIONS EN m
!                   1 = DIMENSIONS EN mm
!     RHOACIER   MASSE VOLUMIQUE DES ACIERS
!     AREINF     COEFF DE PONDER DU RATIO DE DENSITÉ D'ACIER PAR MÈTRE CUBE DE BÉTON
!     ASHEAR     COEFF DE PONDER DU RATIO DE DENSITÉ D'ACIER D'EFFORT TRANCHANT
!     ASTIRR     COEFF DE PONDER DU RATIO DE LONGUEUR DES ÉPINGLES D'ACIER EFF TRANC
!     RHOCRIT    DENSITÉ VOLUMIQUE D'ARMATURE CRITIQUE
!     DATCRIT    FERRAILLAGE D'EFFORT TRANCHANT CRITIQUE
!     LCRIT      LONGUEUR CRITIQUE DES EPINGLE D'ACIERS D'EFFORT TRANCHANT
!     WMAXI      OUVERTURE MAXIMALE DES FISSURES EN FACE INFÉRIEURE (2D)
!     WMAXS      OUVERTURE MAXIMALE DES FISSURES EN FACE SUPÉRIEURE (2D)
!     WMAXYI     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE INFÉRIEURE SUIVANT L'AXE Y (1D)
!     WMAXYS     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE SUPÉRIEURE SUIVANT L'AXE Y (1D)
!     WMAXZI     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE INFÉRIEURE SUIVANT L'AXE Z (1D)
!     WMAXZS     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE SUPÉRIEURE SUIVANT L'AXE Z (1D)
!     SIGELSQP   CONTRAINTE ADMISSIBLE DANS LE BETON À L'ELS QP
!     KT         COEFFICIENT DE DURÉE DE CHARGEMENT
!     PHIXI      DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT X
!     PHIXS      DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT X
!     PHIYI      DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT Y
!     PHIYS      DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT Y
!     PHIZI      DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT Z
!     PHIZS      DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT Z

!     HT         EPAISSEUR DE LA COQUE
!     EFFRTS     TORSEUR DES EFFORTS ET DES MOMENTS (DIM 8)
!
!   PARAMETRES DE SORTIE (RENVOYES A CODE_ASTER)
!
!     DNSITS     DENSITES DE FERRAILLAGE (DIM 5) :
!                   1 A 4 : ACIER LONGITUDINAL (EN M2/M)
!                   5 A 6 : ACIERS TRANSVERSAUX (EN M2/M2)
!     DNSVOL     DENSITE VOLUMIQUE D'ARMATURE (Kg/M3)
!     CONSTRUC   INDICATEUR DE COMPLEXITE DE CONSTRUCTIBILITE (-)
!     IERRL      CODE RETOUR LONGI (0 = OK)
!     IERRT      CODE RETOUR TRNSV (0 = OK)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: cequi, sigs, sigci, sigcs, sigcyi, sigcys, sigczi, sigczs
    real(kind=8) :: alphacc, effrts(8), dnsits(6), dnsmin(4)
    real(kind=8) :: ht, enrobi, enrobs, enrobyi, enrobys, enrobzi, enrobzs
    real(kind=8) :: gammac, gammas, rholmin, rhotmin, slsyme
    real(kind=8) :: facier, fbeton, eys, rhoacier, dnsvol, areinf, ashear
    real(kind=8) :: astirr, rhocrit, datcrit, lcrit, construc
    real(kind=8) :: wmaxi, wmaxs, wmaxyi, wmaxys, wmaxzi, wmaxzs, sigelsqp, kt
    real(kind=8) :: phixi, phixs, phiyi, phiys, phizi, phizs
    real(kind=8) :: reinf, shear, stirrups, thiter, epiter, aphiter
    integer(kind=8) :: ierrl, ierrt, jefge, jfer1, jfer2, itab(7)
    integer(kind=8) :: typcmb, typco, ferrmin, typdiag, ferrsyme, epucisa
    integer(kind=8) ::  clacier, uc, um, ferrminfiss, effechel
    integer(kind=8) :: ino, icmp, iret, comp, meth2D, cond109, precs, nno
    integer(kind=8) :: compress, ferrcomp, typstru, nb
    real(kind=8) :: k, hs, N, M, Asmininf, Asminsup, Act, kc, k1, sigma_c, m1
    real(kind=8) ::  chi0, fctm, unite_m, unite_pa, b
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PFERRA1', 'L', jfer1)
    call jevech('PFERRA2', 'E', jfer2)

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    ht = plateCara%thick
!
    call jevech('PEFFORR', 'L', jefge)
    call tecach('OOO', 'PEFFORR', 'L', iret, nval=7, itab=itab)
    ASSERT(iret .eq. 0)
    nno = itab(3)
    ASSERT(nno .gt. 0 .and. nno .le. 9)
    ASSERT(itab(2) .eq. 8*nno)

!
!       -- CALCUL DE LA CONTRAINTE MOYENNE :
!       ------------------------------------
!
    do icmp = 1, 8
        effrts(icmp) = 0.d0
        do ino = 1, nno
            effrts(icmp) = effrts(icmp)+zr(jefge-1+(ino-1)*8+icmp)/nno
        end do
    end do
!
!       -- RECUPERATION DES DONNEES DE L'UTILISATEUR :
!       ----------------------------------------------
!     FER1_R = 'TYPCOMB','CODIF','METH2D','THITER','EPITER','APHITER',
!                  1        2        3        4        5        6
!              'COND109','TYPSTRU','FERRSYME','SLSYME',
!                  7         8          9        10
!              'FERRCOMP','EPUCISA','FERRMIN','RHOLMIN','RHOTMIN',
!                  11        12         13        14       15
!              'COMPRESS','CEQUI','ENROBI','ENROBS','ENROBYI','ENROBYS',
!                  16        17      18       19       20        21
!              'ENROBZI','ENROBZS','SIGS','SIGCI','SIGCS','SIGCYI','SIGCYS',
!                  22        23      24      25      26      27       28
!              'SIGCZI','SIGCZS','ALPHACC','GAMMAS','GAMMAC','FACIER','EYS',
!                 29       30        31       32       33       34     35
!              'TYPDIAG','FBETON','CLACIER','UC','UM','RHOACIER','AREINF',
!                 36        37       38      39   40      41        42
!              'ASHEAR','ASTIRR','RHOCRIT','DATCRIT','LCRIT','WMAXI','WMAXS',
!                 43       44        45       46       47      48      49
!              'WMAXYI','WMAXYS','WMAXZI','WMAXZS','SIGELSQP','KT',
!                 50       51       52       53        54      55
!              'PHIXI','PHIXS','PHIYI','PHIYS','PHIZI','PHIZS', 'ASMINFIS',
!                 56      57      58      59      60      61       62
!              'ECHELLE'
!                 63
    typcmb = nint(zr(jfer1-1+1))
    typco = nint(zr(jfer1-1+2))
    meth2D = nint(zr(jfer1-1+3))
    thiter = zr(jfer1-1+4)
    epiter = zr(jfer1-1+5)
    aphiter = zr(jfer1-1+6)
    cond109 = nint(zr(jfer1-1+7))
    typstru = nint(zr(jfer1-1+8))
    ferrsyme = nint(zr(jfer1-1+9))
    slsyme = zr(jfer1-1+10)
    ferrcomp = nint(zr(jfer1-1+11))
    epucisa = nint(zr(jfer1-1+12))
    ferrmin = nint(zr(jfer1-1+13))
    rholmin = zr(jfer1-1+14)
    rhotmin = zr(jfer1-1+15)
    compress = int(zr(jfer1-1+16))
    cequi = zr(jfer1-1+17)
    enrobi = zr(jfer1-1+18)
    enrobs = zr(jfer1-1+19)
    enrobyi = zr(jfer1-1+20)
    enrobys = zr(jfer1-1+21)
    enrobzi = zr(jfer1-1+22)
    enrobzs = zr(jfer1-1+23)
    sigs = zr(jfer1-1+24)
    sigci = zr(jfer1-1+25)
    sigcs = zr(jfer1-1+26)
    sigcyi = zr(jfer1-1+27)
    sigcys = zr(jfer1-1+28)
    sigczi = zr(jfer1-1+29)
    sigczs = zr(jfer1-1+30)
    alphacc = zr(jfer1-1+31)
    gammas = zr(jfer1-1+32)
    gammac = zr(jfer1-1+33)
    facier = zr(jfer1-1+34)
    eys = zr(jfer1-1+35)
    typdiag = int(zr(jfer1-1+36))
    fbeton = zr(jfer1-1+37)
    clacier = int(zr(jfer1-1+38))
    uc = int(zr(jfer1-1+39))
    um = int(zr(jfer1-1+40))
    rhoacier = zr(jfer1-1+41)
    areinf = zr(jfer1-1+42)
    ashear = zr(jfer1-1+43)
    astirr = zr(jfer1-1+44)
    rhocrit = zr(jfer1-1+45)
    datcrit = zr(jfer1-1+46)
    lcrit = zr(jfer1-1+47)
    wmaxi = zr(jfer1-1+48)
    wmaxs = zr(jfer1-1+49)
    wmaxyi = zr(jfer1-1+50)
    wmaxys = zr(jfer1-1+51)
    wmaxzi = zr(jfer1-1+52)
    wmaxzs = zr(jfer1-1+53)
    sigelsqp = zr(jfer1-1+54)
    kt = zr(jfer1-1+55)
    phixi = zr(jfer1-1+56)
    phixs = zr(jfer1-1+57)
    phiyi = zr(jfer1-1+58)
    phiys = zr(jfer1-1+59)
    phizi = zr(jfer1-1+60)
    phizs = zr(jfer1-1+61)
    ferrminfiss = nint(zr(jfer1-1+62))
    effechel = nint(zr(jfer1-1+63))

    !Only option '2D'
    if (typstru .eq. 1.d0) then
        call utmess('A', 'CALCULEL_78')
        goto 998
    end if

!
!   -- CALCUL PROPREMENT DIT :
!   --------------------------
!
!   VERIFICATION DE LA COHERENCE DES PARAMETRES
    if (ht .le. enrobi .or. ht .le. enrobs) then
        call utmess('F', 'CALCULEL_72')
    end if
    if (meth2D .eq. 2.d0) then
    if ((epiter .le. 0.d0) .or. (epiter .gt. 0.2)) then
        call utmess('F', 'CALCULEL7_32')
    end if
    end if

!
!       -- INITIALISATION DES VARIABLES DE SORTIES :
!       --------------------------------------------
!
    dnsvol = 0.d0
    construc = 0.d0
    do comp = 1, 6
        dnsits(comp) = 0.d0
    end do
    do comp = 1, 4
        dnsmin(comp) = 0.d0
    end do

    if (meth2D .eq. 1.d0) then
        nb = ceiling(180/thiter)
        precs = ceiling(1/epiter)
        call clcplq(typcmb, typco, nb, precs, &
                    ferrsyme, slsyme, ferrcomp, epucisa, &
                    ferrmin, rholmin, rhotmin, compress, &
                    cequi, enrobi, enrobs, sigs, sigci, sigcs, &
                    alphacc, gammas, gammac, facier, eys, typdiag, &
                    fbeton, clacier, uc, um, &
                    wmaxi, wmaxs, sigelsqp, kt, phixi, phixs, phiyi, phiys, &
                    ht, effrts, dnsits, ierrl, ierrt)
    elseif (meth2D .eq. 2.d0) then
        if (typcmb .eq. 0) then
            if (typco .eq. 1.d0) then
                compress = 0
            end if
            call sandwich(enrobi, enrobs, facier, fbeton, gammas, gammac, &
                          thiter, epiter, aphiter, cond109, &
                          ferrcomp, ferrsyme, slsyme, &
                          epucisa, ferrmin, rholmin, rhotmin, compress, &
                          alphacc, eys, typdiag, clacier, &
                          uc, um, &
                          ht, effrts, dnsits, ierrl, ierrt)
        end if
    end if
!
!       -- CALCUL DE LA DENSITE VOLUMIQUE D'ARMATURE :
!       ----------------------------------------------
!
    if (rhoacier .gt. 0) then
        dnsvol = rhoacier*(dnsits(1)+dnsits(2)+dnsits(3)+dnsits(4) &
                           +dnsits(5)*ht+dnsits(6)*ht)/ht
        if (dnsits(5) .eq. -1.d0 .or. dnsits(6) .eq. -1.d0) then
!           Vrai uniquement pour le calcul du ferraillage transversal au BAEL
!           (pour lequel les aciers d'effort tranchant ne sont pas calculés)
            dnsvol = rhoacier*((dnsits(1)+dnsits(2)+dnsits(3)+dnsits(4)))/ht
        end if
    else
        dnsvol = -1.d0
    end if
!
!       -- CALCUL DE L'INDICATEUR DE CONSTRUCTIBILITE :
!       -----------------------------------------------
!
    if (rhoacier .gt. 0) then
        reinf = areinf*dnsvol/rhocrit
        shear = ashear*dnsits(5)/datcrit
        if (shear .lt. 0.d0) shear = 0.d0
        stirrups = astirr*((dnsits(5)*dnsits(5)+dnsits(6)*dnsits(6))**0.5)
        stirrups = stirrups*(ht-enrobs-enrobi)/(datcrit*lcrit)
        if (stirrups .lt. 0.d0) stirrups = 0.d0
        construc = (reinf+shear+stirrups)/(areinf+ashear+astirr)
    else
        construc = -1.d0
    end if

!       -- GESTION DES ALARMES EMISES :
!       -------------------------------

    if (ierrl .eq. 1001) then
!       ELU : section trop comprimée
        call utmess('A', 'CALCULEL_83')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 10011) then
!       ELU : ferraillage symétrique non possible!
        call utmess('A', 'CALCULEL7_28')
        dnsvol = -1.d0
        construc = -1.d0
    end if
!
    if (ierrl .eq. 1003) then
!       ELS : section trop comprimée
        call utmess('A', 'CALCULEL_84')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 1005) then
!       ELS_QP : section trop comprimée
        call utmess('A', 'CALCULEL_85')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrt .eq. 1002) then
!       ELU BETON TROP CISAILLE : densité transversale fixée à -1 pour l'élément
        call utmess('A', 'CALCULEL_81')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrt .eq. 1004) then
!       ELS BETON TROP CISAILLE : densité transversale fixée à -1 pour l'élément
        call utmess('A', 'CALCULEL_86')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrt .eq. 1007) then
!       ELS_QP BETON TROP CISAILLE : densité transversale fixée à -1 pour l'élément
        call utmess('A', 'CALCULEL_87')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 1006) then
!       ELS QP SOLLICITATION TROP IMPORTANTE : Résolution itérative impossible à l'els qp !
        call utmess('A', 'CALCULEL_77')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 2001) then
!       ELU-Sandwich : équilibre non possible
!       on passe par Capra-Maury...
        call utmess('A', 'CALCULEL7_34')
        !dnsvol = -1.d0
        !construc = -1.d0
    end if

    if (ierrl .eq. 2002) then
!       ELU-Sandwich : section trop comprimée
!       on fixe toutes les densités de ferraillage de l'élément à -1
        call utmess('A', 'CALCULEL7_35')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 2003) then
!       ELU-Sandwich : ferraillage symétrique non possible!
!       on fixe toutes les densités de ferraillage de l'élément à -1
        call utmess('A', 'CALCULEL7_28')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrt .eq. 2004) then
!       ELU-Sandwich BETON TROP CISAILLE : densité transversale fixée à -1 pour l'élément
        call utmess('A', 'CALCULEL7_36')
        dnsvol = -1.d0
        construc = -1.d0
    end if
!

    ! Minimum steel for crack control
    if ((ferrminfiss .eq. 1) .or. ferrminfiss .eq. 2) then
        if (typcmb .eq. 2) then

            if (uc .eq. 0) then
                unite_pa = 1.e6
                unite_m = 1.e-3
            elseif (uc .eq. 1) then
                unite_pa = 1.
                unite_m = 1.
            end if
            ! Convert length to mm and stress to MPa
            facier = facier/unite_pa
            fbeton = fbeton/unite_pa
            ht = ht/unite_m
            b = 1000.d0
            if (fbeton .le. 50.d0) then
                fctm = 0.30*((fbeton)**(2.0/3.0))
            else
                fctm = 2.12*LOG(1.0+(fbeton+8.0)/10.0)
            end if

            ! computation of k
            if (ht .le. 300.d0) then
                k = 1.d0
            else if (ht .ge. 800.d0) then
                k = 0.65d0
            else
                k = (0.65d0-1.d0)/(800.d0-300.d0)*(ht-300.d0)+1.d0
            end if

            ! computation of hs
            hs = min(ht, 1000.d0)

            do comp = 1, 2
                Asmininf = 0.d0
                Asminsup = 0.d0
                N = -effrts(comp)
                M = -effrts(comp+3)/unite_m

                ! Taking into account scale effect
                if (effechel .eq. 1) then

                    if ((ht .le. 6250.d0) .and. (ht .ge. 400.d0) .and. &
                        (fbeton .le. 55.d0) .and. (fbeton .ge. 30.d0)) then
                        !Scale effect is applicable according to DN 2100 RCC-CW

                        kc = 1.d0
                        Act = b*ht/2.d0
                        m1 = 2.5d0*0.1d0-3.6d0*0.001d0*(fbeton+8) &
                             +1.3d0*0.00001d0*(fbeton+8)**2d0

                        chi0 = (6.d0*0.001d0/(ht/1000.d0)**2.d0)**m1
                        Asmininf = kc*k*chi0*fctm*Act/facier
                        Asminsup = Asmininf

                    else
                        !Scale effect is not applicable according to DN 2100 RCC-CW
                        call utmess('A', 'CALCULEL7_38')
                    end if

                elseif ((N .gt. 0.d0) .and. (abs(M/N) .lt. ht/6.d0)) then
                    ! Full compression
                    Asmininf = 0.d0
                    Asminsup = 0.d0
                elseif ((N .lt. 0.d0) .and. (abs(-M/N) .lt. ht/6.d0)) then
                    ! Full tension

                    Act = b*ht
                    kc = 1.d0
                    Asmininf = kc*k*fctm*Act/facier/2.d0
                    Asminsup = Asmininf

                else

                    ! partially compressed
                    !Area of concrete in tension
                    ! to make it simple we consider hcteff = h/2
                    Act = b*ht/2.d0
                    if (N .gt. 0.d0) then
                        k1 = 1.5d0
                    else
                        k1 = 2.d0*hs/(3.d0*ht)
                    end if

                    sigma_c = N/(b*ht)

                    kc = max(0.d0, min(1.d0, 0.4d0*(1-(sigma_c/(k1*(ht/hs)*fctm)))))

                    if (M .gt. 0.d0) then
                        ! lower fiber in tension
                        Asmininf = kc*k*fctm*Act/facier
                        Asminsup = 0.d0
                    elseif (M .lt. 0.d0) then
                        ! upper fiber in tension
                        Asminsup = kc*k*fctm*Act/facier
                        Asmininf = 0.d0
                    else
                        !No miminium steel
                        Asmininf = 0.d0
                        Asminsup = 0.d0
                    end if

                end if

                dnsmin(2*comp-1) = Asmininf*(unite_m)**2.d0
                dnsmin(2*comp) = Asminsup*(unite_m)**2.d0

            end do

            !Taking into account minimum steel
            dnsits(1) = max(dnsits(1), dnsmin(1))
            dnsits(2) = max(dnsits(2), dnsmin(2))
            dnsits(3) = max(dnsits(3), dnsmin(3))
            dnsits(4) = max(dnsits(4), dnsmin(4))

        else
            call utmess('A', 'CALCULEL7_39')

        end if

    end if

!
!       -- STOCKAGE DES RESULTATS DANS FER2 :
!       -------------------------------------
!   FER2_R =  DNSXI DNSXS DNSYI DNSYS DNSXT DNSYT DNSVOL CONSTRUC ASMINXI ASMINXS ASMINYI ASMINYS
!               1     2     3     4     5     6      7       8      9      10     11     12
!
    zr(jfer2-1+1) = dnsits(1)
    zr(jfer2-1+2) = dnsits(3)
    zr(jfer2-1+3) = dnsits(2)
    zr(jfer2-1+4) = dnsits(4)
    zr(jfer2-1+5) = dnsits(5)
    zr(jfer2-1+6) = dnsits(6)
    zr(jfer2-1+7) = dnsvol
    zr(jfer2-1+8) = construc
    zr(jfer2-1+9) = dnsmin(1)
    zr(jfer2-1+10) = dnsmin(2)
    zr(jfer2-1+11) = dnsmin(3)
    zr(jfer2-1+12) = dnsmin(4)
998 continue

end subroutine
