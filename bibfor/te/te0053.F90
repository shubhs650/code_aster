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
subroutine te0053(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vrfplq.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!  BUT: VERIFICATION DE FERRAILLAGE POUR LES ELEMENTS DE PLAQUE
!
! DIAGRAMME D INTERACTION ET CALCUL DE LA MARGE
!              (METHODE DE CAPRA ET MAURY)
!
! --------------------------------------------------------------------------------------------------
!
! PARAMETRES D'ECHANGE ENTRE CODE_ASTER ET VRFPLQ
! (POINT D'ENTREE DU CALCUL DE FERRAILLAGE PAR CAPRA ET MAURY)
!
!   PARAMETRES D'ENTREE (FOURNIS PAR CODE_ASTER)
!
!     TYPCMB     TYPE DE COMBINAISON :
!                   0 = ELU, 1 = ELS, 2 = ELS QP
!     TYPCO      TYPE DE CODIFICATION :
!                   1 = BAEL91, 2 = EUROCODE 2
!     THITER     ANGLE D'ITERATION POUR LA METHODE CAPRA-MAURY,
!     IL S'AGIT DE L'ORIENTATION DES FACETTES)
!     TYPSTRU    TYPE DE STRUCTURE :
!                   0 = 2D, 1 = 1D
!                   0 = NON, 1 = OUI
!     ENROBI     ENROBAGE DES ARMATURES INFERIEURES (2D)
!     ENROBS     ENROBAGE DES ARMATURES SUPERIEURES (2D)
!     SIGS       CONTRAINTE ULTIME DES ACIERS À L'ELS
!     SGICI      CONTRAINTE ULTIME DU BÉTON COMPRIME EN FIBRE INFERIEURE À L'ELS (2D)
!     SGICS      CONTRAINTE ULTIME DU BÉTON COMPRIME EN FIBRE SUPERIEURE À L'ELS (2D)
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
!     HT         EPAISSEUR DE LA COQUE
!     EFFRTS     TORSEUR DES EFFORTS DE CALCUL ET DES MOMENTS (DIM 8)
!     EFFREF     TORSEUR DES EFFORTS DE REFERENCE ET DES MOMENTS (DIM 8)
!
!   PARAMETRES DE SORTIE (RENVOYES A CODE_ASTER)
!
!   Resultats sur la facette critique (correspond a la marge minimale)
!     marge     la marge mecanique de la facette critique
!     Nrd,Mrd   Efforts N,M : projection de la droite des efforts sur le diagramme NM
!     tau_crit  Taux d utilisation sur la facette critique : 1-marge
!     c0_c      Distance C0C entre le point de ref. et le point d etude
!     c0_crd    Distance C0CRD entre le point de ref. et le point du diagramme
!     + tous les parametres necessaires a la production du diagramme d interaction
!     sur la facette critique
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: sigs, sigci, sigcs, marge
    real(kind=8) :: alphacc, effrts(8), ht, enrobi, enrobs, effref(8)
    real(kind=8) :: gammac, gammas, thiter
    real(kind=8) :: facier, fbeton, eys
    real(kind=8) :: dnsxi, dnsyi, dnsxs, dnsys
    real(kind=8) :: cequi
    integer(kind=8) :: jefge, jefge0, jfer0, jfer1, jfer2
    integer(kind=8) :: itab(7), nno, typcmb, typco, typdiag, clacier, uc, um
    integer(kind=8) :: ino, icmp, iret
    integer(kind=8) :: typstru, nb
    real(kind=8) :: dnsinf_crit, dnssup_crit, myNrd_crit, myMrd_crit
    real(kind=8) :: tau_crit, c0c_crit, c0crd_crit
    real(kind=8) :: theta_crit, effn_crit, effm_crit, effn0_crit, effm0_crit, bw
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    ht = plateCara%thick

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

    call jevech('PVFER0', 'L', jfer0)
    call jevech('PVFER1', 'L', jfer1)
    call jevech('PVFER2', 'E', jfer2)

!
    call jevech('PEFFORR', 'L', jefge)
    call tecach('OOO', 'PEFFORR', 'L', iret, nval=7, itab=itab)
    ASSERT(iret .eq. 0)
    nno = itab(3)
    ASSERT(nno .gt. 0 .and. nno .le. 9)
    ASSERT(itab(2) .eq. 8*nno)

    call jevech('PEFFOR0', 'L', jefge0)
    call tecach('OOO', 'PEFFOR0', 'L', iret, nval=7, itab=itab)
    ASSERT(iret .eq. 0)
    ASSERT(nno .eq. itab(3))
!
!       -- CALCUL DE LA CONTRAINTE MOYENNE ET DE LA CONTRAINTE MOYENNE DE REFERENCE :
!       ------------------------------------
!
    do icmp = 1, 8
        effrts(icmp) = 0.d0
        effref(icmp) = 0.d0
        do ino = 1, nno
            effrts(icmp) = effrts(icmp)+zr(jefge-1+(ino-1)*8+icmp)/nno
            effref(icmp) = effref(icmp)+zr(jefge0-1+(ino-1)*8+icmp)/nno
        end do
    end do
!
!       -- RECUPERATION DES DONNEES DE L'UTILISATEUR :
!       ----------------------------------------------
!     VFER1_R = 'TYPCOMB','CODIF', 'THITER', 'TYPSTRU','ENROBI','ENROBS','SIGS'
!                  1        2        3        4          5        6        7
!               'SIGCI', 'SIGCS', 'ALPHACC','GAMMAS','GAMMAC','FACIER'
!                  8        9         10       11       12        13
!                 'EYS', 'TYPDIAG','FBETON','CLACIER','UC',     'UM',
!                  14        15      16        17      18        19
!               'CEQUI'
!                  20
!     VFER2_R = 'DNSXI','DNSXS', 'DNSYI', 'DNSYS'
!                  1       2        3        4
!
    typcmb = nint(zr(jfer1-1+1))
    typco = nint(zr(jfer1-1+2))
    thiter = zr(jfer1-1+3)
    typstru = nint(zr(jfer1-1+4))
    enrobi = zr(jfer1-1+5)
    enrobs = zr(jfer1-1+6)
    sigs = zr(jfer1-1+7)
    sigci = zr(jfer1-1+8)
    sigcs = zr(jfer1-1+9)
    alphacc = zr(jfer1-1+10)
    gammas = zr(jfer1-1+11)
    gammac = zr(jfer1-1+12)
    facier = zr(jfer1-1+13)
    eys = zr(jfer1-1+14)
    typdiag = int(zr(jfer1-1+15))
    fbeton = zr(jfer1-1+16)
    clacier = int(zr(jfer1-1+17))
    uc = int(zr(jfer1-1+18))
    um = int(zr(jfer1-1+19))
    cequi = zr(jfer1-1+20)
    dnsxi = zr(jfer0-1+1)
    dnsxs = zr(jfer0-1+2)
    dnsyi = zr(jfer0-1+3)
    dnsys = zr(jfer0-1+4)

!   --------------------------
!
!   VERIFICATION DE LA COHERENCE DES PARAMETRES
    if (ht/2.d0 .le. enrobi .or. ht/2.d0 .le. enrobs) then
        call utmess('F', 'VERIFERRAILLAGE_10')
    end if

! !
! !       -- INITIALISATION DES VARIABLES DE SORTIES :
! !       --------------------------------------------
! !
    nb = ceiling(180/thiter)
    call vrfplq(typcmb, typco, nb, cequi, enrobi, enrobs, sigs, sigci, &
                sigcs, alphacc, gammas, gammac, facier, eys, &
                typdiag, fbeton, clacier, uc, ht, effrts, &
                effref, dnsxi, dnsxs, dnsyi, dnsys, marge, theta_crit, &
                dnsinf_crit, dnssup_crit, effn0_crit, effm0_crit, effn_crit, &
                effm_crit, myNrd_crit, myMrd_crit, bw, tau_crit, c0c_crit, c0crd_crit)

!
!       -- STOCKAGE DES RESULTATS DANS VFER2 :
!       -------------------------------------
    zr(jfer2-1+1) = marge
    zr(jfer2-1+2) = theta_crit
    zr(jfer2-1+3) = tau_crit
    zr(jfer2-1+4) = c0c_crit
    zr(jfer2-1+5) = c0crd_crit
    zr(jfer2-1+6) = dnsinf_crit
    zr(jfer2-1+7) = dnssup_crit
    zr(jfer2-1+8) = effn0_crit
    zr(jfer2-1+9) = effm0_crit
    zr(jfer2-1+10) = effn_crit
    zr(jfer2-1+11) = effm_crit
    zr(jfer2-1+12) = myNrd_crit
    zr(jfer2-1+13) = myMrd_crit
    zr(jfer2-1+14) = typcmb
    zr(jfer2-1+15) = typco
    zr(jfer2-1+16) = cequi
    zr(jfer2-1+17) = enrobi
    zr(jfer2-1+18) = enrobs
    zr(jfer2-1+19) = sigs
    zr(jfer2-1+20) = sigci
    zr(jfer2-1+21) = sigcs
    zr(jfer2-1+22) = alphacc
    zr(jfer2-1+23) = gammas
    zr(jfer2-1+24) = gammac
    zr(jfer2-1+25) = facier
    zr(jfer2-1+26) = eys
    zr(jfer2-1+27) = typdiag
    zr(jfer2-1+28) = fbeton
    zr(jfer2-1+29) = clacier
    zr(jfer2-1+30) = uc
    zr(jfer2-1+31) = um
    zr(jfer2-1+32) = ht
    zr(jfer2-1+33) = bw

end subroutine
