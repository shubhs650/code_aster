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
subroutine bsthco(plateCara, plateOrie, &
                  nomte, bsigth)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/btsig.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/matrc2.h"
#include "asterfort/moytem.h"
#include "asterfort/promat.h"
#include "asterfort/utmess.h"
#include "asterfort/verifm.h"
#include "asterfort/vexpan.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(out) :: bsigth(51)
!
! --------------------------------------------------------------------------------------------------
!
!      CALCUL DU BSIGMA POUR LES CONTRAINTES THERMIQUES
!      (I.E. BT*D*ALPHA(T-TREF)) POUR LES ELEMENTS DE COQUE (COQUE_3D)
!
! --------------------------------------------------------------------------------------------------
!
!     IN  NOMTE  : NOM DU TYPE D'ELEMENT
!     OUT BSIGTH : BT*SIGMA POUR LES CONTRAINTES THERMIQUES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npge = 3
    real(kind=8), parameter :: zero = 0.0d0, un = 1.0d0, deux = 2.0d0, trois = 3.0d0, quatre = 4.0d0
    integer(kind=8) :: iLayer, jvMaterc, inte, intsn, intsr, jvGeom, lzi
    integer(kind=8) :: lzr, nb1, nb2, nbLayer, npgsn, npgsr, kwgt, iret
    real(kind=8) :: vectTangKpg(2, 3), vectBaseKpg(3, 3)
    real(kind=8) :: hsfm(3, MT_NNOMAX2D), hss(2, 9), hsj1m(3, MT_NNOMAX2D), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, MT_NNOMAX2D), hsj1fx(3, MT_NNOMAX2D), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: epsth(5), sigmth(5), bsigt1(42)
    real(kind=8) :: ksi3s2, kappa, matrElas(5, 5)
    real(kind=8) :: coef, hLayer, eptot, tempMoye
    real(kind=8) :: zic, zmin, epsthe
!
! --------------------------------------------------------------------------------------------------
!
    epsth = zero
    sigmth = zero
    bsigt1 = zero
    bsigth = zero

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get properties of shell
    nbLayer = plateCara%nbLayer

! - Very strange (see issue35983)
    nbLayer = 1

    eptot = plateCara%thick
    kappa = plateCara%shearCoef
    zmin = -eptot/deux
    hLayer = eptot/nbLayer

! - RECUPERATION DES COORDONNEES DES NOEUDS DANS LA GEOMETRIE INITIALE
    call jevech('PGEOMER', 'L', jvGeom)

! - RECUPERATION DU MATERIAU
    call jevech('PMATERC', 'L', jvMaterc)

    kwgt = 0
    do iLayer = 1, nbLayer
        do inte = 1, npge
! ---      POSITION DANS L'EPAISSEUR :
            if (inte .eq. 1) then
                zic = zmin+(iLayer-1)*hLayer
                coef = un/trois
            else if (inte .eq. 2) then
                zic = zmin+hLayer/deux+(iLayer-1)*hLayer
                coef = quatre/trois
            else if (inte .eq. 3) then
                zic = zmin+hLayer+(iLayer-1)*hLayer
                coef = un/trois
            end if

! ---      COORDONNEE ISOPARAMETRIQUE DANS L'EPAISSEUR DIVISEE PAR 2
            ksi3s2 = zic/hLayer
!
! ---      CALCUL POUR L'INTEGRATION REDUITE DES PARTIES MEMBRANE
! ---      BTDM ET CISAILLEMENT BTDS DE LA MATRICE B
            do intsr = 1, npgsr
!
! ---       .D'UNE PART :
! ---        DETERMINATION DES REPERES LOCAUX AUX POINTS D'INTEGRATION
! ---        DANS LA CONFIGURATION INITIALE
! ---        VECTG DESIGNE LES VECTEURS COVARIANTS DANS LE PLAN MOYEN
! ---              EN CHAQUE POINT D'INTEGRATION
! ---        VECTT DESIGNE LES REPERES LOCAUX ORTHORNORMES EN CHAQUE
! ---        POINT D'INTEGRATION DANS LA CONFIGURATION INITIALE
! ---       .D'AUTRE-PART :
! ---        SOIT H LA MATRICE DE PASSAGE DU TENSEUR DE GREEN-LAGRANGE
! ---        DU REPERE LOCAL AU REPERE GLOBAL
! ---        SOIT S LA MATRICE CONSTANTE TELLE QUE [S]*(DU/DX)
! ---        REPRESENTE LA PARTIE LINEAIRE DU TENSEUR DE GREEN-LAGRANGE
! ---        ON CALCULE LES PRODUITS [HSFM] = [H]*[S] POUR LA PARTIE
! ---                                MEMBRANE-FLEXION
! ---                                [HSS]  = [H] * [S] POUR LA PARTIE
! ---                                CISAILLEMENT
                call mahsms(plateOrie, &
                            0, nb1, &
                            zr(jvGeom), ksi3s2, intsr, &
                            zr(lzr), hLayer, &
                            vectBaseKpg, vectTangKpg, &
                            hsfm, hss)

! ---       MULTIPLICATION DES MATRICES [HSFM] ET [HSS] PAR L'INVERSE
! ---       DE LA MATRICE JACOBIENNE [JM1]:
! ---       [HSJ1M] = [HSFM]*[JM1] , [HSJ1S] = [HSS]*[JM1]
                call hsj1ms(hLayer, vectTangKpg, vectBaseKpg, hsfm, hss, &
                            hsj1m, hsj1s)
!
! ---       CALCUL POUR L'INTEGRATION REDUITE DES PARTIES MEMBRANE
! ---       BTDM ET CISAILLEMENT BTDS DE LA MATRICE B :
!           -----------------------------------------
                call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                            hLayer, plateOrie%vectTang, hsj1m, hsj1s, btdm, &
                            btds)
            end do
!
! ---      CALCUL POUR L'INTEGRATION NORMALE DE LA PARTIE FLEXION
! ---      BTDFN DE LA MATRICE B
            do intsn = 1, npgsn
!
! ---       .D'UNE PART :
! ---        DETERMINATION DES REPERES LOCAUX AUX POINTS D'INTEGRATION
! ---        DANS LA CONFIGURATION INITIALE
! ---        VECTG DESIGNE LES VECTEURS COVARIANTS DANS LE PLAN MOYEN
! ---              EN CHAQUE POINT D'INTEGRATION
! ---        VECTT DESIGNE LES REPERES LOCAUX ORTHORNORMES EN CHAQUE
! ---        POINT D'INTEGRATION DANS LA CONFIGURATION INITIALE
! ---       .D'AUTRE-PART :
! ---        SOIT H LA MATRICE DE PASSAGE DU TENSEUR DE GREEN-LAGRANGE
! ---        DU REPERE LOCAL AU REPERE GLOBAL
! ---        SOIT S LA MATRICE CONSTANTE TELLE QUE [S]*(DU/DX)
! ---        REPRESENTE LA PARTIE LINEAIRE DU TENSEUR DE GREEN-LAGRANGE
! ---        ON CALCULE LE PRODUIT [HSF] = [H]*[S] POUR LA PARTIE
! ---                              FLEXION
                call mahsf(plateOrie, &
                           1, nb1, &
                           zr(jvGeom), ksi3s2, intsn, &
                           zr(lzr), hLayer, &
                           vectBaseKpg, vectTangKpg, &
                           hsf)
!
! ---       MULTIPLICATION DE LA MATRICE [HSF] PAR L'INVERSE
! ---       DE LA MATRICE JACOBIENNE [JM1]:
! ---       [HSJ1FX] = [HSF]*[JM1]
                call hsj1f(intsn, zr(lzr), hLayer, vectTangKpg, vectBaseKpg, &
                           hsf, kwgt, hsj1fx, wgt)
!
! ---       PRODUIT DU POIDS DU POINT DE GAUSS DANS L'EPAISSEUR PAR WGT
! ---       QUI EST LE PRODUIT DU POIDS DU POINT DE GAUSS COURANT
! ---       DANS L'EPAISSEUR PAR LE JACOBIEN EN CE POINT :
!           --------------------------------------------
                wgt = coef*wgt
!
! ---       CALCUL POUR L'INTEGRATION NORMALE DE LA PARTIE FLEXION
! ---       BTDF DE LA MATRICE B :
!           --------------------
                call btdfn(1, nb1, nb2, ksi3s2, intsn, &
                           zr(lzr), hLayer, plateOrie%vectTang, hsj1fx, btdf)
!
! ---       CALCUL DE LA MATRICE B [BTILD] PAR INTEGRATION SELECTIVE
! ---       ET INSERTION DES PARTIES [BTDM] ET [BDTS] ET INSERTION
! ---       DE LA PARTIE [BTDF]  :
!           -------------------
                call btdmsn(1, nb1, intsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)
!
! ---       EVALUATION DES DEFORMATIONS THERMIQUES :
!           ======================================
                call verifm('RIGI', inte, 3, '+', zi(jvMaterc), epsthe)
                call moytem('RIGI', inte, 3, '+', tempMoye, iret)
                epsth(1:2) = epsthe

! ------------- Elastic matrix
                call matrc2(plateOrie, vectBaseKpg, tempMoye, kappa, matrElas)
!
! ---       CALCUL DES CONTRAINTES THERMIQUES SIGMTH(5) :
!           -------------------------------------------
                call promat(matrElas, 5, 5, 5, epsth, &
                            5, 5, 1, sigmth)
!
! ---       CALCUL DES FORCES INTERNES DUES AUX CONTRAINTES THERMIQUES :
!           ----------------------------------------------------------
                call btsig(5*nb1+2, 5, wgt, btild, sigmth, &
                           bsigt1)
            end do
        end do
    end do
!
! --- EXPANSION DE BSIGT1 DANS BSIGTH :
!     -------------------------------
    call vexpan(nb1, bsigt1, bsigth)
    bsigth(6*nb1+1) = bsigt1(5*nb1+1)
    bsigth(6*nb1+2) = bsigt1(5*nb1+2)
!
end subroutine
