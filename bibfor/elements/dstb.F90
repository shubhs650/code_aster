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
subroutine dstb(plateCara, plateOrie, &
                carat3, igau, jacgau, bmat)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/bcoqaf.h"
#include "asterfort/dstbfa.h"
#include "asterfort/dstbfb.h"
#include "asterfort/dstcis.h"
#include "asterfort/dsxhft.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxtbm.h"
#include "asterfort/elrefe_info.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8), intent(in) :: carat3(*)
    integer(kind=8), intent(in) :: igau
    real(kind=8), intent(out) :: bmat(8, 1), jacgau
! --- CALCUL DE LA MATRICE (B) RELIANT LES DEFORMATIONS DU PREMIER
! --- ORDRE AUX DEPLACEMENTS AU POINT D'INTEGRATION D'INDICE IGAU
! --- POUR UN ELEMENT DE TYPE DST
! --- (I.E. (EPS_1) = (B)*(UN))
! --- D'AUTRE_PART, ON CALCULE LE PRODUIT NOTE JACGAU = JACOBIEN*POIDS
!     ------------------------------------------------------------------
!     IN  PGL(3,3)      : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                         LOCAL
!     IN  IGAU          : INDICE DU POINT D'INTEGRATION
!     OUT JACGAU        : PRODUIT JACOBIEN*POIDS AU POINT D'INTEGRATION
!                         COURANT
!     OUT BMAT(8,1)     : MATRICE (B) AU POINT D'INTEGRATION COURANT
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: multic, i, j, k
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bfb(3, 9), bfa(3, 3), bfn(3, 9), bm(3, 6), bf(3, 9)
    real(kind=8) :: bca(2, 3), bcn(2, 9), bc(2, 9)
    real(kind=8) :: hft2(2, 6), an(3, 9)
    real(kind=8) :: qsi, eta
    aster_logical :: coupmf
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
! --- COORDONNEES DU POINT D'INTEGRATION COURANT :
!     ------------------------------------------
    qsi = zr(icoopg-1+ndim*(igau-1)+1)
    eta = zr(icoopg-1+ndim*(igau-1)+2)
!
! --- PRODUIT JACOBIEN*POIDS
!     ----------------------
    jacgau = zr(ipoids+igau-1)*carat3(7)
!
! --- CALCUL DES MATRICES DE HOOKE DE FLEXION, MEMBRANE,
! --- MEMBRANE-FLEXION, CISAILLEMENT, CISAILLEMENT INVERSE
!     ----------------------------------------------------
    call dxmate(plateCara, plateOrie, &
                'RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)
!
! --- CALCUL DE LA MATRICE NOTEE (HF.T2) PAR BATOZ RELIANT LES
! --- EFFORTS TRANCHANTS (T) AU VECTEUR DES DERIVEES DES COURBURES
! --- PAR RAPPORT AUX COORDONNEES PARAMETRIQUES ,
! --- SOIT (T) = (HF.T2)*(D (D_BETA/D_QSI)/D_QSI)
!     -----------------------------------------------------------
    call dsxhft(df, carat3(9), hft2)
!
! --- CALCUL DE LA MATRICE (AN) RELIANT LES INCONNUES NOTEES (ALFA)
! --- PAR BATOZ AUX INCONNUES (UN) = (...,W_I,BETAX_I,BETAY_I,...)
! --- SOIT (ALFA) = (AN)*(UN)
! --- ON RAPPELLE QUE CES INCONNUES SONT DEFINIES AU MILIEU DES COTES
! --- DE TELLE MANIERE QUE LA COMPOSANTE DES ROTATIONS LE LONG DES
! --- COTES EST UNE FONCTION QUADRATIQUE DE L'ABSCISSE CURVILIGNE
!     -----------------------------------------------------------
    call dstcis(dci, carat3, hft2, bca, an)
!
! --- CALCUL DE LA MATRICE B_MEMBRANE NOTEE, ICI, (BM)
!     ------------------------------------------------
    call dxtbm(carat3(9), bm)
!
! --- CALCUL DE LA MATRICE B_FLEXION RELATIVE AUX INCONNUES W, BETAX
! --- ET BETAY ET NE TENANT PAS COMPTE DE L'INTERPOLATION QUADRATIQUE
! --- DES ROTATIONS EN FONCTION DES INCONNUES ALFA.
! --- CETTE MATRICE EST NOTEE (BF_BETA) PAR BATOZ ET ICI (BFB)
!     --------------------------------------------------------
    call dstbfb(carat3(9), bfb)
!
! --- CALCUL DE LA PARTIE DE LA MATRICE B_FLEXION RELATIVE
! --- AUX INCONNUES ALFA
! --- CETTE MATRICE EST NOTEE (BF_ALFA) PAR BATOZ ET ICI (BFA)
!     --------------------------------------------------------
    call dstbfa(qsi, eta, carat3, bfa)
!
! --- CALCUL DE LA MATRICE B_FLEXION COMPLETE, NOTEE (BF) ET RELATIVE
! --- AUX SEULES INCONNUES (UN) :
! --- (BF) = (BFB) + (BFA)*(AN)
!     -------------------------
    do i = 1, 3
        do j = 1, 9
            bfn(i, j) = 0.d0
        end do
    end do
!
    do i = 1, 3
        do j = 1, 9
            do k = 1, 3
                bfn(i, j) = bfn(i, j)+bfa(i, k)*an(k, j)
            end do
            bf(i, j) = bfb(i, j)+bfn(i, j)
        end do
    end do
!
! --- CALCUL DE LA MATRICE B_CISAILLEMENT, NOTEE (BC) ET RELATIVE
! --- AUX  INCONNUES (UN) , AVEC LES NOTATIONS DE BATOZ
! --- (T) = (BCA)*(ALFA) OU T EST LE VECTEUR DES EFFORTS TRANCHANTS.
! --- (T) = (BCA)*(AN)*(UN)
! --- D'AUTRE-PART (GAMMA) = (DCI)*(T) OU GAMMA EST LE VECTEUR
! --- DES DEFORMATIONS DE CISAILLEMENT TRANSVERSE , D'OU :
! --- (BC) = (DCI)*(BCA)*(AN)
!     -----------------------
    do i = 1, 2
        do j = 1, 9
            bc(i, j) = 0.d0
            bcn(i, j) = 0.d0
        end do
    end do
!
    do i = 1, 2
        do j = 1, 9
            do k = 1, 3
                bcn(i, j) = bcn(i, j)+bca(i, k)*an(k, j)
            end do
        end do
    end do
!
    do i = 1, 2
        do j = 1, 9
            do k = 1, 2
                bc(i, j) = bc(i, j)+dci(i, k)*bcn(k, j)
            end do
        end do
    end do
!
! --- AFFECTATION DE LA MATRICE B COMPLETE, NOTEE (BMAT)
! --- AVEC LES MATRICES (BM), (BF) ET (BC)
!     ------------------------------------
    call bcoqaf(bm, bf, bc, nno, bmat)
!
end subroutine
