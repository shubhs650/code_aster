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
subroutine dsqb(plateCara, plateOrie, &
                caraq4, xyzl, igau, jacgau, &
                bmat)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/bcoqaf.h"
#include "asterfort/dsqbfa.h"
#include "asterfort/dsqbfb.h"
#include "asterfort/dsqcis.h"
#include "asterfort/dsqdis.h"
#include "asterfort/dsxhft.h"
#include "asterfort/dxhmft.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jquad4.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8), intent(in) :: caraq4(*), xyzl(3, 1)
    integer(kind=8), intent(in) :: igau
    real(kind=8), intent(out) :: bmat(8, 1), jacgau
! --- CALCUL DE LA MATRICE (B) RELIANT LES DEFORMATIONS DU PREMIER
! --- ORDRE AUX DEPLACEMENTS AU POINT D'INTEGRATION D'INDICE IGAU
! --- POUR UN ELEMENT DE TYPE DSQ
! --- (I.E. (EPS_1) = (B)*(UN))
! --- D'AUTRE_PART, ON CALCULE LE PRODUIT NOTE JACGAU = JACOBIEN*POIDS
!     ------------------------------------------------------------------
!     IN  XYZL(3,NNO)   : COORDONNEES DES CONNECTIVITES DE L'ELEMENT
!                         DANS LE REPERE LOCAL DE L'ELEMENT
!     IN  PGL(3,3)      : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                         LOCAL
!     IN  IGAU          : INDICE DU POINT D'INTEGRATION
!     OUT JACGAU        : PRODUIT JACOBIEN*POIDS AU POINT D'INTEGRATION
!                         COURANT
!     OUT BMAT(8,1)     : MATRICE (B) AU POINT D'INTEGRATION COURANT
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bfb(3, 12), bfa(3, 4), bfn(3, 12), bf(3, 12)
    real(kind=8) :: bcb(2, 12), bca(2, 4), bcn(2, 12), bc(2, 12), bcm(2, 8)
    real(kind=8) :: hft2(2, 6), an(4, 12), hmft2(2, 6)
    real(kind=8) :: bm(3, 8), qsi, eta, jacob(5)
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: multic, i, j, k
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
! --- CALCUL DES MATRICES DE HOOKE DE FLEXION, MEMBRANE,
! --- MEMBRANE-FLEXION, CISAILLEMENT, CISAILLEMENT INVERSE
!     ----------------------------------------------------
    call dxmate(plateCara, plateOrie, &
                'RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)
!
! --- CALCUL DE LA MATRICE (AN) RELIANT LES INCONNUES NOTEES (ALFA)
! --- PAR BATOZ AUX INCONNUES (UN) = (...,W_I,BETAX_I,BETAY_I,...)
! --- SOIT (ALFA) = (AN)*(UN)
! --- ON RAPPELLE QUE CES INCONNUES SONT DEFINIES AU MILIEU DES COTES
! --- DE TELLE MANIERE QUE LA COMPOSANTE DES ROTATIONS LE LONG DES
! --- COTES EST UNE FONCTION QUADRATIQUE DE L'ABSCISSE CURVILIGNE
!     -----------------------------------------------------------
    call dsqdis(xyzl, caraq4, df, dci, an)
!
! --- CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT AU POINT
! --- D'INTEGRATION D'INDICE IGAU
!     ---------------------------
    call jquad4(xyzl, qsi, eta, jacob)
!
! --- PRODUIT JACOBIEN*POIDS
!     ----------------------
    jacgau = zr(ipoids+igau-1)*jacob(1)
!
! --- CALCUL DE LA MATRICE B_MEMBRANE NOTEE, ICI, (BM)
!     ------------------------------------------------
    call dxqbm(qsi, eta, jacob(2), bm)
!
! --- CALCUL DE LA MATRICE B_FLEXION RELATIVE AUX INCONNUES W, BETAX
! --- ET BETAY ET NE TENANT PAS COMPTE DE L'INTERPOLATION QUADRATIQUE
! --- DES ROTATIONS EN FONCTION DES INCONNUES ALFA.
! --- CETTE MATRICE EST NOTEE (BF_BETA) PAR BATOZ ET ICI (BFB)
!     --------------------------------------------------------
    call dsqbfb(qsi, eta, jacob(2), bfb)
!
! --- CALCUL DE LA PARTIE DE LA MATRICE B_FLEXION RELATIVE
! --- AUX INCONNUES ALFA
! --- CETTE MATRICE EST NOTEE (BF_ALFA) PAR BATOZ ET ICI (BFA)
!     --------------------------------------------------------
    call dsqbfa(qsi, eta, jacob(2), caraq4, bfa)
!
! --- CALCUL DE LA MATRICE B_FLEXION COMPLETE, NOTEE (BF) ET RELATIVE
! --- AUX SEULES INCONNUES (UN) :
! --- (BF) = (BFB) + (BFA)*(AN)
!     -------------------------
    do i = 1, 3
        do j = 1, 12
            bfn(i, j) = 0.d0
        end do
    end do
!
    do i = 1, 3
        do j = 1, 12
            do k = 1, 4
                bfn(i, j) = bfn(i, j)+bfa(i, k)*an(k, j)
            end do
            bf(i, j) = bfb(i, j)+bfn(i, j)
        end do
    end do
!
! --- CALCUL DE LA MATRICE B_CISAILLEMENT, NOTEE (BC) ET RELATIVE
! --- AUX  INCONNUES (UN) , AVEC LES NOTATIONS DE BATOZ
! --- (T) = (BCB)*(UN) + (BCA)*(ALFA) OU T EST LE VECTEUR DES EFFORTS
! --- TRANCHANTS. SOIT (T) = ((BCB)+(BCA)*(AN))*(UN)
! --- D'AUTRE-PART (GAMMA) = (DCI)*(T) OU GAMMA EST LE VECTEUR
! --- DES DEFORMATIONS DE CISAILLEMENT TRANSVERSE , D'OU :
! --- (BC) = (DCI)*((BCB)+(BCA)*(AN))
!     -------------------------------
!
! --- CALCUL DU PRODUIT (HF)*(T2)
!     ---------------------------
    call dsxhft(df, jacob(2), hft2)
!
! --- CALCUL DU PRODUIT HMF.T2 :
!     ------------------------
    call dxhmft(dmf, jacob(2), hmft2)
!
! --- CALCUL DES MATRICES (BCB) ET (BCA) ET (BCM) :
!     -------------------------------------------
    call dsqcis(qsi, eta, caraq4, hmft2, hft2, &
                bcm, bcb, bca)
!
    do i = 1, 2
        do j = 1, 12
            bc(i, j) = 0.d0
            bcn(i, j) = 0.d0
        end do
    end do
!
    do i = 1, 2
        do j = 1, 12
            do k = 1, 4
                bcn(i, j) = bcn(i, j)+bca(i, k)*an(k, j)
            end do
        end do
    end do
!
    do i = 1, 2
        do j = 1, 12
            do k = 1, 2
                bc(i, j) = bc(i, j)+dci(i, k)*(bcb(k, j)+bcn(k, j))
            end do
        end do
    end do
!
!
! --- AFFECTATION DE LA MATRICE B COMPLETE, NOTEE (BMAT)
! --- AVEC LES MATRICES (BM), (BF) ET (BC)
!     ------------------------------------
    call bcoqaf(bm, bf, bc, nno, bmat)
!
end subroutine
