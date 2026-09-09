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
subroutine dxbsig(plateCara, plateOrie, &
                  nomte, optionZ, &
                  xyzl, pgl, sigma, &
                  bsigma)
!
    use plateGeom_module, only: isPlateTria, isPlateQuad, isPlateQ4GG
    use plate_type
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/btsir.h"
#include "asterfort/btsig.h"
#include "asterfort/dxbmat.h"
#include "asterfort/gquad4.h"
#include "asterfort/gtria3.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvlg.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: nomte
    character(len=*), intent(in) :: optionZ
    real(kind=8), intent(in) :: xyzl(3, 1), pgl(3, 3)
    real(kind=8), intent(in) :: sigma(*)
    real(kind=8), intent(out) :: bsigma(*)
!
! --------------------------------------------------------------------------------------------------
!
! --- CALCUL DES FORCES INTERNES B*SIGMA AUX NOEUDS DE L'ELEMENT
! --- DUES AU CHAMP DE CONTRAINTES SIGMA DEFINI AUX POINTS
! --- D'INTEGRATION POUR LES ELEMENTS : DST, DKT, DSQ, DKQ, Q4G
!
! --------------------------------------------------------------------------------------------------
!
!     IN  NOMTE         : NOM DU TYPE D'ELEMENT
!     IN  XYZL(3,NNO)   : COORDONNEES DES CONNECTIVITES DE L'ELEMENT
!                         DANS LE REPERE LOCAL DE L'ELEMENT
!     IN  PGL(3,3)      : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                         LOCAL
!     IN  SIGMA(1)      : CONTRAINTES GENERALISEES DEFINIES AUX POINTS
!                         D'INTEGRATION DE L'ELEMENT
!     OUT BSIGMA(1)     : FORCES INTERNES AUX NOEUDS DE L'ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbsig = 8, lgligb = 24
    integer(kind=8) :: i, kpg, nno, npg
    real(kind=8) :: bsivar
    real(kind=8), parameter :: zero = 0.d0
    character(len=16) :: option
    real(kind=8) :: bmat(nbsig, lgligb)
    real(kind=8) :: bsiloc(lgligb), jacgau, cara(25)
!
! --------------------------------------------------------------------------------------------------
!
    option = optionZ
    bsiloc = zero
    bsigma(1:lgligb) = zero
    bmat = zero
!
    if (isPlateQuad(plateCara)) then
        npg = 4
        nno = 4
        call gquad4(xyzl, cara)
    elseif (isPlateTria(plateCara)) then
        npg = 3
        nno = 3
        if (isPlateQ4GG(plateCara)) then
            npg = 1
        end if
        call gtria3(xyzl, cara)
    else
        ASSERT(ASTER_FALSE)
    end if

! - CALCUL DE SOMME_ELEMENT(BT_SIGMA)
    do kpg = 1, npg
!  --   CALCUL DE LA MATRICE B RELIANT LES DEFORMATIONS DU
!  --   PREMIER ORDRE AUX DEPLACEMENTS AU POINT D'INTEGRATION
!  --   COURANT : (EPS_1) = (B)*(UN)
        call dxbmat(plateCara, plateOrie, &
                    nomte, cara, xyzl, kpg, &
                    jacgau, bmat)

!  --   CALCUL DU PRODUIT (BT)*(SIGMA)*JACOBIEN*POIDS
        if (optionZ .eq. 'FORC_NODA') then
            call btsig(lgligb, nbsig, jacgau, bmat, sigma(1+8*(kpg-1)), &
                       bsiloc)
        elseif (optionZ .eq. 'REFE_FORC_NODA') then
            call btsir(lgligb, nbsig, jacgau, bmat, sigma(1+8*(kpg-1)), &
                       bsiloc)
        else
            ASSERT(ASTER_FALSE)
        end if
    end do

! - PERMUTATION DES COMPOSANTES EN BETA_X ET BETA_Y EN TETA_Y ET -TETA_X
    do i = 1, nno
        bsivar = bsiloc(4+6*(i-1))
        bsiloc(4+6*(i-1)) = -bsiloc(5+6*(i-1))
        bsiloc(5+6*(i-1)) = bsivar
    end do

! - PASSAGE DU VECTEUR(BT_SIGMA) DU REPERE LOCAL AU REPERE GLOBAL :
! - (LE RESULTAT EST ICI LE VECTEUR BSIGMA)
    call utpvlg(nno, 6, pgl, bsiloc, bsigma)
!
end subroutine
