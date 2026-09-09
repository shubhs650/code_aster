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
subroutine dxbmat(plateCara, plateOrie, &
                  nomte, cara, xyzl, igau, &
                  jacgau, bmat)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dkqb.h"
#include "asterfort/dktb.h"
#include "asterfort/dsqb.h"
#include "asterfort/dstb.h"
#include "asterfort/q4gb.h"
#include "asterfort/t3gb.h"
#include "asterfort/utmess.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: cara(*), xyzl(3, *)
    integer(kind=8), intent(in) :: igau
    real(kind=8), intent(out) :: jacgau, bmat(8, *)
!
! --------------------------------------------------------------------------------------------------
!
! --- CALCUL DE LA MATRICE (B) RELIANT LES DEFORMATIONS DU PREMIER
! --- ORDRE AUX DEPLACEMENTS AU POINT D'INTEGRATION D'INDICE IGAU
! --- POUR LES ELEMENTS : DST, DKT, DSQ, DKQ, Q4G
! --- (I.E. (EPS_1) = (B)*(UN))
! --- D'AUTRE_PART, ON CALCULE LE PRODUIT NOTE JACGAU = JACOBIEN*POIDS
!
! --------------------------------------------------------------------------------------------------
!
!     IN  NOMTE         : NOM DU TYPE D'ELEMENT
!     IN  XYZL(3,NBNO)  : COORDONNEES DES CONNECTIVITES DE L'ELEMENT
!                         DANS LE REPERE LOCAL DE L'ELEMENT
!     IN  PGL(3,3)      : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                         LOCAL
!     IN  IGAU          : INDICE DU POINT D'INTEGRATION
!     OUT JACGAU        : PRODUIT JACOBIEN*POIDS AU POINT D'INTEGRATION
!                         COURANT
!     OUT BMAT(6,1)     : MATRICE (B) AU POINT D'INTEGRATION COURANT
!
! --------------------------------------------------------------------------------------------------
!
    if (nomte .eq. 'MEDKTR3 ' .or. nomte .eq. 'MEDKTG3 ') then
        call dktb(cara, igau, jacgau, bmat)

    else if (nomte .eq. 'MEDSTR3 ') then
        call dstb(plateCara, plateOrie, &
                  cara, igau, jacgau, bmat)

    else if (nomte .eq. 'MEDKQU4 ' .or. nomte .eq. 'MEDKQG4 ') then
        call dkqb(cara, xyzl, igau, jacgau, bmat)

    else if (nomte .eq. 'MEDSQU4 ') then
        call dsqb(plateCara, plateOrie, &
                  cara, xyzl, igau, jacgau, &
                  bmat)

    else if (nomte .eq. 'MEQ4QU4 ' .or. nomte .eq. 'MEQ4GG4 ') then
        call q4gb(cara, xyzl, igau, jacgau, bmat)

    else if (nomte .eq. 'MET3GG3 ' .or. nomte .eq. 'MET3TR3 ') then
        call t3gb(cara, xyzl, bmat)
        jacgau = cara(8)

    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
