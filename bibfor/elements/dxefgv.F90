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
subroutine dxefgv(plateCara, plateOrie, &
                  nomte, optionZ, xyzl, pgl, depl, effgt)
!
    use plate_type
    implicit none
!
#include "asterfort/dxefg2.h"
#include "asterfort/dxefgm.h"
#include "asterfort/dxefgt.h"
#include "asterfort/dxefn2.h"
#include "asterfort/dxefnt.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16) :: nomte
    character(len=*) :: optionZ
    real(kind=8) :: xyzl(3, 1), pgl(3, 3)
    real(kind=8) :: depl(*)
    real(kind=8) :: effgt(*)
!
! --------------------------------------------------------------------------------------------------
!
! EFFORTS GENERALISES 'VRAIS' (I.E. EFFOR_MECA - EFFOR_THERM)
! AUX POINTS D'INTEGRATION POUR LES ELEMENTS COQUES A
! FACETTES PLANES : DST, DKT, DSQ, DKQ, Q4G
!
! --------------------------------------------------------------------------------------------------
!
!     IN  NOMTE        : NOM DU TYPE D'ELEMENT
!     IN  OPTION       : NOM DE L'OPTION
!     IN  XYZL(3,NNO)  : COORDONNEES DES CONNECTIVITES DE L'ELEMENT
!                        DANS LE REPERE LOCAL DE L'ELEMENT
!     IN  PGL(3,3)     : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                        LOCAL
!     IN  DEPL(1)      : VECTEUR DES DEPLACEMENTS AUX NOEUDS
!     IN  TSUP(1)      : TEMPERATURES AUX NOEUDS DU PLAN SUPERIEUR
!                        DE LA COQUE
!     IN  TINF(1)      : TEMPERATURES AUX NOEUDS DU PLAN INFERIEUR
!                        DE LA COQUE
!     IN  TMOY(1)      : TEMPERATURES AUX NOEUDS DU PLAN MOYEN
!                        DE LA COQUE
!     OUT EFFGT(1)     : EFFORTS  GENERALISES VRAIS AUX POINTS
!                        D'INTEGRATION (I.E.
!                           EFFORTS_MECA - EFFORTS_THERM)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i
    character(len=16) :: option
    real(kind=8) :: sigmTher(32)
!
! --------------------------------------------------------------------------------------------------
!
    option = optionZ

! - CALCUL DES EFFORTS GENERALISES
    call dxefgm(plateCara, plateOrie, &
                nomte, option, xyzl, depl, effgt)

! - CALCUL DES EFFORTS GENERALISES D'ORIGINE THERMIQUE
    if (option(8:9) .eq. 'GA') then
        if (nomte .eq. 'MEDKQG4' .or. nomte .eq. 'MEDKTG3') then
            call dxefg2(plateCara, plateOrie, &
                        pgl, sigmTher)
        else
            call dxefgt(plateCara, plateOrie, &
                        sigmTher)
        end if
    else if (option(8:9) .eq. 'NO') then
        if (nomte .eq. 'MEDKQG4' .or. nomte .eq. 'MEDKTG3') then
            call dxefn2(plateCara, plateOrie, &
                        nomte, pgl, sigmTher)
        else
            call dxefnt(plateCara, plateOrie, &
                        sigmTher)
        end if
    end if

! - CALCUL DES EFFORTS GENERALISES 'VRAIS'
    do i = 1, 32
        effgt(i) = effgt(i)-sigmTher(i)
    end do
!
end subroutine
