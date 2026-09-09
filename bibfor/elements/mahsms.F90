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

subroutine mahsms(plateOrie, &
                  ind1, nb1, &
                  nodeCoor, ksi3s2, intsr, &
                  desr, epais, &
                  vectBaseKpg, vectTangKpg, &
                  hsfm, hss)
!
    use plate_type
    implicit none
!
#include "asterfort/hfmss.h"
#include "asterfort/vectgt.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: ind1, nb1
    real(kind=8), intent(in) :: nodeCoor(3, *), ksi3s2
    integer(kind=8), intent(in) :: intsr
    real(kind=8), intent(in) :: desr(*), epais
    real(kind=8), intent(out) :: vectBaseKpg(3, 3), vectTangKpg(2, 3)
    real(kind=8), intent(out) :: hsfm(3, 9), hss(2, 9)
!
! --------------------------------------------------------------------------------------------------
!
!     CONSTRUCTION DU VECTEUR N AUX PTS D'INTEGRATION REDUIT
!     (POUR CHAQUE INTSR, STOCKAGE DANS VECTT)
!
!     ET
!
!     CONSTRUCTION DES VECTEURS GA AUX PTS D'INTEGRATION REDUIT
!     (POUR CHAQUE INTSR, STOCKAGE DANS VECTG)
!
!     ET
!
!     CONSTRUCTION DES VECTEURS TA AUX PTS D'INTEGRATION REDUIT (T3=N)
!     (POUR CHAQUE INTSR, STOCKAGE DANS VECTT)
!
!     IND1= 0     0 : CALCULS AUX PTS D'INTEGRATION REDUIT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ind2
!
! --------------------------------------------------------------------------------------------------
!
    call vectgt(plateOrie, ind1, nb1, &
                nodeCoor, ksi3s2, intsr, &
                epais, desr, &
                vectBaseKpg, vectTangKpg)
!
!     CONSTRUCTION DE HSM = HFM * S:(3,9) AUX PTS D'INTEGRATION REDUITS
!
!     ET
!
!     CONSTRUCTION DE HSS = HS * S :(2,9) AUX PTS D'INTEGRATION REDUITS
!
!     IND2= 1  --->  CALCUL DE HSS ( 0 SINON )
!
    ind2 = 1
!
    call hfmss(ind2, vectBaseKpg, hsfm, hss)
!
end subroutine
