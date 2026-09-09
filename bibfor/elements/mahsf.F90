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

subroutine mahsf(plateOrie, &
                 ind1, nb1, &
                 nodeCoor, ksi3s2, intsn, &
                 desr, epais, &
                 vectBaseKpg, vectTangKpg, &
                 hsf)
!
    use plate_type
    implicit none
!
#include "asterfort/hfmss.h"
#include "asterfort/vectgt.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: ind1, nb1
    real(kind=8), intent(in) :: nodeCoor(3, *)
    real(kind=8), intent(in) ::  ksi3s2
    integer(kind=8), intent(in) :: intsn
    real(kind=8), intent(in) :: desr(*), epais
    real(kind=8), intent(out) :: vectBaseKpg(3, 3), vectTangKpg(2, 3)
    real(kind=8), intent(out) :: hsf(3, 9)
!
! --------------------------------------------------------------------------------------------------
!
!     CONSTRUCTION DU VECTEUR N AUX PTS D'INTEGRATION NORMAL
!     (POUR CHAQUE INTSN, STOCKAGE DANS VECTT)
!
!     ET
!
!     CONSTRUCTION DES VECTEURS GA AUX PTS D'INTEGRATION NORMAL
!     (POUR CHAQUE INTSN, STOCKAGE DANS VECTG)
!
!     ET
!
!     CONSTRUCTION DES VECTEURS TA AUX PTS D'INTEGRATION NORMAL (T3=N)
!     (POUR CHAQUE INTSN, STOCKAGE DANS VECTT)
!
!     IND1= 1      1 : CALCULS AUX PTS D'INTEGRATION NORMAL
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ind2 = 0
    real(kind=8) :: hss(2, 9)
!
! --------------------------------------------------------------------------------------------------
!
    call vectgt(plateOrie, ind1, nb1, &
                nodeCoor, ksi3s2, intsn, &
                epais, desr, &
                vectBaseKpg, vectTangKpg)

!   CONSTRUCTION DE HSM = HFM * S :(3,9) AUX PTS D'INTEGRATION NORMAL
    call hfmss(ind2, vectBaseKpg, hsf, hss)
!
end subroutine
