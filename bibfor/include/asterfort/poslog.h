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
#include "asterf_types.h"
!
interface
    subroutine poslog(BEHInteg, &
                      lCorr, lMatr, lSigm, lVari, &
                      tlogPrev, tlogCurr, fPrev, &
                      lgpg, vip, ndim, &
                      fCurr, dtde, sigm, cplan, &
                      instp, gn, &
                      lamb, logl, sigmCurr, dsidep, pk2Prev, &
                      pk2Curr, codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        aster_logical, intent(in) :: lCorr, lMatr, lSigm, lVari
        aster_logical, intent(in) :: cplan
        real(kind=8), intent(in) :: tlogPrev(6)
        real(kind=8), intent(in) :: tlogCurr(6)
        real(kind=8), intent(in) :: fPrev(3, 3)
        real(kind=8), intent(in) :: fCurr(3, 3)
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: lgpg
        real(kind=8), intent(out) :: vip(lgpg)
        real(kind=8), intent(in) :: dtde(6, 6)
        real(kind=8), intent(in) :: sigm(2*ndim)
        real(kind=8), intent(in) :: instp
        real(kind=8), intent(in) :: gn(3, 3)
        real(kind=8), intent(in) :: lamb(3)
        real(kind=8), intent(in) :: logl(3)
        real(kind=8), intent(out) :: sigmCurr(2*ndim)
        real(kind=8), intent(out) :: dsidep(6, 6)
        real(kind=8), intent(out) :: pk2Prev(6)
        real(kind=8), intent(out) :: pk2Curr(6)
        integer(kind=8), intent(out) :: codret
    end subroutine poslog
end interface
