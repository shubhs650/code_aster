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
    subroutine nmvarc_prep(poum, model, caraElem, materCode, varcRefe, &
                           compor, exis_temp, &
                           nbFieldInMax, nbFieldIn, lpain, lchin, &
                           nbFieldOutMax, nbFieldOut, lpaout, lchout, &
                           sigmPrev, variPrev, varcPrev, varcCurr)
        character(len=1), intent(in) :: poum
        character(len=24), intent(in) :: model
        character(len=24), intent(in) :: materCode
        character(len=24), intent(in) :: varcRefe
        character(len=24), intent(in) :: caraElem
        character(len=24), intent(in) :: compor
        aster_logical, intent(in) :: exis_temp
        integer(kind=8), intent(in) :: nbFieldInMax
        character(len=8), intent(inout) :: lpain(nbFieldInMax)
        character(len=19), intent(inout) :: lchin(nbFieldInMax)
        integer(kind=8), intent(out) :: nbFieldIn
        integer(kind=8), intent(in) :: nbFieldOutMax
        character(len=8), intent(inout) :: lpaout(nbFieldOutMax)
        character(len=19), intent(inout) :: lchout(nbFieldOutMax)
        integer(kind=8), intent(out) :: nbFieldOut
        character(len=19), intent(in) :: sigmPrev, variPrev
        character(len=19), intent(in) :: varcPrev, varcCurr
    end subroutine nmvarc_prep
end interface
