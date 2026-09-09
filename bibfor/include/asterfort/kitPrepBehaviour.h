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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine kitPrepBehaviour(compor, nvi_tot, comporFlua, comporPlas)
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        integer(kind=8), intent(in) :: nvi_tot
        character(len=16), intent(out) :: comporFlua(COMPOR_SIZE)
        character(len=16), intent(out) :: comporPlas(COMPOR_SIZE)
    end subroutine kitPrepBehaviour
end interface
