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
    subroutine compEnergyPotential(optionZ, &
                                   modelZ, materCodeZ, caraElemZ, comporZ, &
                                   chdispZ, chharmZ, chgeomZ, &
                                   chtimeZ, chvarcZ, chvrefZ, &
                                   l_temp, chtempZ, &
                                   ligrelZ, jvBaseZ, epotElemZ, &
                                   codret)
        character(len=*), intent(in) :: optionZ
        character(len=*), intent(in) :: modelZ, materCodeZ, caraElemZ, comporZ
        character(len=*), intent(in) :: chdispZ, chharmZ, chgeomZ, chtimeZ
        character(len=*), intent(in) :: chvarcZ, chvrefZ
        aster_logical, intent(in) :: l_temp
        character(len=*), intent(in) :: chtempZ
        character(len=*), intent(in) :: ligrelZ, jvBaseZ, epotElemZ
        integer(kind=8), intent(out) :: codret
    end subroutine compEnergyPotential
end interface
