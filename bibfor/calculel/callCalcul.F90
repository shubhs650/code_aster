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
function callCalcul(optionZ)
!
    implicit none
!
#include "asterf_types.h"
!
    character(len=*), intent(in) :: optionZ
    aster_logical :: callCalcul
!
! --------------------------------------------------------------------------------------------------
!
! May I call CALC_CHAMP with this option
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: option
!
! --------------------------------------------------------------------------------------------------
!
    callCalcul = ASTER_TRUE
    option = optionZ

    if ((option .eq. 'ERTH_ELEM') .or. (option .eq. 'ERTH_ELNO') .or. &
        (option .eq. 'ERME_ELEM') .or. (option .eq. 'ERME_ELNO') .or. &
        (option .eq. 'QIRE_ELEM') .or. (option .eq. 'QIRE_ELNO') .or. &
        (option .eq. 'SIZ1_NOEU') .or. (option .eq. 'SIZ2_NOEU') .or. &
        (option .eq. 'ERZ1_ELEM') .or. (option .eq. 'ERZ2_ELEM') .or. &
        (option .eq. 'QIZ1_ELEM') .or. (option .eq. 'QIZ2_ELEM') .or. &
        (option .eq. 'SING_ELEM') .or. (option .eq. 'SING_ELNO')) then
        callCalcul = ASTER_FALSE
    end if
!
end function
