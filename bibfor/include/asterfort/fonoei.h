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
    subroutine fonoei(ds_thm, &
                      ndim, dimdef, dimcon, &
                      addeme, &
                      addep1, addep2, addlh1, adcome, &
                      adcop1, adcop2, congem, &
                      r)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        integer(kind=8) :: dimcon
        integer(kind=8) :: dimdef
        integer(kind=8) :: ndim
        integer(kind=8) :: addeme
        integer(kind=8) :: addep1
        integer(kind=8) :: addep2
        integer(kind=8) :: addlh1
        integer(kind=8) :: adcome
        integer(kind=8) :: adcp11
        integer(kind=8) :: adcop1
        integer(kind=8) :: adcop2
        real(kind=8) :: congem(dimcon)
        real(kind=8) :: r(dimdef)
    end subroutine fonoei
end interface
