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
    subroutine lcmohr(ndim, &
                      typmod, imate, carcri, option, &
                      dstrai0, stresm0, stres, vim, vip, &
                      dsidep, codret)
        integer(kind=8) :: ndim
        character(len=8) :: typmod(*)
        integer(kind=8) :: imate
        character(len=16) :: option
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8) :: dstrai0(6)
        real(kind=8) :: stresm0(6)
        real(kind=8) :: stres(6)
        real(kind=8) :: vim(*)
        real(kind=8) :: vip(*)
        real(kind=8) :: dsidep(6, 6)
        integer(kind=8) :: codret
    end subroutine lcmohr
end interface
