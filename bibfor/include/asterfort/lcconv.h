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
    subroutine lcconv(rela_comp, yd, dy, ddy, &
                      nr, itmax, toler, iter, intg, &
                      nmat, mater, r, rini, epstr, &
                      typess, essai, cutLevel, nvi, &
                      vinf, &
                      iret)
        integer(kind=8) :: nvi
        integer(kind=8) :: nmat
        integer(kind=8) :: nr
        character(len=16), intent(in) :: rela_comp
        real(kind=8) :: yd(*)
        real(kind=8) :: dy(*)
        real(kind=8) :: ddy(*)
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        integer(kind=8) :: iter
        integer(kind=8) :: intg
        real(kind=8) :: mater(nmat, 2)
        real(kind=8) :: r(*)
        real(kind=8) :: rini(*)
        real(kind=8) :: epstr(6)
        integer(kind=8) :: typess
        real(kind=8) :: essai
        integer(kind=8) :: cutLevel
        real(kind=8) :: vinf(nvi)
        real(kind=8) :: vind1(nvi)
        integer(kind=8) :: indi(7)
        aster_logical :: bnews(3)
        aster_logical :: mtrac
        integer(kind=8) :: iret
    end subroutine lcconv
end interface
