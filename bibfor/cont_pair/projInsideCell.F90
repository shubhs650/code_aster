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

subroutine projInsideCell(pair_tole, elem_dime, elem_code, &
                          poin_coor, iret)
!
    implicit none
!
#include "asterfort/assert.h"
!
!
    real(kind=8), intent(in) :: pair_tole
    integer(kind=8), intent(in) :: elem_dime
    character(len=8), intent(in) :: elem_code
    real(kind=8), intent(in) :: poin_coor(2)
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
! Contact - Pairing segment to segment
!
! Test point is inside parametric cell
!
! --------------------------------------------------------------------------------------------------
!
! In  pair_tole        : tolerance for pairing
! In  elem_dime        : dimension of elements
! In  elem_code        : code of elements
! In  poin_coor        : parametric coordinates of points
! Out iret             : 0- inside / 1- outside
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: xpt, ypt
!
! --------------------------------------------------------------------------------------------------
!
    iret = 1
    if (elem_code .eq. 'SE2' .or. elem_code .eq. 'SE3' &
        .or. elem_code .eq. 'SE4') then
        xpt = poin_coor(1)
        if (xpt .ge. (-1.d0-pair_tole) .and. &
            xpt .le. (1.d0+pair_tole)) then
            iret = 0
        end if
    elseif (elem_code .eq. 'TR3' .or. elem_code .eq. 'TR6' &
            .or. elem_code .eq. 'TR7' .or. elem_code .eq. 'TR1') then
        xpt = poin_coor(1)
        ypt = poin_coor(2)
        if (xpt .ge. -pair_tole .and. &
            ypt .ge. -pair_tole .and. &
            (ypt+xpt) .le. (1.d0+pair_tole)) then
            iret = 0
        end if
    elseif (elem_code .eq. 'QU4' .or. elem_code .eq. 'QU8' &
            .or. elem_code .eq. 'QU9' .or. elem_code .eq. 'Q12') then
        xpt = poin_coor(1)
        ypt = poin_coor(2)
        if (xpt .ge. (-1.d0-pair_tole) .and. &
            xpt .le. (1.d0+pair_tole) .and. &
            ypt .ge. (-1.d0-pair_tole) .and. &
            ypt .le. (1.d0+pair_tole)) then
            iret = 0
        end if
    else
        ASSERT(.false.)
    end if
!
!
end subroutine
