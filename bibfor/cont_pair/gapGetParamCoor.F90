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
subroutine gapGetParamCoor(elem_code, para_coor)
!
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
!
    character(len=8), intent(in) :: elem_code
    real(kind=8), intent(out) :: para_coor(2, MT_NNOMAX2D)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Gap computation
!
! Get vertices of element in parametric space
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_dime        : dimension of current element
! In  elem_code        : code of current element
! Out para_coor        : parametric coordinates for vertices of current element
!
! --------------------------------------------------------------------------------------------------
!

    para_coor(:, :) = 0.d0
    if (elem_code .eq. 'SE3') then
        para_coor(1, 1) = -1.d0
        para_coor(1, 2) = 1.d0
        para_coor(1, 3) = 0.d0
    elseif (elem_code .eq. 'SE2') then
        para_coor(1, 1) = -1.d0
        para_coor(1, 2) = 1.d0
    elseif (elem_code .eq. 'SE4') then
        para_coor(1, 1) = -1.d0
        para_coor(1, 2) = 1.d0
        para_coor(1, 3) = -1.d0/3.d0
        para_coor(1, 4) = 1.d0/3.d0
    elseif (elem_code .eq. 'TR3') then
        para_coor(1, 1) = 0.d0
        para_coor(2, 1) = 0.d0
        para_coor(1, 2) = 1.d0
        para_coor(2, 2) = 0.d0
        para_coor(1, 3) = 0.d0
        para_coor(2, 3) = 1.d0
    elseif (elem_code .eq. 'TR6') then
        para_coor(1, 1) = 0.d0
        para_coor(2, 1) = 0.d0
        para_coor(1, 2) = 1.d0
        para_coor(2, 2) = 0.d0
        para_coor(1, 3) = 0.d0
        para_coor(2, 3) = 1.d0
        para_coor(1, 4) = 0.5d0
        para_coor(2, 4) = 0.d0
        para_coor(1, 5) = 0.5d0
        para_coor(2, 5) = 0.5d0
        para_coor(1, 6) = 0.d0
        para_coor(2, 6) = 0.5d0
    else if (elem_code .eq. 'QU4') then
        para_coor(1, 1) = -1.0d0
        para_coor(2, 1) = -1.0d0
        para_coor(1, 2) = 1.d0
        para_coor(2, 2) = -1.d0
        para_coor(1, 3) = 1.d0
        para_coor(2, 3) = 1.d0
        para_coor(1, 4) = -1.d0
        para_coor(2, 4) = 1.d0
    else if (elem_code .eq. 'QU8') then
        para_coor(1, 1) = -1.0d0
        para_coor(2, 1) = -1.0d0
        para_coor(1, 2) = 1.d0
        para_coor(2, 2) = -1.d0
        para_coor(1, 3) = 1.d0
        para_coor(2, 3) = 1.d0
        para_coor(1, 4) = -1.d0
        para_coor(2, 4) = 1.d0
        para_coor(1, 5) = 0.d0
        para_coor(2, 5) = -1.d0
        para_coor(1, 6) = 1.d0
        para_coor(2, 6) = 0.d0
        para_coor(1, 7) = 0.d0
        para_coor(2, 7) = 1.d0
        para_coor(1, 8) = -1.d0
        para_coor(2, 8) = 0.d0
    else if (elem_code .eq. 'QU9') then
        para_coor(1, 1) = -1.0d0
        para_coor(2, 1) = -1.0d0
        para_coor(1, 2) = 1.d0
        para_coor(2, 2) = -1.d0
        para_coor(1, 3) = 1.d0
        para_coor(2, 3) = 1.d0
        para_coor(1, 4) = -1.d0
        para_coor(2, 4) = 1.d0
        para_coor(1, 5) = 0.d0
        para_coor(2, 5) = -1.d0
        para_coor(1, 6) = 1.d0
        para_coor(2, 6) = 0.d0
        para_coor(1, 7) = 0.d0
        para_coor(2, 7) = 1.d0
        para_coor(1, 8) = -1.d0
        para_coor(2, 8) = 0.d0
        para_coor(1, 9) = 0.d0
        para_coor(2, 9) = 0.d0
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
