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

subroutine ptinma(elem_nbnode, elem_dime, elem_code, elem_coor, pair_tole, &
                  poin_coorx, poin_coory, test, cor_inte_ori)
!
    implicit none
!
#include "asterfort/assert.h"
!
    integer(kind=8), intent(in) :: elem_nbnode
    integer(kind=8), intent(in) :: elem_dime
    character(len=8), intent(in) :: elem_code
    real(kind=8), intent(in) :: elem_coor(2, elem_nbnode)
    real(kind=8), intent(in) :: pair_tole
    real(kind=8), intent(in) :: poin_coorx
    real(kind=8), intent(in) :: poin_coory
    integer(kind=8), intent(out) :: test
    real(kind=8), intent(out), optional ::cor_inte_ori(2)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Test if point is inside element
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_nbnode      : number of nodes of element
! In  elem_dime        : dimension of current element
! In  elem_code        : code of current element
! In  elem_coor        : coordinates of nodes for current element
! In  pair_tole        : tolerance for pairing
! In  poin_coorx       : x coordinate of point
! In  poin_coory       : y coordinate of point
! Out test             : flag for position of point about element
!                        -1 - Error (aligned points)
!                         0 - Point not in element
!                        +1 - Point is in element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: v0(2), v1(2), v2(2), d00, d10, d11, m, u, v
    real(kind=8) :: d02, d12, xpmin, xpmax
!
! --------------------------------------------------------------------------------------------------
!
    test = -1
    ASSERT(elem_code .eq. 'TR3' .or. elem_code .eq. 'SE2' .or. elem_code .eq. 'QU4')
    if (present(cor_inte_ori)) then
        cor_inte_ori = 0.d0
    end if
    if (elem_dime .eq. 3) then
!
! ----- Vectorial basis for element
!
        v0(1) = elem_coor(1, 2)-elem_coor(1, 1)
        v0(2) = elem_coor(2, 2)-elem_coor(2, 1)
        v1(1) = elem_coor(1, 3)-elem_coor(1, 1)
        v1(2) = elem_coor(2, 3)-elem_coor(2, 1)
        d00 = v0(1)*v0(1)+v0(2)*v0(2)
        d10 = v0(1)*v1(1)+v0(2)*v1(2)
        d11 = v1(1)*v1(1)+v1(2)*v1(2)
        m = (d00*d11-d10*d10)
!
! ----- Degenerated vectorial basis for element (colinear vectors) => exit
!
        if (abs(m) .le. pair_tole) then
            test = -1
            goto 99
        end if
!
! ----- Coordinates of point in element's basis
!
        v2(1) = poin_coorx-elem_coor(1, 1)
        v2(2) = poin_coory-elem_coor(2, 1)
        d02 = v0(1)*v2(1)+v0(2)*v2(2)
        d12 = v1(1)*v2(1)+v1(2)*v2(2)
!
! ----- Point is in element => exit
!
        if (sqrt(v2(1)**2+v2(2)**2) .le. 0.d0+pair_tole) then
            test = 1
            if (present(cor_inte_ori)) then
                if (elem_code .eq. 'TR3') then
                    cor_inte_ori(1) = 0.d0
                    cor_inte_ori(2) = 0.d0
                elseif (elem_code .eq. 'QU4') then
                    cor_inte_ori(1) = -1.d0
                    cor_inte_ori(2) = -1.d0
                end if
            end if
            goto 99
        end if
!
! ----- Extension with pair_tole
!
        u = 1/m*(d11*d02-d10*d12)
        v = 1/m*(d00*d12-d10*d02)
        if (u .ge. (0.d0-pair_tole) .and. &
            v .ge. (0.d0-pair_tole) .and. &
            (u+v) .le. (1.d0+pair_tole)) then
            test = 1
            if (present(cor_inte_ori)) then
                if (elem_code .eq. 'TR3') then
                    cor_inte_ori(1) = 0.d0+u
                    cor_inte_ori(2) = 0.d0+v
                elseif (elem_code .eq. 'QU4') then
                    cor_inte_ori(1) = -1.d0+u*2+v*2
                    cor_inte_ori(2) = -1.d0+v*2
                end if
            end if
            goto 99
        else
            test = 0
        end if
!
! ----- Coordinates for test next(QUAD4)
!
        if (elem_code .eq. 'QU4') then
!
! --------- Vectorial basis for element
!
            v0(1) = elem_coor(1, 3)-elem_coor(1, 1)
            v0(2) = elem_coor(2, 3)-elem_coor(2, 1)
            v1(1) = elem_coor(1, 4)-elem_coor(1, 1)
            v1(2) = elem_coor(2, 4)-elem_coor(2, 1)
            d00 = v0(1)*v0(1)+v0(2)*v0(2)
            d10 = v0(1)*v1(1)+v0(2)*v1(2)
            d11 = v1(1)*v1(1)+v1(2)*v1(2)
            m = (d00*d11-d10*d10)
!
! --------- Degenerated vectorial basis for element (colinear vectors) => exit
!
            if (abs(m) .le. pair_tole) then
                test = -1
                goto 99
            end if
!
! --------- Coordinates of point in element's basis
!
            v2(1) = poin_coorx-elem_coor(1, 1)
            v2(2) = poin_coory-elem_coor(2, 1)
            d02 = v0(1)*v2(1)+v0(2)*v2(2)
            d12 = v1(1)*v2(1)+v1(2)*v2(2)
!
! --------- Point is in element => exit
!
            if (sqrt(v2(1)**2+v2(2)**2) .le. 0.d0+pair_tole) then
                test = 1
                if (present(cor_inte_ori)) then
                    if (elem_code .eq. 'QU4') then
                        cor_inte_ori(1) = -1.d0
                        cor_inte_ori(2) = -1.d0
                    end if
                end if
                goto 99
            end if
!
! --------- Extension with pair_tole
!
            u = 1/m*(d11*d02-d10*d12)
            v = 1/m*(d00*d12-d10*d02)
            if (u .ge. (0.d0-pair_tole) .and. &
                v .ge. (0.d0-pair_tole) .and. &
                (u+v) .le. (1.d0+pair_tole)) then
                test = 1
                if (present(cor_inte_ori)) then
                    if (elem_code .eq. 'QU4') then
                        cor_inte_ori(1) = -1.d0+u*2
                        cor_inte_ori(2) = -1.d0+v*2+u*2
                    end if
                end if
                goto 99
            else
                test = 0
            end if
        end if
    elseif (elem_dime .eq. 2) then
        xpmin = min(elem_coor(1, 1), elem_coor(1, 2))
        xpmax = max(elem_coor(1, 1), elem_coor(1, 2))
        if (poin_coorx .ge. (xpmin-pair_tole) .and. &
            poin_coorx .le. (xpmax+pair_tole)) then
            test = 1
            if (present(cor_inte_ori)) then
                cor_inte_ori(1) = 2.0*((poin_coorx-elem_coor(1, 1))/ &
                                       (elem_coor(1, 2)-elem_coor(1, 1)))-1.d0
            end if
        else
            test = 0
        end if
    else
        ASSERT(.false.)
    end if
99  continue

end subroutine
