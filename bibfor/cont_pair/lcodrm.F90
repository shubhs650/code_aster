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

subroutine lcodrm(elem_dime, pair_tole, nb_poin_inte, poin_inte, poin_inte2)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/ordr8.h"

! aslint: disable=W1306
!
    integer(kind=8), intent(in) :: elem_dime
    real(kind=8), intent(in) :: pair_tole
    integer(kind=8), intent(inout) :: nb_poin_inte
    real(kind=8), intent(inout) :: poin_inte(2, nb_poin_inte)
    real(kind=8), optional, intent(inout) :: poin_inte2(2, nb_poin_inte)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Suppression des doublons et classement des noeuds d'intersections
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_dime        : dimension of current element
! In  pair_tole        : tolerance for pairing
! IO  nb_poin_inte     : number of intersection points
! IO  poin_inte        : list of intersection points
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: poin_inte_sort(2, 16)
    real(kind=8) :: poin_inte_sort2(2, 16)
    real(kind=8) :: angle(nb_poin_inte), bary(2)
    real(kind=8) :: v(2), norm
    integer(kind=8) :: i_poin_inte, angle_sorted(nb_poin_inte), list_poin_next(nb_poin_inte)
    integer(kind=8) :: nb_inte_new
!
! --------------------------------------------------------------------------------------------------
!
    bary(1:2) = 0.d0
    v(1:2) = 0.d0
    nb_inte_new = 0
    if ((elem_dime-1) .eq. 2) then
!
! ----- Coordinates of barycenter
!
        do i_poin_inte = 1, nb_poin_inte
            bary(1) = bary(1)+poin_inte(1, i_poin_inte)/real(nb_poin_inte)
            bary(2) = bary(2)+poin_inte(2, i_poin_inte)/real(nb_poin_inte)
        end do
!
! ----- Compute angles
!
        do i_poin_inte = 1, nb_poin_inte
            v(1) = poin_inte(1, i_poin_inte)-bary(1)
            v(2) = poin_inte(2, i_poin_inte)-bary(2)
            angle(i_poin_inte) = atan2(v(1), v(2))
        end do
!
! ----- Sort angles
!
        call ordr8(angle, nb_poin_inte, angle_sorted)
!
! ----- Set index of next points
!
        do i_poin_inte = 2, nb_poin_inte
            list_poin_next(i_poin_inte-1) = i_poin_inte
        end do
        list_poin_next(nb_poin_inte) = 1
!
! ----- Sort
!
        nb_inte_new = 0
        do i_poin_inte = 1, nb_poin_inte
            norm = sqrt((angle(angle_sorted(i_poin_inte))- &
                         angle(angle_sorted(list_poin_next(i_poin_inte))))**2)
            if (norm .gt. 10*pair_tole) then
                nb_inte_new = nb_inte_new+1
                poin_inte_sort(1, nb_inte_new) = poin_inte(1, angle_sorted(i_poin_inte))
                poin_inte_sort(2, nb_inte_new) = poin_inte(2, angle_sorted(i_poin_inte))
                if (present(poin_inte2)) then
                    poin_inte_sort2(1, nb_inte_new) = poin_inte2(1, angle_sorted(i_poin_inte))
                    poin_inte_sort2(2, nb_inte_new) = poin_inte2(2, angle_sorted(i_poin_inte))
                end if
            end if
        end do
        nb_poin_inte = nb_inte_new
    elseif ((elem_dime-1) .eq. 1) then
        poin_inte_sort(1, 1) = poin_inte(1, 1)
        poin_inte_sort(1, 2) = poin_inte(1, 1)
        if (present(poin_inte2)) then
            poin_inte_sort2(1, 1) = poin_inte2(1, 1)
            poin_inte_sort2(1, 2) = poin_inte2(1, 1)
        end if
        do i_poin_inte = 2, nb_poin_inte
            if (poin_inte(1, i_poin_inte) .le. poin_inte_sort(1, 1) .and. &
                poin_inte(1, i_poin_inte) .ge. (-1.d0)) then
                poin_inte_sort(1, 1) = poin_inte(1, i_poin_inte)
                if (present(poin_inte2)) then
                    poin_inte_sort2(1, 1) = poin_inte2(1, i_poin_inte)
                end if
            elseif (poin_inte(1, i_poin_inte) .ge. poin_inte_sort(1, 2) .and. &
                    poin_inte(1, i_poin_inte) .le. (1.d0)) then
                poin_inte_sort(1, 2) = poin_inte(1, i_poin_inte)
                if (present(poin_inte2)) then
                    poin_inte_sort2(1, 2) = poin_inte2(1, i_poin_inte)
                end if
            end if
        end do
        nb_poin_inte = 2
    else
        ASSERT(.false.)
    end if
!
! - Copy
!
    poin_inte(:, :) = 0.d0
    if (present(poin_inte2)) then
        poin_inte2(:, :) = 0.d0
    end if
    do i_poin_inte = 1, nb_poin_inte
        poin_inte(1, i_poin_inte) = poin_inte_sort(1, i_poin_inte)
        if (present(poin_inte2)) then
            poin_inte2(1, i_poin_inte) = poin_inte_sort2(1, i_poin_inte)
        end if
        if (elem_dime == 3) then
            poin_inte(2, i_poin_inte) = poin_inte_sort(2, i_poin_inte)
            if (present(poin_inte2)) then
                poin_inte2(2, i_poin_inte) = poin_inte_sort2(2, i_poin_inte)
            end if
        end if
    end do
!
end subroutine
