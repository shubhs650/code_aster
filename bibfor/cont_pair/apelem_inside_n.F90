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

subroutine apelem_inside_n(pair_tole, elem_dime, elem_code, &
                           nb_poin_coor, poin_coor, &
                           nb_poin_inte, poin_inte, inte_neigh, &
                           poin_coor_ori, poin_inte_ori)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/projInsideCell.h"
!
!
    real(kind=8), intent(in) :: pair_tole
    integer(kind=8), intent(in) :: elem_dime
    character(len=8), intent(in) :: elem_code
    integer(kind=8), intent(in) :: nb_poin_coor
    real(kind=8), intent(in) :: poin_coor(2, 4)
    real(kind=8), intent(in) :: poin_coor_ori(2, 4)
    integer(kind=8), intent(inout) :: nb_poin_inte
    real(kind=8), intent(inout) :: poin_inte(2, 16)
    real(kind=8), intent(inout) :: poin_inte_ori(2, 16)
    integer(kind=8), intent(inout) :: inte_neigh(4)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Add points in list of intersection point if inside element
!
! --------------------------------------------------------------------------------------------------
!
! In  pair_tole        : tolerance for pairing
! In  elem_dime        : dimension of elements
! In  elem_code        : code of element
! In  nb_poin_coor     : number of points
! In  poin_coor        : parametric coordinates of points
! IO  nb_poin_inte     : number of intersection points
! IO  poin_inte        : list of intersection points
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_node, iret
    integer(kind=8) :: prev(nb_poin_coor)
!
! --------------------------------------------------------------------------------------------------
!
    do i_node = 2, nb_poin_coor
        prev(i_node) = i_node-1
    end do
    prev(1) = nb_poin_coor
    if (elem_dime .eq. 2) then
        do i_node = 1, nb_poin_coor
            call projInsideCell(pair_tole, elem_dime, elem_code, poin_coor(1:2, i_node), iret)
            if (iret == 0) then
                nb_poin_inte = nb_poin_inte+1
                ASSERT(nb_poin_inte .le. 16)
                poin_inte(1, nb_poin_inte) = poin_coor(1, i_node)
                poin_inte_ori(1, nb_poin_inte) = poin_coor_ori(1, i_node)
                inte_neigh(i_node) = 1
            end if
        end do
    elseif (elem_dime .eq. 3) then
        do i_node = 1, nb_poin_coor
            call projInsideCell(pair_tole, elem_dime, elem_code, poin_coor(1:2, i_node), iret)
            if (iret == 0) then
                nb_poin_inte = nb_poin_inte+1
                ASSERT(nb_poin_inte .le. 16)
                poin_inte(1:2, nb_poin_inte) = poin_coor(1:2, i_node)
                poin_inte_ori(1:2, nb_poin_inte) = poin_coor_ori(1:2, i_node)
                !print*, 'nbptinte', nb_poin_inte
                !print*, "inte1_1_true", poin_inte(1:2,nb_poin_inte)
                ! print*, "inte1_1", poin_inte_ori(1:2,nb_poin_inte)
                inte_neigh(i_node) = 1
                inte_neigh(prev(i_node)) = 1
            end if
        end do
    else
        ASSERT(.false.)
    end if
!
!
end subroutine
