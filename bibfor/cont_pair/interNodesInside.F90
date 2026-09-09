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
subroutine interNodesInside(proj_tole, elem_dime, &
                            elem_mast_code, elem_slave_code, &
                            proj_coor, nb_node_proj, &
                            nb_poin_inte, poin_inte, inte_neigh, &
                            poin_inte_ori, iret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/apelem_getvertex_n.h"
#include "asterfort/apelem_inside_n.h"
#include "asterfort/ptinma.h"
!
    real(kind=8), intent(in) :: proj_tole
    integer(kind=8), intent(in) :: elem_dime
    character(len=8), intent(in) :: elem_mast_code, elem_slave_code
    real(kind=8), intent(in) :: proj_coor(2, 9)
    integer(kind=8), intent(in) :: nb_node_proj
    integer(kind=8), intent(out) :: inte_neigh(4), nb_poin_inte
    real(kind=8), intent(inout) :: poin_inte(2, 16)
    real(kind=8), intent(inout) :: poin_inte_ori(2, 16)
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Search nodes that are inside cells
!
! --------------------------------------------------------------------------------------------------
!
! In  proj_tole        : tolerance for projection
! In  elem_dime        : dimension of elements
! In  elem_mast_nbnode : number of nodes of master element
! In  elem_slav_nbnode : number of nodes for slave element
! In  elem_slav_coor   : coordinates of slave element
! In  elem_mast_code   : code of master element
! Out proj_coor        : projection of slave nodes on master element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: min_sl, max_sl
    real(kind=8) :: elem_mast_line_coop(2, 4)
    real(kind=8) :: elem_slave_line_coop(2, 4)
    integer(kind=8) :: elem_mast_line_nbnode, i_node, test, elem_slav_line_nbnode
    character(len=8) :: elem_mast_line_code, elem_slav_line_code
    real(kind=8) :: xpt, ypt, cor_inte_ori(2), t1

!
! - Compute intersection in parametric master space
!
    nb_poin_inte = 0
    poin_inte = 0.d0
    inte_neigh = 0
!
    if (elem_dime == 2) then
        ASSERT(nb_node_proj == 2)
        cor_inte_ori(1) = -1.0
        cor_inte_ori(2) = 1.0
! - Add projection of slave nodes inside master cell (parametric space)
        do i_node = 1, nb_node_proj
            xpt = proj_coor(1, i_node)
            if (xpt > (-1.d0-proj_tole) .and. xpt < (1.d0+proj_tole)) then
                nb_poin_inte = nb_poin_inte+1
                poin_inte(1, nb_poin_inte) = xpt
                poin_inte_ori(1, nb_poin_inte) = cor_inte_ori(i_node)
                inte_neigh(i_node) = 1
            end if
        end do
!
! - Add master nodes if they are inside projected slave cell
!
        min_sl = min(proj_coor(1, 1), proj_coor(1, 2))
        max_sl = max(proj_coor(1, 1), proj_coor(1, 2))
!
        if ((min_sl-proj_tole) <= -1.d0 .and. -1.d0 <= (max_sl+proj_tole)) then
            nb_poin_inte = nb_poin_inte+1
            poin_inte(1, nb_poin_inte) = -1.d0
            if (abs(-1.d0-proj_coor(1, 1)) .lt. proj_tole) then
                poin_inte_ori(1, nb_poin_inte) = cor_inte_ori(1)
            else
                t1 = (-1.d0-proj_coor(1, 1))/(proj_coor(1, 2)-proj_coor(1, 1))
                poin_inte_ori(1, nb_poin_inte) = 2*t1-1
            end if
        end if
        if ((min_sl-proj_tole) <= 1.d0 .and. 1.d0 <= (max_sl+proj_tole)) then
            nb_poin_inte = nb_poin_inte+1
            poin_inte(1, nb_poin_inte) = 1.d0
            if (abs(1.d0-proj_coor(1, 1)) .lt. proj_tole) then
                poin_inte_ori(1, nb_poin_inte) = cor_inte_ori(1)
            else
                t1 = (1.d0-proj_coor(1, 1))/(proj_coor(1, 2)-proj_coor(1, 1))
                poin_inte_ori(1, nb_poin_inte) = 2*t1-1
            end if
        end if
!
    elseif (elem_dime == 3) then
!
! - Get parametric coordinates of master nodes (linear)
!
        call apelem_getvertex_n(elem_dime, elem_mast_code, &
                                elem_mast_line_coop, elem_mast_line_nbnode, &
                                elem_mast_line_code)
        call apelem_getvertex_n(elem_dime, elem_slave_code, &
                                elem_slave_line_coop, elem_slav_line_nbnode, &
                                elem_slav_line_code)

!
! - Add projection of slave nodes inside master cell (parametric space)
!
        call apelem_inside_n(proj_tole, elem_dime, elem_mast_line_code, &
                             nb_node_proj, proj_coor, nb_poin_inte, poin_inte, &
                             inte_neigh, elem_slave_line_coop, poin_inte_ori)
!
! - Add master nodes if they are inside projected slave cell
!
        do i_node = 1, elem_mast_line_nbnode
! ----- Current parametric coordinates of master node
            xpt = elem_mast_line_coop(1, i_node)
            ypt = elem_mast_line_coop(2, i_node)
! ----- Test if master node is inside projected slave cell
            call ptinma(elem_slav_line_nbnode, elem_dime, elem_slav_line_code, &
                        proj_coor, proj_tole, xpt, ypt, test, cor_inte_ori)
            if (test == 1) then
                nb_poin_inte = nb_poin_inte+1
                !print*, 'nbptinte', nb_poin_inte
                poin_inte(1:2, nb_poin_inte) = [xpt, ypt]
                poin_inte_ori(1:2, nb_poin_inte) = [cor_inte_ori(1), cor_inte_ori(2)]
                ! print*, "inte1_2_true", poin_inte(1:2,nb_poin_inte)
                ! print*, "inte1_2", poin_inte_ori(1:2,nb_poin_inte)
            else if (test == -1) then
                iret = 1
            end if
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
