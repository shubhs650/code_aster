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
subroutine interNodesEdge(proj_tole, elem_dime, &
                          elem_mast_code, elem_slave_code, &
                          proj_coor, nb_node_proj, &
                          nb_poin_inte, poin_inte, inte_neigh, &
                          poin_inte_ori)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/apelem_getvertex_n.h"
#include "asterfort/gapGetParamCoor.h"
#include "MeshTypes_type.h"
#include "asterfort/insema.h"

!
    real(kind=8), intent(in) :: proj_tole
    integer(kind=8), intent(in) :: elem_dime
    character(len=8), intent(in) :: elem_mast_code, elem_slave_code
    real(kind=8), intent(in) :: proj_coor(2, 9)
    integer(kind=8), intent(in) :: nb_node_proj
    integer(kind=8), intent(inout) :: inte_neigh(4), nb_poin_inte
    real(kind=8), intent(out) :: poin_inte(2, 16)
    real(kind=8), intent(out) :: poin_inte_ori(2, 16)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Intersection of edges in 3D
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
    real(kind=8) :: elem_mast_line_coop(2, 4)

    integer(kind=8) :: elem_mast_line_nbnode, i_node
    integer(kind=8) :: list_next(8), test, i_add, nb_int_add
    character(len=8) :: elem_mast_line_code, elem_code
    real(kind=8) :: xp1, yp1, xp2, yp2, t1, t2, para_coor_ori(2, MT_NNOMAX2D)
    real(kind=8) :: xp1_ori, yp1_ori, xp2_ori, yp2_ori

    t1 = 0.d0
    t2 = 0.d0
!
! - Compute intersection in parametric master space
!
    if (elem_dime == 3) then
        ASSERT(nb_node_proj <= 8)
!
! - Get parametric coordinates of master nodes (linear)
!
        call apelem_getvertex_n(elem_dime, elem_mast_code, &
                                elem_mast_line_coop, elem_mast_line_nbnode, &
                                elem_mast_line_code)
!
! - Set index of next nodes
!
        list_next = 0
        if (elem_slave_code == "TR3") then
            list_next(1:3) = [2, 3, 1]
            elem_code = "TR3"
            call gapGetParamCoor(elem_code, para_coor_ori)
        elseif (elem_slave_code == "TR6") then
            list_next(1:3) = [2, 3, 1]
            elem_code = "TR3"
            call gapGetParamCoor(elem_code, para_coor_ori)
        elseif (elem_slave_code == "QU8" .or. elem_slave_code == "QU9" .or. &
                elem_slave_code == "QU4") then
            list_next(1:4) = [2, 3, 4, 1]
            elem_code = "QU4"
            call gapGetParamCoor(elem_code, para_coor_ori)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! - Intersection of edges
!
        do i_node = 1, nb_node_proj
!
! --------- Segment from edge of projected slave cell
!
            xp1 = proj_coor(1, i_node)
            yp1 = proj_coor(2, i_node)
            xp2 = proj_coor(1, list_next(i_node))
            yp2 = proj_coor(2, list_next(i_node))
            !print*, "x_1",  xp1, "y_1", yp1
            ! print*, "x_2",  xp2, "y_2", yp2
            xp1_ori = para_coor_ori(1, i_node)
            yp1_ori = para_coor_ori(2, i_node)
            xp2_ori = para_coor_ori(1, list_next(i_node))
            yp2_ori = para_coor_ori(2, list_next(i_node))
            !print*, "x_1_ori",  xp1_ori, "y_1_ori", yp1_ori
            !print*, "x_2_ori",  xp2_ori, "y_2_ori", yp2_ori

!
! --------- Compute intersection between edge of master and projected slave cells
!
            test = nb_poin_inte
            call insema(elem_mast_line_nbnode, elem_dime, elem_mast_line_coop, proj_tole, &
                        xp1, yp1, xp2, yp2, nb_poin_inte, poin_inte)
            if (nb_poin_inte .gt. test) then
                inte_neigh(i_node) = 1
                nb_int_add = nb_poin_inte-test
                do i_add = 1, nb_int_add
                    t1 = 0.d0
                    t2 = 0.d0
                    !print*, 'nbptinte', test+i_add
                    !print*, 'test', abs(poin_inte(1,test+i_add)-xp1)
                    if (abs(poin_inte(1, test+i_add)-xp1) .gt. proj_tole) then
                        t1 = (xp2-xp1)/(poin_inte(1, test+i_add)-xp1)
                        poin_inte_ori(1, test+i_add) = (1.0/t1)*(xp2_ori-xp1_ori)+xp1_ori
                    else
                        poin_inte_ori(1, test+i_add) = xp1_ori
                    end if
                    if (abs(poin_inte(2, test+i_add)-yp1) .gt. proj_tole) then
                        t2 = (yp2-yp1)/(poin_inte(2, test+i_add)-yp1)
                        poin_inte_ori(2, test+i_add) = (1.0/t2)*(yp2_ori-yp1_ori)+yp1_ori
                        poin_inte_ori(1, test+i_add) = (1.0/t2)*(xp2_ori-xp1_ori)+xp1_ori
                    else
                        if (t1 .lt. proj_tole) then
                            poin_inte_ori(2, test+i_add) = yp1_ori
                        else
                            poin_inte_ori(2, test+i_add) = (1.0/t1)*(yp2_ori-yp1_ori)+yp1_ori
                        end if
                    end if

                    !print*, "inte2_1_true", poin_inte(1:2,test+i_add)
                    !print*, "inte2_1", poin_inte_ori(1:2,test+i_add)
                    !print*, 't1', t1
                    !print*, 't2', t2
                end do
            end if
        end do
    end if
!
end subroutine
