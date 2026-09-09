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
subroutine prjint_ray(proj_tole, dist_ratio, elem_dime, &
                      elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                      elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                      poin_inte_ma, poin_inte_es, inte_weight, nb_poin_inte, &
                      inte_neigh_, ierror_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apinte_weight.h"
#include "asterfort/assert.h"
#include "asterfort/cfadju.h"
#include "MeshTypes_type.h"
#include "asterfort/lcodrm.h"
#include "asterfort/projMaAndCheck.h"
#include "asterfort/interNodesInside.h"
#include "asterfort/interNodesEdge.h"
#include "asterfort/projInsideCell.h"
!
    real(kind=8), intent(in) :: proj_tole, dist_ratio
    integer(kind=8), intent(in) :: elem_dime
    integer(kind=8), intent(in) :: elem_mast_nbnode
    real(kind=8), intent(in) :: elem_mast_coor(3, MT_NNOMAX2D)
    character(len=8), intent(in) :: elem_mast_code
    integer(kind=8), intent(in) :: elem_slav_nbnode
    real(kind=8), intent(in) :: elem_slav_coor(3, MT_NNOMAX2D)
    character(len=8), intent(in) :: elem_slav_code
    real(kind=8), intent(out) :: poin_inte_ma(2, 8)
    real(kind=8), intent(out) :: poin_inte_es(2, 8)
    real(kind=8), intent(out) :: inte_weight
    integer(kind=8), intent(out) :: nb_poin_inte
    integer(kind=8), optional, intent(inout) :: inte_neigh_(4)
    integer(kind=8), optional, intent(inout) :: ierror_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Projection/intersection of elements in slave parametric space
!
! --------------------------------------------------------------------------------------------------
!
! In  proj_tole        : tolerance for projection
! In  elem_dime        : dimension of elements
! In  elem_mast_nbnode : number of nodes of master cell
! In  elem_mast_coor   : coordinates of master cell
! In  elem_mast_code   : code of master cell
! In  elem_slav_nbnode : number of nodes for slave cell
! In  elem_slav_coor   : coordinates of slave cell
! In  elem_slav_code   : code of slave cell
! Out poin_inte_es        : list (sorted) of intersection points
! Out inte_weight      : total weight of intersection
! Out nb_poin_inte     : number of intersection points
! IO  inte_neigh       : activation of neighbours of intersection
! Out iret : 0 - Ok and intersection is not void
!            1 - Ok and intersection is void
!            2 - Error during projection
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: debug = ASTER_FALSE
    aster_logical :: error
    real(kind=8) :: proj_coop(2, 9), coor_test(2)
    integer(kind=8) :: i_node, iret, nb_node_proj, test
    integer(kind=8) :: inte_neigh(4), nb_poin_inte_ma, nb_poin_inte_es
    real(kind=8) :: poin_inte_ma_tmp(2, 16)
    real(kind=8) :: poin_inte_es_tmp(2, 16)
    character(len=8) :: elin_mast_code
    integer(kind=8) :: elin_mast_nbnode
!
! --------------------------------------------------------------------------------------------------
!
    nb_poin_inte = 0
    nb_poin_inte_ma = 0
    nb_poin_inte_es = 0
    inte_weight = 0.d0
    iret = 0
    poin_inte_es_tmp = 0.d0
    poin_inte_ma_tmp = 0.d0
    poin_inte_es = 0.d0
    poin_inte_ma = 0.d0
    inte_neigh = 0
    error = ASTER_FALSE
!
    if (debug) then
        write (*, *) ". START Projection/intersection with raytracing"
    end if
!
    if (elem_mast_code .eq. "TR6" .or. elem_mast_code .eq. "TR7" &
        .or. elem_mast_code .eq. "TR1") then
        elin_mast_code = "TR3"
        elin_mast_nbnode = 3
    elseif (elem_mast_code .eq. "QU8" .or. elem_mast_code .eq. "QU9" &
            .or. elem_mast_code .eq. "Q12") then
        elin_mast_code = "QU4"
        elin_mast_nbnode = 4
    elseif (elem_mast_code .eq. "SE3" .or. elem_mast_code .eq. "SE4") then
        elin_mast_code = "SE2"
        elin_mast_nbnode = 2
    else
        elin_mast_code = elem_mast_code
        elin_mast_nbnode = elem_mast_nbnode
    end if

!
! - Projection on master cell
!
    call projMaAndCheck(proj_tole, dist_ratio, elem_dime, &
                        elin_mast_nbnode, elem_mast_coor, elin_mast_code, &
                        elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                        proj_coop, nb_node_proj, iret)
!
    if (iret == 2 .or. iret == 1) then
        error = ASTER_TRUE
        goto 100
    end if
!
! - Compute intersection in parametric master space
!
    call interNodesInside(proj_tole, elem_dime, &
                          elem_mast_code, elem_slav_code, &
                          proj_coop, nb_node_proj, nb_poin_inte_ma, poin_inte_ma_tmp, &
                          inte_neigh, poin_inte_es_tmp, iret)
!
    call interNodesEdge(proj_tole, elem_dime, &
                        elem_mast_code, elem_slav_code, &
                        proj_coop, nb_node_proj, nb_poin_inte_ma, poin_inte_ma_tmp, inte_neigh, &
                        poin_inte_es_tmp)
!
    if (nb_poin_inte_ma == 0 .or. iret == 1) then
        error = ASTER_TRUE
        iret = 1
        goto 100
    end if
!
    ASSERT(nb_poin_inte_ma .le. 16)
!
! - Sort list of intersection points
!
    if ((nb_poin_inte_ma .gt. 2 .and. elem_dime == 3) .or. &
        (nb_poin_inte_ma .ge. 2 .and. elem_dime == 2)) then
        nb_poin_inte_es = nb_poin_inte_ma
        call lcodrm(elem_dime, proj_tole, nb_poin_inte_ma, poin_inte_ma_tmp, poin_inte_es_tmp)
        nb_poin_inte_es = nb_poin_inte_ma
    end if
!
    if (nb_poin_inte_es == 0 .or. nb_poin_inte_es > 8) then
        iret = 1
        go to 99
    end if
!
! - All nodes have to be inside slave cell
!
    nb_poin_inte = nb_poin_inte_es
    ASSERT(nb_poin_inte <= 8)
    do i_node = 1, nb_poin_inte
        poin_inte_es(1, i_node) = poin_inte_es_tmp(1, i_node)
        poin_inte_ma(1, i_node) = poin_inte_ma_tmp(1, i_node)
        if (elem_dime .eq. 3) then
            poin_inte_es(2, i_node) = poin_inte_es_tmp(2, i_node)
            poin_inte_ma(2, i_node) = poin_inte_ma_tmp(2, i_node)
        end if
        !----- Test if point is inside element
        call cfadju(elem_slav_code, poin_inte_es(1, i_node), poin_inte_es(2, i_node), &
                    proj_tole, test)
        coor_test(:) = poin_inte_es_tmp(:, i_node)
        call projInsideCell(proj_tole*100, elem_dime, elem_slav_code, coor_test, test)
        if (test == 1) then
            print *, "projection to far: TO FIX"
            print *, "tole", proj_tole
            print *, "Node", i_node, "mail", elem_slav_code
            print *, "X", poin_inte_es_tmp(1, i_node)
            print *, "Y", poin_inte_es_tmp(2, i_node)
            ASSERT(.false.)
        end if
    end do
    !print*,"POINTINTEMAPAIR", poin_inte_ma(1,1:2)
    if (debug) then
        write (*, *) ". END Projection/intersection with raytracing"
    end if
!
! - Compute weight of intersection
!
    call apinte_weight(elem_dime, nb_poin_inte, poin_inte_es_tmp, inte_weight)
!
! - Error
!
100 continue
    if (error) then
        nb_poin_inte = 0
        poin_inte_es = 0.d0
        inte_neigh = 0
        inte_weight = 0.d0
    end if
!
! - No intersection exit
!
99  continue
    if (debug) then
        ASSERT(iret < 2)
    end if
!
! - Copy
!
    if (present(inte_neigh_)) then
        inte_neigh_(1:4) = inte_neigh(1:4)
    end if
    if (present(ierror_)) then
        ierror_ = iret
    end if
!
end subroutine
