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
subroutine projMaAndCheck(proj_tole, dist_ratio, elem_dime, &
                          elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                          elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                          proj_coor, nb_node_proj, iret)
!
    use contact_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apinte_chck2.h"
#include "asterfort/apinte_norm.h"
#include "MeshTypes_type.h"
#include "asterfort/apinte_prma_n.h"
!
    real(kind=8), intent(in) :: proj_tole, dist_ratio
    integer(kind=8), intent(in) :: elem_dime
    integer(kind=8), intent(in) :: elem_mast_nbnode
    real(kind=8), intent(in) :: elem_mast_coor(3, MT_NNOMAX2D)
    integer(kind=8), intent(in) :: elem_slav_nbnode
    real(kind=8), intent(in) :: elem_slav_coor(3, MT_NNOMAX2D)
    character(len=8), intent(in) :: elem_mast_code, elem_slav_code
    real(kind=8), intent(out) :: proj_coor(2, 9)
    integer(kind=8), intent(out) :: iret, nb_node_proj
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Project slave nodes in master element parametric space
!
! --------------------------------------------------------------------------------------------------
!
! In  proj_tole        : tolerance for projection
! In  elem_dime        : dimension of elements
! In  elem_mast_nbnode : number of nodes of master element
! In  elem_mast_coor   : coordinates of master element
! In  elem_slav_nbnode : number of nodes for slave element
! In  elem_slav_coor   : coordinates of slave element
! In  elem_mast_code   : code of master element
! Out proj_coor        : projection of slave nodes on master element
! Out iret : 0 - Ok and intersection not void
!            1 - Ok but intersection is void
!            2 - Error during projection
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: mast_norm(3), slav_norm(3), ps
    real(kind=8) :: bar_ma(3), bar_sl(3), hf_ma, hf_sl
    aster_logical :: l_inter
!
    iret = 0
    proj_coor = 0.d0
    nb_node_proj = 0
!
! - Compute norms at barycenter
!
    call apinte_norm(elem_dime, &
                     elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                     elem_slav_coor, elem_slav_code, &
                     mast_norm, slav_norm)
!
! - Linear normal are orthonormal - exit
!
    ps = mast_norm(1)*slav_norm(1)+mast_norm(2)*slav_norm(2)+mast_norm(3)*slav_norm(3)
    if (abs(ps) <= proj_tole) then
        iret = 1
        go to 99
    end if
!
! - If distance between barycenter is to high -> exit
!
    if (dist_ratio > 0) then
        hf_ma = diameter(elem_mast_nbnode, elem_mast_coor)
        bar_ma = barycenter(elem_mast_nbnode, elem_mast_coor)
        hf_sl = diameter(elem_slav_nbnode, elem_slav_coor)
        bar_sl = barycenter(elem_slav_nbnode, elem_slav_coor)
        !
        if (norm2(bar_ma-bar_sl) >= 2*dist_ratio*max(hf_ma, hf_sl)) then
            iret = 1
            go to 99
        end if
    end if
!
! - Project slave nodes in master cell parametric space using raytracing
!
    call apinte_prma_n(proj_tole, elem_dime, &
                       elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                       elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                       proj_coor, nb_node_proj, iret)
!
! - Check if intersection is void or not
!
    call apinte_chck2(proj_tole, elem_dime, &
                      nb_node_proj, elem_slav_coor, &
                      elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                      mast_norm, slav_norm, &
                      proj_coor, l_inter)
    if (.not. l_inter) then
        iret = 1
        goto 99
    end if
!
99  continue
!
end subroutine
