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
subroutine apinte_chck2(proj_tole, elem_dime, &
                        elem_sside_nbnode, elem_sside_coor, &
                        elem_pside_nbnode, elem_pside_coor, elem_pside_code, &
                        norm_pside, norm_sside, &
                        proj_coor, l_inter)
!
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/reerel.h"
!
    real(kind=8), intent(in) :: proj_tole
    integer(kind=8), intent(in) :: elem_dime
    integer(kind=8), intent(in) :: elem_sside_nbnode
    real(kind=8), intent(in) :: elem_sside_coor(3, MT_NNOMAX2D)
    integer(kind=8), intent(in) :: elem_pside_nbnode
    real(kind=8), intent(in) :: elem_pside_coor(3, MT_NNOMAX2D)
    real(kind=8), intent(in) :: norm_pside(3)
    real(kind=8), intent(in) :: norm_sside(3)
    character(len=8), intent(in) :: elem_pside_code
    real(kind=8), intent(in) :: proj_coor(2, 4)
    aster_logical, intent(out) :: l_inter
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Check if intersection is void or not
!
! --------------------------------------------------------------------------------------------------
!
! In  proj_tole        : tolerance for projection
! In  elem_dime        : dimension of elements
! In  elem_sside_nbnode : number of nodes of start element
! In  elem_sside_coor   : coordinates of start element
! In  elem_slav_nbnode : number of nodes for target element
! In  elem_slav_coor   : coordinates of  target element
! In  elem_slav_code   : code of target element
! In  proj_coor        : projection of start nodes on target element (param)
! Out l_inter          : .true. if intersection is non-void
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: debug
    real(kind=8) :: no_sside_coor_1(3), no_sside_coor_2(3)
    real(kind=8) :: no_pside_coor_1(3), no_pside_coor_2(3)
    real(kind=8) :: ksi1(2), ksi2(2)
    real(kind=8) :: vect_pside(3), vect_sside(3), vect_stpj(3)
    real(kind=8) :: sig, sp_ns_vsp, sp_np_vsp
    integer(kind=8) :: i_node, i_dime, listNodenext(4)
!
! --------------------------------------------------------------------------------------------------
!
    debug = ASTER_FALSE
    l_inter = ASTER_TRUE
    no_sside_coor_1 = 0.d0
    no_sside_coor_2 = 0.d0
    no_pside_coor_1 = 0.d0
    no_pside_coor_2 = 0.d0
    vect_pside = 0.d0
    vect_sside = 0.d0
    vect_stpj = 0.d0

    do i_node = 1, elem_sside_nbnode-1
        listNodenext(i_node) = i_node+1
    end do
    listNodenext(elem_sside_nbnode) = 1
!
    do i_node = 1, elem_sside_nbnode
! ----- Get coordinates of start side nodes
        no_sside_coor_1(1:3) = 0.d0
        do i_dime = 1, elem_dime
            no_sside_coor_1(i_dime) = elem_sside_coor(i_dime, i_node)
        end do
        no_sside_coor_2(1:3) = 0.d0
        do i_dime = 1, elem_dime
            no_sside_coor_2(i_dime) = elem_sside_coor(i_dime, listNodenext(i_node))
        end do
! ----- Parametric coordinates of projection
        ksi1(1) = proj_coor(1, i_node)
        if (elem_dime .eq. 3) then
            ksi1(2) = proj_coor(2, i_node)
        end if
        ksi2(1) = proj_coor(1, listNodenext(i_node))
        if (elem_dime .eq. 3) then
            ksi2(2) = proj_coor(2, listNodenext(i_node))
        end if
! ----- Real coordinate of the the projection
        call reerel(elem_pside_code, elem_pside_nbnode, 3, elem_pside_coor, ksi1, &
                    no_pside_coor_1)
        call reerel(elem_pside_code, elem_pside_nbnode, 3, elem_pside_coor, ksi2, &
                    no_pside_coor_2)

! ----- Compute vector start side vect_sside
        vect_sside(:) = no_sside_coor_2(:)-no_sside_coor_1(:)

! ----- Compute vector projected side vect_pside
        vect_pside(:) = no_pside_coor_2(:)-no_pside_coor_1(:)

! ----- Compute vector start side to projected side
        vect_stpj(:) = no_pside_coor_1(:)-no_sside_coor_1(:)

! ----- Sign of colinear product vect_pside . vect_sside
        sig = 0.d0
        if (elem_dime .eq. 3) then
            sig = vect_pside(1)*vect_sside(1)+ &
                  vect_pside(2)*vect_sside(2)+ &
                  vect_pside(3)*vect_sside(3)
        elseif (elem_dime .eq. 2) then
            sig = vect_pside(1)*vect_sside(1)+ &
                  vect_pside(2)*vect_sside(2)
        else
            ASSERT(ASTER_FALSE)
        end if
! ----- Sign of colinear product  vect stpj . norm_start and vect stpj . norm_proj
        sp_ns_vsp = 0.d0
        sp_np_vsp = 0.d0
        if (elem_dime .eq. 3) then
            sp_ns_vsp = norm_sside(1)*vect_stpj(1)+ &
                        norm_sside(2)*vect_stpj(2)+ &
                        norm_sside(3)*vect_stpj(3)
            sp_np_vsp = norm_pside(1)*vect_stpj(1)+ &
                        norm_pside(2)*vect_stpj(2)+ &
                        norm_pside(3)*vect_stpj(3)
        elseif (elem_dime .eq. 2) then
            sp_ns_vsp = norm_sside(1)*vect_stpj(1)+ &
                        norm_sside(2)*vect_stpj(2)
            sp_np_vsp = norm_pside(1)*vect_stpj(1)+ &
                        norm_pside(2)*vect_stpj(2)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- check sign, negative => no intersection
        if (sig .lt. 0-proj_tole) then
            l_inter = ASTER_FALSE
            !print *, 'Problème check 1'
            exit
        elseif (sp_ns_vsp .lt. 0-proj_tole .and. sp_np_vsp .lt. 0-proj_tole) then
            !print*, "check"
            !print*, "no_pside_coor_1", no_pside_coor_1
            !print*, "no_sside_coor_1", no_sside_coor_1
            !print*, "vsp", vect_stpj(1:2)
            !print*, "np", norm_pside(1:2)
            !print*, "ns", norm_sside(1:2)
            !print*, "sp_ns_vsp", sp_ns_vsp
            !print*, "sp_np_vsp", sp_np_vsp
            !print*, "sig", sig
            !print *, 'Problème check 2'
            l_inter = ASTER_FALSE
            exit
        elseif (sp_ns_vsp .gt. 0+proj_tole .and. sp_np_vsp .gt. 0+proj_tole) then
            !print*, "no_pside_coor_1", no_pside_coor_1
            !print*, "no_sside_coor_1", no_sside_coor_1
            !print*, "check"
            !print*, "vsp", vect_stpj(1:2)
            !print*, "np", norm_pside(1:2)
            !print*, "ns", norm_sside(1:2)
            !print*, "sp_ns_vsp", sp_ns_vsp
            !print*, "sp_np_vsp", sp_np_vsp
            !print*, "sig", sig
            !print *, 'Problème check 2'
            l_inter = ASTER_FALSE
            exit
        end if

    end do
!
end subroutine
