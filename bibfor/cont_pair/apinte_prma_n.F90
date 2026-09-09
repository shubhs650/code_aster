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
subroutine apinte_prma_n(proj_tole, elem_dime, &
                         elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                         elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                         proj_coor, nb_node_proj, iret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/gapGetParamCoor.h"
#include "MeshTypes_type.h"
#include "asterfort/mmnewd.h"
!
    real(kind=8), intent(in) :: proj_tole
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
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: debug = ASTER_FALSE
    integer(kind=8) :: i_node, i_dime, elem_mast_line_nbnode
    real(kind=8) :: nosl_coor(3), ksi1_line, ksi2_line
    real(kind=8) :: ksi1, ksi2, tau1(3), tau2(3), para_coor(2, MT_NNOMAX2D)
    real(kind=8) :: norm_slav(3)
    character(len=8) :: elem_mast_line_code

!
! --------------------------------------------------------------------------------------------------
!
    if (debug) then
        write (*, *) ".. Project slaves nodes in master element parametric space"
    end if
    ASSERT(elem_slav_nbnode <= size(proj_coor, 2))
    proj_coor = 0.d0
!
! --- Projection only vertex of slave cell
!
    nb_node_proj = elem_slav_nbnode
    if (elem_slav_code == "SE3" .or. elem_slav_code == "SE4") then
        nb_node_proj = 2
    elseif (elem_slav_code == "QU8" .or. elem_slav_code == "QU9" &
            .or. elem_slav_code == "Q12") then
        nb_node_proj = 4
    elseif (elem_slav_code == "TR6" .or. elem_slav_code == "TR7" &
            .or. elem_slav_code == "TR1") then
        nb_node_proj = 3
    end if

    if (elem_mast_code(1:2) == "SE") then
        elem_mast_line_code = "SE2"
        elem_mast_line_nbnode = 2
    elseif (elem_mast_code(1:2) == "TR") then
        elem_mast_line_code = "TR3"
        elem_mast_line_nbnode = 3
    elseif (elem_mast_code(1:2) == "QU" .or. elem_mast_code == "Q12") then
        elem_mast_line_code = "QU4"
        elem_mast_line_nbnode = 4
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Get parametric slave coordinates
!
    call gapGetParamCoor(elem_slav_code, para_coor)
!
    do i_node = 1, nb_node_proj
! ----- Get coordinates of slaves nodes (real space)
        nosl_coor = 0.d0
        do i_dime = 1, elem_dime
            nosl_coor(i_dime) = elem_slav_coor(i_dime, i_node)
        end do
!
! --------- Compute outward normal
!
        call apnorm(elem_slav_nbnode, elem_slav_code, elem_dime, elem_slav_coor, &
                    para_coor(1, i_node), para_coor(2, i_node), norm_slav)
!
! ----- Projection on master element
        call mmnewd(elem_mast_code, elem_mast_nbnode, elem_dime, &
                    elem_mast_coor, nosl_coor, 75, &
                    proj_tole, norm_slav, &
                    ksi1, ksi2, &
                    tau1, tau2, iret)
        if (iret == 1) then
!
! ----- Try with linearization
            if (debug) then
                write (*, *) "... Error in projection: try with linear cell"
            end if

            call mmnewd(elem_mast_line_code, elem_mast_line_nbnode, elem_dime, &
                        elem_mast_coor, nosl_coor, 75, &
                        proj_tole, norm_slav, &
                        ksi1_line, ksi2_line, &
                        tau1, tau2, iret)
            call mmnewd(elem_mast_code, elem_mast_nbnode, elem_dime, &
                        elem_mast_coor, nosl_coor, 75, &
                        proj_tole, norm_slav, &
                        ksi1, ksi2, &
                        tau1, tau2, iret, &
                        ksi1_init=ksi1_line, ksi2_init=ksi2_line)
!
            if (iret == 1) then
                proj_coor = 0.d0
                nb_node_proj = 0
                iret = 2
                exit
            end if
        end if
!
! ----- Get parametric coordinates of projection
        proj_coor(1, i_node) = ksi1
        if (elem_dime .eq. 3) then
            proj_coor(2, i_node) = ksi2
        end if

!
        if (debug) then
            print *, "ES REF: ", para_coor(1:2, i_node)
            print *, "NORM: ", norm_slav
            print *, "PROJ: ", proj_coor(1:2, i_node)
        end if
    end do
!
end subroutine
