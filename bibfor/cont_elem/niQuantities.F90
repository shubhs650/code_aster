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
subroutine niQuantities(geom, param)
!
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    type(ContactGeom), intent(inout) :: geom
    type(ContactParameters), intent(inout) :: param
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Get physical quantities
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_fric
    integer(kind=8) :: i_node_volu, i_node_mast, i_dime, elem_dime, nb_node_volu
    integer(kind=8) :: jv_geom, jv_disp_incr, jv_disp, jv_geom_c, index, j_time
    integer(kind=8) :: jv_cont, jv_frot, jcont, i_node_slav, map, iret, itab(8)
    real(kind=8) :: depl_mast_incr(3, MT_NNOMAX2D), depl_volu_incr(3, MT_NNOMAX3D)
    real(kind=8) :: depl_mast_prev(3, MT_NNOMAX2D), depl_volu_prev(3, MT_NNOMAX3D)
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PGEOMCR', 'L', jv_geom_c)
    call jevech('PDEPL_P', 'L', jv_disp_incr)
    call jevech('PDEPL_M', 'L', jv_disp)
!
! - Initializations
!
    depl_mast_prev = 0.d0
    depl_volu_prev = 0.d0
    depl_mast_incr = 0.d0
    depl_volu_incr = 0.d0
!
    l_fric = lteatt('FROTTEMENT', 'OUI')
    elem_dime = geom%elem_dime
    nb_node_volu = geom%nb_node_volu
!
! - Slave nodes
!
    index = 0
    do i_node_volu = 1, nb_node_volu
        do i_dime = 1, elem_dime
            geom%coor_volu_init(i_dime, i_node_volu) = zr(jv_geom-1+index+i_dime)
            geom%coor_volu_pair(i_dime, i_node_volu) = zr(jv_geom_c-1+index+i_dime)
            depl_volu_prev(i_dime, i_node_volu) = zr(jv_disp-1+index+i_dime)
            depl_volu_incr(i_dime, i_node_volu) = zr(jv_disp_incr-1+index+i_dime)
        end do
!
        index = index+elem_dime
!
    end do
!
! - Master nodes
!
    do i_node_mast = 1, geom%nb_node_mast
        do i_dime = 1, elem_dime
            geom%coor_mast_init(i_dime, i_node_mast) = zr(jv_geom-1+index+i_dime)
            geom%coor_mast_pair(i_dime, i_node_mast) = zr(jv_geom_c-1+index+i_dime)
            depl_mast_prev(i_dime, i_node_mast) = zr(jv_disp-1+index+i_dime)
            depl_mast_incr(i_dime, i_node_mast) = zr(jv_disp_incr-1+index+i_dime)
        end do
        index = index+elem_dime
    end do
!
    geom%depl_volu_curr = depl_volu_prev+depl_volu_incr
    geom%depl_mast_curr = depl_mast_prev+depl_mast_incr
!
    geom%coor_volu_prev = geom%coor_volu_init+depl_volu_prev
    geom%coor_mast_prev = geom%coor_mast_init+depl_mast_prev
    geom%coor_volu_curr = geom%coor_volu_init+geom%depl_volu_curr
    geom%coor_mast_curr = geom%coor_mast_init+geom%depl_mast_curr
!
! - Map Volu to Surf
!
    call jevech('PCONFR', 'L', jcont)
    ASSERT(nint(zr(jcont+50)) == geom%nb_node_slav)

    do i_node_slav = 1, geom%nb_node_slav
        map = nint(zr(jcont+50+i_node_slav))
        geom%mapVolu2Surf(i_node_slav) = map
        geom%coor_slav_init(1:elem_dime, i_node_slav) = geom%coor_volu_init(1:elem_dime, map)
        geom%coor_slav_prev(1:elem_dime, i_node_slav) = geom%coor_volu_prev(1:elem_dime, map)
        geom%coor_slav_curr(1:elem_dime, i_node_slav) = geom%coor_volu_curr(1:elem_dime, map)
        geom%coor_slav_pair(1:elem_dime, i_node_slav) = geom%coor_volu_pair(1:elem_dime, map)
        geom%depl_slav_curr(1:elem_dime, i_node_slav) = geom%depl_volu_curr(1:elem_dime, map)
    end do
!
! - Times
!
    call jevech('PINSTMR', 'L', j_time)
    geom%time_prev = zr(j_time)
    call jevech('PINSTPR', 'L', j_time)
    geom%time_curr = zr(j_time)
!
! - COEF_CONT and COEF_FROT
!
    call tecach('OON', 'PCCONTR', 'L', iret, nval=8, itab=itab)
    ASSERT(iret == 3)
    jv_cont = itab(1)
    if (l_fric) then
        call tecach('OON', 'PCFROTR', 'L', iret, nval=8, itab=itab)
        ASSERT(iret == 3)
        jv_frot = itab(1)
    end if
    do i_node_slav = 1, geom%nb_node_slav
        map = geom%mapVolu2Surf(i_node_slav)
        ASSERT(zl(itab(8)-1+map))
        param%coef_cont(i_node_slav) = zr(jv_cont-1+map)
        if (l_fric) then
            param%coef_fric(i_node_slav) = zr(jv_frot-1+map)
        end if
    end do

!
end subroutine
