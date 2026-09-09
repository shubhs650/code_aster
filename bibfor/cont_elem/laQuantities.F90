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
subroutine laQuantities(geom, param)
!
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
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
    integer(kind=8) :: i_node_slav, i_node_mast, i_dime, nb_lagr, nb_lagr_c, elem_dime, nb_node_slav
    integer(kind=8) :: jv_geom, jv_disp_incr, jv_disp, jv_geom_c, index
    integer(kind=8) :: jv_cont, jv_frot
    real(kind=8) :: depl_mast_incr(3, MT_NNOMAX2D), depl_slav_incr(3, MT_NNOMAX2D)
    real(kind=8) :: depl_mast_prev(3, MT_NNOMAX2D), depl_slav_prev(3, MT_NNOMAX2D)
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
    depl_slav_prev = 0.d0
    depl_mast_incr = 0.d0
    depl_slav_incr = 0.d0
!
    l_fric = lteatt('FROTTEMENT', 'OUI')
    elem_dime = geom%elem_dime
    nb_node_slav = geom%nb_node_slav
!
! - Slave nodes
!
    nb_lagr = 0
    nb_lagr_c = 0
    index = 0
    do i_node_slav = 1, nb_node_slav
        do i_dime = 1, elem_dime
            geom%coor_slav_init(i_dime, i_node_slav) = zr(jv_geom-1+index+i_dime)
            geom%coor_slav_pair(i_dime, i_node_slav) = zr(jv_geom_c-1+index+i_dime)
            depl_slav_prev(i_dime, i_node_slav) = zr(jv_disp-1+index+i_dime+nb_lagr)
            depl_slav_incr(i_dime, i_node_slav) = zr(jv_disp_incr-1+index+i_dime+nb_lagr)
        end do
!
        index = index+elem_dime
!
        if (geom%indi_lagc(i_node_slav) > 0) then
            nb_lagr_c = nb_lagr_c+1
            geom%lagc_slav_curr(nb_lagr_c) = &
                zr(jv_disp-1+index+nb_lagr+1)+zr(jv_disp_incr-1+index+nb_lagr+1)
            if (l_fric) then
                geom%lagf_slav_curr(1, nb_lagr_c) = &
                    zr(jv_disp-1+index+nb_lagr+2)+zr(jv_disp_incr-1+index+nb_lagr+2)
                if (geom%elem_dime == 3) then
                    geom%lagf_slav_curr(2, nb_lagr_c) = &
                        zr(jv_disp-1+index+nb_lagr+3)+zr(jv_disp_incr-1+index+nb_lagr+3)
                end if
            end if
            nb_lagr = nb_lagr+geom%indi_lagc(i_node_slav)
        end if
    end do
!
! - Master nodes
!
    do i_node_mast = 1, geom%nb_node_mast
        do i_dime = 1, elem_dime
            geom%coor_mast_init(i_dime, i_node_mast) = zr(jv_geom-1+index+i_dime)
            geom%coor_mast_pair(i_dime, i_node_mast) = zr(jv_geom_c-1+index+i_dime)
            depl_mast_prev(i_dime, i_node_mast) = zr(jv_disp-1+index+nb_lagr+i_dime)
            depl_mast_incr(i_dime, i_node_mast) = zr(jv_disp_incr-1+index+nb_lagr+i_dime)
        end do
        index = index+elem_dime
    end do
!
    geom%depl_slav_curr = depl_slav_prev+depl_slav_incr
    geom%depl_mast_curr = depl_mast_prev+depl_mast_incr
!
    geom%coor_slav_prev = geom%coor_slav_init+depl_slav_prev
    geom%coor_mast_prev = geom%coor_mast_init+depl_mast_prev
    geom%coor_slav_curr = geom%coor_slav_init+geom%depl_slav_curr
    geom%coor_mast_curr = geom%coor_mast_init+geom%depl_mast_curr
!
! - Times
!
    ! call jevech('PINSTMR', 'L', j_time)
    ! geom%time_prev = zr(j_time)
    ! call jevech('PINSTPR', 'L', j_time)
    ! geom%time_curr = zr(j_time)
!
! - COEF_CONT and COEF_FROT
!
    call jevech('PCCONTR', 'L', jv_cont)
    if (l_fric) then
        call jevech('PCFROTR', 'L', jv_frot)
    end if
!
    do i_node_slav = 1, geom%nb_node_slav
        param%coef_cont(i_node_slav) = zr(jv_cont-1+i_node_slav)
        if (l_fric) then
            param%coef_fric(i_node_slav) = zr(jv_frot-1+i_node_slav)
        end if
    end do
!
end subroutine
