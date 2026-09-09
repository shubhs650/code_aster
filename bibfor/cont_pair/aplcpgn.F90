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
subroutine aplcpgn(mesh, newgeo, &
                   mastConxInvName, mastNeighName, slavNeighName, &
                   pair_tole, dist_ratio, &
                   nb_elem_mast, list_elem_mast, nb_elem_slav, list_elem_slav, list_node_mast, &
                   nb_node_mast, meshPairing)
!
    use mesh_type
    use mesh_cell_module
    use MeshPairing_module
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/ap_infast_n.h"
#include "asterfort/apcoor.h"
#include "asterfort/aptype.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/mesh_pairing_type.h"
#include "asterfort/prjint_ray.h"
#include "asterfort/testvois.h"
#include "asterfort/utmess.h"
#include "MeshTypes_type.h"
#include "Contact_type.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: newgeo
    character(len=24), intent(in) :: mastConxInvName
    character(len=24), intent(in) :: mastNeighName, slavNeighName
    real(kind=8), intent(in) :: pair_tole, dist_ratio
    integer(kind=8), intent(in) :: nb_elem_slav
    integer(kind=8), intent(in) :: nb_elem_mast
    integer(kind=8), intent(in) :: nb_node_mast
    integer(kind=8), intent(in) :: list_elem_mast(nb_elem_mast)
    integer(kind=8), intent(in) :: list_elem_slav(nb_elem_slav)
    integer(kind=8), intent(in) :: list_node_mast(nb_node_mast)
    type(MESH_PAIRING), intent(inout) :: meshPairing
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Pairing by PANG method
!
! NO OUTPUT (only for verbose debug)
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  newgeo           : name of field for geometry update from initial coordinates of nodes
! In  zone             : name of contact zone sd
! In  pair_tole        : tolerance for pairing
! In  nb_elem_mast     : number of master elements on current zone
! In  nb_elem_slav     : number of slave elements on current zone
! In  pair_method      : pairing method => "RAPIDE" or "ROBUSTE"
! In  list_elem_mast   : name of datastructure for list of master elements on current zone
! In  list_elem_slav   : name of datastructure for list of slave elements on current zone
! IO  nb_pair_zone     : number of contact elements
! IO  list_pair_zone   : list of contact elements
! IO  list_nbptit_zone :
! IO  list_ptitsl_zone : list of intersection points in slave elements
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, vali(2)
    integer(kind=8) :: elem_slav_nbnode, elem_slav_nume, elem_slav_dime, elem_slav_indx
    integer(kind=8) :: elem_mast_nbnode, elem_mast_nume, elem_mast_dime, elem_mast_indx
    character(len=8) :: elem_mast_code, elem_slav_code
    character(len=8) :: elem_slav_type, elem_mast_type
    real(kind=8) :: elem_mast_coor(3, MT_NNOMAX2D), elem_slav_coor(3, MT_NNOMAX2D)
    integer(kind=8) :: nb_pair, nb_poin_inte, i_dim, i_pt
    integer(kind=8) :: i_mast_neigh, i_slav_start, i_mast_start, i_find_mast
    integer(kind=8) :: i_slav_neigh
    real(kind=8) :: inte_weight
    real(kind=8) :: poin_inte_sl(SIZE_MAX_INTE_SL), poin_inte_sl2(2, MAX_NB_INTE)
    real(kind=8) :: poin_inte_ma(SIZE_MAX_INTE_SL), point_inte_tmp(SIZE_MAX_INTE_SL)
    character(len=8) :: elem_slav_name, elem_name
    integer(kind=8) :: nb_slav_start, nb_find_mast, nb_mast_start
    integer(kind=8) :: elem_start, elem_nume
    integer(kind=8) :: slav_indx_mini, mast_indx_mini, slav_indx_maxi, mast_indx_maxi
    integer(kind=8) :: elem_neigh_indx, mast_find_indx, elem_slav_neigh, elem_mast_neigh
    aster_logical :: l_recup, debug, pair_exist
    integer(kind=8), pointer :: mast_find_flag(:) => null()
    integer(kind=8), pointer :: elem_slav_flag(:) => null()
    integer(kind=8), pointer :: v_sdappa_slne(:) => null()
    integer(kind=8), pointer :: v_sdappa_mane(:) => null()
    integer(kind=8) :: list_slav_master(4)
    integer(kind=8) :: nb_mast_neigh, nb_slav_neigh
    integer(kind=8) :: inte_neigh(4)
    integer(kind=8) :: jv_geom, elem_type_nume
    real(kind=8) :: list_slav_weight(4), weight_test, tole_weight
    integer(kind=8), pointer :: v_mesh_typmail(:) => null()
    integer(kind=8), pointer :: v_mesh_connex(:) => null()
    integer(kind=8), pointer :: v_connex_lcum(:) => null()
    real(kind=8), pointer :: li_pt_inte_sl(:) => null()
    real(kind=8), pointer :: li_pt_inte_ma(:) => null()
    integer(kind=8), pointer :: list_pair(:) => null()
    integer(kind=8), pointer :: li_nb_pt_inte_sl(:) => null()
    integer(kind=8), pointer :: list_find_mast(:) => null()
    integer(kind=8), pointer :: elem_slav_start(:) => null()
    integer(kind=8), pointer :: elem_mast_start(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - some initializations
!
    debug = meshPairing%debug
    pair_exist = ASTER_TRUE
    inte_neigh(1:4) = 0
    list_slav_master(1:4) = 0
    list_slav_weight(1:4) = 0.d0
    nb_pair = 0

    mast_indx_maxi = maxval(list_elem_mast)
    slav_indx_maxi = maxval(list_elem_slav)
    mast_indx_mini = minval(list_elem_mast)
    slav_indx_mini = minval(list_elem_slav)
!
! - Access to updated geometry
!
    call jeveuo(newgeo(1:19)//'.VALE', 'L', jv_geom)
!
! - Access to mesh
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=v_mesh_typmail)
    call jeveuo(mesh//'.CONNEX', 'L', vi=v_mesh_connex)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=v_connex_lcum)
!
! - Objects for flags
!
    AS_ALLOCATE(vi=elem_slav_flag, size=slav_indx_maxi+1-slav_indx_mini)
    AS_ALLOCATE(vi=mast_find_flag, size=mast_indx_maxi+1-mast_indx_mini)
    AS_ALLOCATE(vr=li_pt_inte_sl, size=nb_elem_mast*nb_elem_slav*SIZE_MAX_INTE_SL)
    AS_ALLOCATE(vr=li_pt_inte_ma, size=nb_elem_mast*nb_elem_slav*SIZE_MAX_INTE_SL)
    AS_ALLOCATE(vi=list_pair, size=2*nb_elem_slav*nb_elem_mast)
    AS_ALLOCATE(vi=li_nb_pt_inte_sl, size=nb_elem_slav*nb_elem_mast)
    AS_ALLOCATE(vi=list_find_mast, size=nb_elem_mast)
    AS_ALLOCATE(vi=elem_slav_start, size=nb_elem_slav)
    AS_ALLOCATE(vi=elem_mast_start, size=nb_elem_slav)
    list_find_mast = 0

! - Object for neighbours (inverse connectivity)
    call jeveuo(mastNeighName, 'L', vi=v_sdappa_mane)
    call jeveuo(slavNeighName, 'L', vi=v_sdappa_slne)

! - while loop on the existence of a pair slave-master
    do while (pair_exist)
        if (debug) then
            WRITE (6, *) "Get cells for starting search"
            WRITE (6, *) "============================="
        end if

! ----- Search by computing the minimum distance between the barycenters
        call ap_infast_n(mesh, newgeo, pair_tole, dist_ratio, nb_elem_mast, &
                         list_elem_mast, nb_elem_slav, list_elem_slav, elem_slav_flag, &
                         nb_mast_start, elem_mast_start, nb_slav_start, elem_slav_start, &
                         mastConxInvName, list_node_mast, nb_node_mast)
        if (debug) then
            WRITE (6, *) "  Starting points (Master): ", nb_mast_start, elem_mast_start
            WRITE (6, *) "  Starting points (Slave) : ", nb_slav_start, elem_slav_start
        end if

        if (nb_slav_start == 0) then
            pair_exist = ASTER_FALSE
        end if

        do while (nb_slav_start > 0)
!
! ----- Get slave element start
!
            elem_slav_nume = elem_slav_start(1)
            elem_slav_indx = elem_slav_nume+1-slav_indx_mini
            elem_type_nume = v_mesh_typmail(elem_slav_nume)
            call jenuno(jexnum('&CATA.TM.NOMTM', elem_type_nume), elem_slav_type)
!
! ----- Shift list of slave element start
!
            do i_slav_start = 1, nb_slav_start-1
                elem_slav_start(i_slav_start) = elem_slav_start(i_slav_start+1)
            end do
            nb_slav_start = nb_slav_start-1
!
! ----- Get informations about slave element
!
            call aptype(elem_slav_type, &
                        elem_slav_nbnode, elem_slav_code, elem_slav_dime)
!
! ----- Get coordinates of slave element
!
            call apcoor(v_mesh_connex, v_connex_lcum, jv_geom, &
                        elem_slav_nume, elem_slav_nbnode, elem_slav_dime, &
                        elem_slav_coor)

            if (debug) then
                elem_slav_name = int_to_char8(elem_slav_nume)
                write (*, *) "Current slave element      : ", elem_slav_nume
                write (6, *) " Coordinates (global frame): ", &
                    elem_slav_coor(1:3, 1:elem_slav_nbnode)
            end if
!
! ----- Number of neighbours
!
            if (elem_slav_dime == 2) then
                nb_slav_neigh = 2
            elseif (elem_slav_code == 'TR3' .or. &
                    elem_slav_code == 'TR6' .or. &
                    elem_slav_code == 'TR7') then
                nb_slav_neigh = 3
            elseif (elem_slav_code == 'QU4' .or. &
                    elem_slav_code == 'QU8' .or. &
                    elem_slav_code == 'QU9') then
                nb_slav_neigh = 4
            else
                ASSERT(ASTER_FALSE)
            end if
!
            if (debug) then
                write (6, *) "Potential number of neighbours: ", nb_slav_neigh
                do i_slav_neigh = 1, nb_slav_neigh
                    elem_nume = v_sdappa_slne((elem_slav_indx-1)*4+i_slav_neigh)
                    if (elem_nume .ne. 0) then
                        elem_name = int_to_char8(elem_nume)
                    else
                        elem_name = 'None'
                    end if
                    write (6, *) "Neighbour (", i_slav_neigh, "): ", elem_nume
                end do
            end if
!
            list_slav_master(1:nb_slav_neigh) = 0
            list_slav_weight(1:4) = 0.d0
!
! ----- Get master element to start
!
            elem_start = elem_mast_start(1)
            mast_find_indx = elem_start+1-mast_indx_mini
            if (debug) then
                WRITE (6, *) "Avant décalage"
                WRITE (6, *) " => ", nb_mast_start, elem_mast_start(1:nb_mast_start)
            end if

! ----- Shift list of master element start
!
            do i_mast_start = 1, nb_mast_start-1
                elem_mast_start(i_mast_start) = elem_mast_start(i_mast_start+1)
            end do
            nb_mast_start = nb_mast_start-1
            if (debug) then
                WRITE (6, *) "Après décalage"
                WRITE (6, *) " => ", nb_mast_start, elem_mast_start(1:nb_mast_start)
            end if
!
! ----- Management of list of master elements: first element to seek
!
            list_find_mast(1) = elem_start
            nb_find_mast = 1
            mast_find_flag(mast_find_indx) = 1
            if (debug) then
                write (6, *) "Master cell to start: ", elem_start
            end if
!
! ----- Initialization list of contact pairs
!
            l_recup = ASTER_TRUE

! -----     Loop on master elements => Look for the master elements
!
            do while (nb_find_mast > 0)
!
                inte_weight = 0.d0
                nb_poin_inte = 0
                poin_inte_sl = 0.d0
!
! ------------- Get master element
!
                elem_mast_nume = list_find_mast(1)
                elem_mast_indx = elem_mast_nume+1-mast_indx_mini
                elem_type_nume = v_mesh_typmail(elem_mast_nume)
                call jenuno(jexnum('&CATA.TM.NOMTM', elem_type_nume), elem_mast_type)
!
                if (debug) then
                    elem_name = int_to_char8(elem_mast_nume)
                    write (*, *) "Current master element: ", elem_mast_nume
                    write (6, *) " Coordinates (global frame): "
                end if
!
! ------------- Shift list of master elements (on supprime de la liste)
!
                do i_find_mast = 1, nb_find_mast-1
                    list_find_mast(i_find_mast) = list_find_mast(i_find_mast+1)
                end do
                nb_find_mast = nb_find_mast-1
!
! --------- Get informations about master element
!
                call aptype(elem_mast_type, &
                            elem_mast_nbnode, elem_mast_code, elem_mast_dime)
!
! --------- Get coordinates of master element
!
                call apcoor(v_mesh_connex, v_connex_lcum, jv_geom, &
                            elem_mast_nume, elem_mast_nbnode, elem_mast_dime, &
                            elem_mast_coor)
!
! --------- Projection/intersection of elements in slave parametric space
!
                call prjint_ray(pair_tole, dist_ratio, elem_slav_dime, &
                                elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                                elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                                poin_inte_ma, poin_inte_sl, inte_weight, nb_poin_inte, &
                                inte_neigh_=inte_neigh, ierror_=iret)
!
                if (iret == 2) then
                    vali(1) = elem_slav_nume
                    vali(2) = elem_mast_nume
                    call utmess('A', 'CONTACT4_6', ni=2, vali=vali)
!
                    inte_weight = 0.d0
                    nb_poin_inte = 0.d0
                    go to 101
                end if

                ASSERT(nb_poin_inte .le. 8)
!
101             continue
!
! --------- Add element paired
!
                if (inte_weight > pair_tole .and. iret == 0) then
                    if (elem_slav_dime == 2) then
                        ! In 2D, we have to reshape since function have fixed size now
                        point_inte_tmp = poin_inte_sl
                        poin_inte_sl = 0.d0
                        poin_inte_sl(1:SIZE_MAX_INTE_SL/2) = point_inte_tmp(1:SIZE_MAX_INTE_SL:2)
                        point_inte_tmp = poin_inte_ma
                        poin_inte_ma = 0.d0
                        poin_inte_ma(1:SIZE_MAX_INTE_SL/2) = point_inte_tmp(1:SIZE_MAX_INTE_SL:2)
                    end if

                    nb_pair = nb_pair+1
                    if (debug) then
                        WRITE (6, *) "Add pair: ", nb_pair, &
                            "(", elem_slav_nume, "-", elem_mast_nume, ")"
                        WRITE (6, *) "  Nb points integrations    : ", nb_poin_inte
                        WRITE (6, *) "  Coef. points integrations : ", poin_inte_sl
                    end if

                    ASSERT(nb_pair .le. nb_elem_slav*nb_elem_mast)
                    list_pair(2*(nb_pair-1)+1) = elem_slav_nume
                    list_pair(2*(nb_pair-1)+2) = elem_mast_nume
                    li_nb_pt_inte_sl(nb_pair) = nb_poin_inte
                    ASSERT(nb_poin_inte .le. 8)
                    li_pt_inte_sl((nb_pair-1)*SIZE_MAX_INTE_SL+1: &
                                  (nb_pair-1)*SIZE_MAX_INTE_SL+SIZE_MAX_INTE_SL) = poin_inte_sl
                    li_pt_inte_ma((nb_pair-1)*SIZE_MAX_INTE_SL+1: &
                                  (nb_pair-1)*SIZE_MAX_INTE_SL+SIZE_MAX_INTE_SL) = poin_inte_ma
                    !print*,"LIPTMA_APLC", li_pt_inte_ma((nb_pair-1)*SIZE_MAX_INTE_SL+1:&
                    !            (nb_pair-1)*SIZE_MAX_INTE_SL+2)

                    poin_inte_sl2 = 0.d0
                    if (elem_slav_dime == 2) then
                        do i_pt = 1, nb_poin_inte
                            poin_inte_sl2(1, i_pt) = poin_inte_sl(i_pt)
                        end do
                    else
                        do i_pt = 1, nb_poin_inte
                            do i_dim = 1, 2
                                poin_inte_sl2(i_dim, i_pt) = poin_inte_sl(2*i_pt+i_dim)
                            end do
                        end do
                    end if
                    call pairAdd(elem_slav_nume, elem_mast_nume, &
                                 nb_poin_inte, poin_inte_sl2, &
                                 meshPairing)
                end if

                if (debug) then
                    write (*, *) ". Contact pair: ", elem_slav_nume, elem_mast_nume, &
                        " (weight: ", inte_weight, ", nb point inter: ", nb_poin_inte, ")"
                end if
!
! --------- Find neighbour of current master element
!
                if (inte_weight > pair_tole .or. l_recup) then
!
! ------------- Number of neighbours
!
                    if (elem_mast_code == 'SE2' .or. elem_mast_code == 'SE3' &
                        .or. elem_mast_code == 'SE4') then
                        nb_mast_neigh = 2
                        tole_weight = 0.5
                    elseif (elem_mast_code == 'TR3' .or. elem_mast_code == 'TR6' &
                            .or. elem_mast_code == 'TR7' .or. elem_mast_code == 'TR1') then
                        nb_mast_neigh = 3
                        tole_weight = 0.05
                    elseif (elem_mast_code == 'QU4' .or. elem_mast_code == 'QU8' .or. &
                            elem_mast_code == 'QU9' .or. elem_mast_code == 'Q12') then
                        nb_mast_neigh = 4
                        tole_weight = 0.4
                    else
                        ASSERT(ASTER_FALSE)
                    end if
!
! ------------- Prepare next master element
!
                    ! WRITE (6, *) "Prepare next master element"
                    do i_mast_neigh = 1, nb_mast_neigh
                        elem_mast_neigh = v_sdappa_mane((elem_mast_indx-1)*4+i_mast_neigh)
                        elem_neigh_indx = elem_mast_neigh+1-mast_indx_mini
                        if (elem_mast_neigh .ne. 0) then
                            if (mast_find_flag(elem_neigh_indx) == 0) then
                                list_find_mast(nb_find_mast+1) = elem_mast_neigh
                                nb_find_mast = nb_find_mast+1
                                if (debug) then
                                    WRITE (6, *) " => added: ", nb_find_mast, &
                                        list_find_mast(nb_find_mast)
                                end if
                                mast_find_flag(elem_neigh_indx) = 1
                            end if
                        end if
                    end do

!
! ------------- Prepare next slave element: higher weight
!
                    if (debug) then
                        WRITE (6, *) "Prepare next slave element"
                    end if
                    do i_slav_neigh = 1, nb_slav_neigh
                        elem_slav_neigh = v_sdappa_slne((elem_slav_indx-1)*4+i_slav_neigh)
                        elem_neigh_indx = elem_slav_neigh+1-slav_indx_mini
                        if (debug) then
                            WRITE (6, *) " < Index:", i_slav_neigh
                            WRITE (6, *) " < cellNeighNume:", elem_slav_neigh
                            WRITE (6, *) " < cellNeighIndx: ", elem_neigh_indx
                            WRITE (6, *) " < cellSlavFlag: ", elem_slav_flag(elem_neigh_indx)
                            WRITE (6, *) " < inteNeigh: ", inte_neigh(i_slav_neigh)
                        end if
                        if (elem_slav_neigh .ne. 0) then
                            if (inte_neigh(i_slav_neigh) == 1 &
                                .and. elem_slav_flag(elem_neigh_indx) .ne. 1 &
                                .and. list_slav_weight(i_slav_neigh) .lt. tole_weight) then
                                weight_test = 0.d0
                                ! IS IT NECESSARY WITH RAY_TRACING ?
                                !call testvois(jv_geom       , elem_slav_type,&
                                !                elem_mast_coor, elem_mast_code, elem_slav_nume,&
                                !                pair_tole     , weight_test,    v_mesh_connex ,&
                                !                v_connex_lcum)
                                !if (weight_test > list_slav_weight(i_slav_neigh).and.&
                                !    weight_test > pair_tole) then
                                list_slav_master(i_slav_neigh) = elem_mast_nume
                                if (debug) then
                                    WRITE (6, *) " => added: ", elem_mast_nume
                                end if
                                !   list_slav_weight(i_slav_neigh) = weight_test
                            end if
                        end if
                    end do
                    l_recup = ASTER_FALSE
                end if

            end do
!
! ----- Next elements
!
            if (debug) then
                write (*, *) 'Prepare next elements - Nb: ', nb_slav_neigh
            end if
            do i_slav_neigh = 1, nb_slav_neigh
                elem_slav_neigh = v_sdappa_slne((elem_slav_indx-1)*4+i_slav_neigh)
                elem_neigh_indx = elem_slav_neigh+1-slav_indx_mini
                if (debug) then
                    write (*, *) 'Next elements - Current: ', i_slav_neigh, elem_slav_neigh, &
                        list_slav_master(i_slav_neigh), elem_slav_flag(elem_neigh_indx)
                end if
                if (elem_slav_neigh .ne. 0) then
                    if (list_slav_master(i_slav_neigh) .ne. 0 .and. &
                        elem_slav_flag(elem_neigh_indx) .ne. 1) then
                        elem_slav_start(nb_slav_start+1) = elem_slav_neigh
                        nb_slav_start = nb_slav_start+1
                        elem_slav_flag(elem_neigh_indx) = 1
                        elem_mast_start(nb_mast_start+1) = list_slav_master(i_slav_neigh)
                        nb_mast_start = nb_mast_start+1
                    end if
                end if
            end do
            mast_find_flag(1:mast_indx_maxi+1-mast_indx_mini) = 0

        end do

        !pair_exist = ASTER_FALSE
    end do

! - DEALLOCATE
    AS_DEALLOCATE(vi=mast_find_flag)
    AS_DEALLOCATE(vi=elem_slav_flag)
    AS_DEALLOCATE(vr=li_pt_inte_sl)
    AS_DEALLOCATE(vr=li_pt_inte_ma)
    AS_DEALLOCATE(vi=list_pair)
    AS_DEALLOCATE(vi=li_nb_pt_inte_sl)
    AS_DEALLOCATE(vi=list_find_mast)
    AS_DEALLOCATE(vi=elem_slav_start)
    AS_DEALLOCATE(vi=elem_mast_start)

    call jedema()
end subroutine
