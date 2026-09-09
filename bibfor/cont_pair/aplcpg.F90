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

! aslint: disable=W1306
!
subroutine aplcpg(mesh, newgeo, sdappa, i_zone, pair_tole, &
                  nb_elem_mast, list_elem_mast, nb_elem_slav, list_elem_slav, &
                  nb_pair_zone, list_pair_zone, list_nbptit_zone, list_ptitsl_zone, &
                  list_ptitma_zone, li_ptgausma_zone, i_proc, nb_proc, pair_method)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "asterfort/aprtpm.h"
#include "asterfort/jexatr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/apcoor.h"
#include "asterfort/aptype.h"
#include "asterfort/prjint.h"
#include "asterfort/clpoma.h"
#include "asterfort/assert.h"
#include "asterfort/ap_infast.h"
#include "asterfort/apprin.h"
#include "asterfort/wkvect.h"
#include "asterfort/codent.h"
#include "asterfort/testvois.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/apsave_pair.h"
#include "asterfort/apsave_patch.h"
#include "MeshTypes_type.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: newgeo
    character(len=19), intent(in) :: sdappa
    integer(kind=8), intent(in) :: i_zone
    real(kind=8), intent(in) :: pair_tole
    integer(kind=8), intent(in) :: nb_elem_slav
    integer(kind=8), intent(in) :: nb_elem_mast
    integer(kind=8), intent(in) :: list_elem_mast(nb_elem_mast)
    integer(kind=8), intent(in) :: list_elem_slav(nb_elem_slav)
    integer(kind=8), intent(inout) :: nb_pair_zone
    integer(kind=8), pointer :: list_pair_zone(:)
    integer(kind=8), pointer :: list_nbptit_zone(:)
    real(kind=8), pointer :: list_ptitsl_zone(:)
    real(kind=8), pointer :: list_ptitma_zone(:)
    real(kind=8), pointer :: li_ptgausma_zone(:)
    integer(kind=8), intent(in) :: i_proc
    integer(kind=8), intent(in) :: nb_proc
    character(len=24), intent(in) :: pair_method
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Pairing by PANG method
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  newgeo           : name of field for geometry update from initial coordinates of nodes
! In  sdappa           : name of pairing datastructure
! In  i_zone           : index of contact zone
! In  pair_tole        : tolerance for pairing
! In  nb_elem_mast     : number of master elements on current zone
! In  nb_elem_slav     : number of slave elements on current zone
! In  list_elem_mast   : name of datastructure for list of master elements on current zone
! In  list_elem_slav   : name of datastructure for list of slave elements on current zone
! IO  nb_pair_zone     : number of contact elements
! IO  list_pair_zone   : list of contact elements
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter ::  debug = ASTER_FALSE
    integer(kind=8) :: nbpatch_t, iret, vali(2)
    integer(kind=8) :: list_pair(nb_elem_mast), li_nb_pt_inte_sl(nb_elem_mast)
    real(kind=8) :: li_pt_inte_sl(nb_elem_mast*16), li_pt_inte_ma(nb_elem_mast*16)
    real(kind=8) :: li_pt_gaus_ma(nb_elem_mast*72)
    integer(kind=8) :: elem_slav_nbnode, cellSlavNume, elem_slav_dime, cellSlavIndx
    integer(kind=8) :: elem_mast_nbnode, cellMastNume, elem_mast_dime, cellMastIndx
    character(len=8) :: elem_mast_code, elem_slav_code
    character(len=8) :: elem_slav_type, elem_mast_type
    real(kind=8) :: elem_mast_coor(3, MT_NNOMAX2D), elem_slav_coor(3, MT_NNOMAX2D)
    integer(kind=8) :: nb_pair, nb_poin_inte
    integer(kind=8) :: iMastNeigh, ISlavStart, iMastStart, iCell
    integer(kind=8) :: iSlavNeigh
    integer(kind=8) :: patch_indx, nb_next_alloc
    real(kind=8) :: inteArea, elem_slav_weight
    real(kind=8) :: poin_inte_sl(32)
    real(kind=8) :: poin_inte_ma(32)
    real(kind=8) :: poin_gauss_ma(74)
    integer(kind=8) ::  elin_mast_nbnode
    integer(kind=8) ::  elin_slav_nbnode
    character(len=8) :: elin_mast_code, elin_slav_code, elem_name
    integer(kind=8) :: nbSlavStart, nbMastPaired, nbMastStart
    integer(kind=8) :: cellMastPaired(nb_elem_mast)
    integer(kind=8) :: elem_start, cellSlavStart(nb_elem_slav), cellMastStart(nb_elem_slav)
    integer(kind=8) :: cellNeighIndx, cellNeighNume
    integer(kind=8) :: slavIndxMini, mastIndxMini, slavIndxMaxi, mastIndxMaxi
    integer(kind=8) :: mast_find_indx
    aster_logical :: l_recup
    integer(kind=8), pointer :: cellMastFlag(:) => null()
    integer(kind=8), pointer :: elem_mast_flag(:) => null()
    integer(kind=8), pointer :: cellSlavFlag(:) => null()
    character(len=8) :: knuzo
    character(len=24) :: sdappa_slne, sdappa_mane
    integer(kind=8), pointer :: meshSlavNeigh(:) => null()
    integer(kind=8), pointer :: meshMastNeigh(:) => null()
    integer(kind=8) :: list_slav_master(4)
    integer(kind=8) :: nbMastNeigh, nbSlavNeigh
    integer(kind=8) :: inteNeigh(4)
    integer(kind=8) :: jv_geom, elem_type_nume
    real(kind=8) :: list_slav_weight(4), weight_test, tole_weight
    character(len=24) :: njv_weight_t, njv_nb_pair_zmpi
    real(kind=8), pointer :: patch_weight_t(:) => null()
    integer(kind=8), pointer :: v_mesh_comapa(:) => null()
    integer(kind=8), pointer :: v_mesh_typmail(:) => null()
    integer(kind=8), pointer :: nb_pair_zmpi(:) => null()
    integer(kind=8), pointer :: list_pair_zmpi(:) => null()
    integer(kind=8), pointer :: li_nbptsl_zmpi(:) => null()
    real(kind=8), pointer :: li_ptintsl_zmpi(:) => null()
    real(kind=8), pointer :: li_ptintma_zmpi(:) => null()
    real(kind=8), pointer :: li_ptgausma_zmpi(:) => null()
    integer(kind=8), pointer :: v_mesh_connex(:) => null()
    integer(kind=8), pointer :: v_connex_lcum(:) => null()
    character(len=16), pointer :: valk(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    inteNeigh(1:4) = 0
    list_slav_master(1:4) = 0
    list_slav_weight(1:4) = 0.d0
    call jelira(mesh//'.PATCH', 'NUTIOC', nbpatch_t)
    nbpatch_t = nbpatch_t-1
    njv_weight_t = sdappa(1:19)//'.PWT '
    call wkvect(njv_weight_t, "V V R", nbpatch_t, vr=patch_weight_t)
    njv_nb_pair_zmpi = sdappa(1:19)//'.NAPP'
    call wkvect(njv_nb_pair_zmpi, "V V I", nb_proc, vi=nb_pair_zmpi)
    call jeexin(sdappa(1:19)//'.MPID', iret)
    if (iret .eq. 0) then
        call wkvect(sdappa(1:19)//'.MPID', 'V V K16', 1, vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
        call wkvect(sdappa(1:19)//'.MPIE', 'V V K16', 1, vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
    else
        call jeveuo(sdappa(1:19)//'.MPID', 'E', vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
        call jeveuo(sdappa(1:19)//'.MPIE', 'E', vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
    end if
    mastIndxMaxi = maxval(list_elem_mast)
    slavIndxMaxi = maxval(list_elem_slav)
    mastIndxMini = minval(list_elem_mast)
    slavIndxMini = minval(list_elem_slav)
!
! - Access to updated geometry
!
    call jeveuo(newgeo(1:19)//'.VALE', 'L', jv_geom)
!
! - Access to mesh
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=v_mesh_typmail)
    call jeveuo(mesh//'.COMAPA', 'L', vi=v_mesh_comapa)
    call jeveuo(mesh//'.CONNEX', 'L', vi=v_mesh_connex)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=v_connex_lcum)
!
! - Objects for flags
!
    AS_ALLOCATE(vi=cellSlavFlag, size=slavIndxMaxi+1-slavIndxMini)
    AS_ALLOCATE(vi=cellMastFlag, size=mastIndxMaxi+1-mastIndxMini)
    AS_ALLOCATE(vi=elem_mast_flag, size=mastIndxMaxi+1-mastIndxMini)
    cellMastPaired(1:nb_elem_mast) = 0
!
! - Object for neighbours (inverse connectivity)
!
    ASSERT(i_zone .le. 100)
    call codent(i_zone-1, 'G', knuzo)
    sdappa_mane = sdappa(1:19)//'.MN'//knuzo(1:2)
    sdappa_slne = sdappa(1:19)//'.EN'//knuzo(1:2)
    call jeveuo(sdappa_mane, 'L', vi=meshMastNeigh)
    call jeveuo(sdappa_slne, 'L', vi=meshSlavNeigh)
!
! - Find initial elements for pairing by PANG method
!
120 continue
    if (debug) then
        write (*, *) 'Recherche mailles de départ'
    end if
    if (pair_method .eq. "RAPIDE") then
        call ap_infast(mesh, newgeo, pair_tole, nb_elem_mast, &
                       list_elem_mast, nb_elem_slav, list_elem_slav, cellSlavFlag, &
                       nbMastStart, cellMastStart, nbSlavStart, cellSlavStart, &
                       sdappa, i_zone)
    elseif (pair_method .eq. "ROBUSTE") then
        call apprin(mesh, newgeo, pair_tole, nb_elem_mast, &
                    list_elem_mast, nb_elem_slav, list_elem_slav, cellSlavFlag, &
                    nbMastStart, cellMastStart, nbSlavStart, cellSlavStart)
    end if
    if (debug) then
        if (nbSlavStart .eq. 0) then
            write (*, *) ". No more slave start element "
        else
            write (*, *) ". Start slave element: ", cellSlavStart(1)
        end if
        if (nbMastStart .eq. 0) then
            write (*, *) ". No more master start element "
        else
            write (*, *) ". Start master element: ", cellMastStart(1)
        end if
    end if
    if (nbSlavStart .eq. 0) then
        goto 110
    end if
!
! - Pairing
!
    if (debug) then
        write (*, *) 'Boucle appariement PANG'
    end if
    do while (nbSlavStart .gt. 0)
!
! ----- Get slave element start
!
        cellSlavNume = cellSlavStart(1)
        cellSlavIndx = cellSlavNume+1-slavIndxMini
        elem_type_nume = v_mesh_typmail(cellSlavNume)
        call jenuno(jexnum('&CATA.TM.NOMTM', elem_type_nume), elem_slav_type)
!
! ----- Shift list of slave element start
!
        do ISlavStart = 1, nbSlavStart-1
            cellSlavStart(ISlavStart) = cellSlavStart(ISlavStart+1)
        end do
        nbSlavStart = nbSlavStart-1
!
! ----- Get current patch
!
        patch_indx = v_mesh_comapa(cellSlavNume)
        if (debug) then
            write (*, *) "Current patch: ", patch_indx
        end if
!
! ----- Get informations about slave element
!
        call aptype(elem_slav_type, &
                    elem_slav_nbnode, elem_slav_code, elem_slav_dime)
!
! ----- Get coordinates of slave element
!
        call apcoor(v_mesh_connex, v_connex_lcum, jv_geom, &
                    cellSlavNume, elem_slav_nbnode, elem_slav_dime, &
                    elem_slav_coor)
        if (debug) then
            write (*, *) "Current slave element: ", cellSlavNume, &
                '(type : ', elem_slav_code, ')'
        end if
!
! ----- Cut element in linearized sub-elements
!
        if (elem_slav_code .eq. "TR6") then
            elin_slav_code = "TR3"
            elin_slav_nbnode = 3
        elseif (elem_slav_code .eq. "QU8" .or. elem_slav_code .eq. "QU9") then
            elin_slav_code = "QU4"
            elin_slav_nbnode = 4
        elseif (elem_slav_code .eq. "SE3") then
            elin_slav_code = "SE2"
            elin_slav_nbnode = 2
        else
            elin_slav_code = elem_slav_code
            elin_slav_nbnode = elem_slav_nbnode
        end if
!
! ----- Compute weight of element
!
        call clpoma(elem_slav_dime, elem_slav_code, elem_slav_coor, elem_slav_nbnode, &
                    elem_slav_weight)
!
! ----- Total weight for patch
!
        patch_weight_t(patch_indx) = patch_weight_t(patch_indx)+elem_slav_weight
!
! ----- Number of neighbours
!
        if (elem_slav_dime .eq. 2) then
            nbSlavNeigh = 2
        elseif (elem_slav_code .eq. 'TR3' .or. &
                elem_slav_code .eq. 'TR6') then
            nbSlavNeigh = 3
        elseif (elem_slav_code .eq. 'QU4' .or. &
                elem_slav_code .eq. 'QU8' .or. &
                elem_slav_code .eq. 'QU9') then
            nbSlavNeigh = 4
        else
            ASSERT(.false.)
        end if
!
! ----- Access to neighbours
!
        if (debug) then
            do iSlavNeigh = 1, nbSlavNeigh
                cellNeighIndx = 4*(cellSlavIndx-1)+iSlavNeigh
                if (cellNeighIndx .ne. 0) then
                    cellNeighNume = meshSlavNeigh(cellNeighIndx)
                    if (cellNeighNume .ne. 0) then
                        elem_name = int_to_char8(cellNeighNume)
                    else
                        elem_name = 'None'
                    end if
                end if
                write (*, *) "Current slave element neighbours: ", elem_name
            end do
        end if
        list_slav_master(1:nbSlavNeigh) = 0
        list_slav_weight(1:4) = 0.d0
!
! ----- Get master element to start
!
        elem_start = cellMastStart(1)
        mast_find_indx = elem_start+1-mastIndxMini

! ----- Shift list of master element start
        do iMastStart = 1, nbMastStart-1
            cellMastStart(iMastStart) = cellMastStart(iMastStart+1)
        end do
        nbMastStart = nbMastStart-1

! ----- Management of list of master elements: first element to seek
        cellMastPaired(1) = elem_start
        nbMastPaired = 1
        cellMastFlag(mast_find_indx) = 1

! ----- Initialization list of contact pairs
        list_pair = 0
        li_pt_inte_sl = 0.0
        li_nb_pt_inte_sl = 0
        nb_pair = 0
        l_recup = .true.

! ----- Loop on master elements
        do while (nbMastPaired .gt. 0)
!
            inteArea = 0.d0

! --------- Get master element
            cellMastNume = cellMastPaired(1)
            cellMastIndx = cellMastNume+1-mastIndxMini
            elem_type_nume = v_mesh_typmail(cellMastNume)
            call jenuno(jexnum('&CATA.TM.NOMTM', elem_type_nume), elem_mast_type)
            if (debug) then
                write (*, *) "Current master element: ", cellMastNume, &
                    '(type : ', elem_mast_type, ')'
            end if

! --------- Access to neighbours
            if (debug) then
                do iMastNeigh = 1, 4
                    cellNeighNume = meshMastNeigh((cellMastIndx-1)*4+iMastNeigh)
                    if (cellNeighNume .ne. 0) then
                        elem_name = int_to_char8(cellNeighNume)
                    else
                        elem_name = 'None'
                    end if
                    write (*, *) "Current master element neighbours: ", elem_name
                end do
            end if

! --------- Shift list of master elements (on supprime de la liste)
            do iCell = 1, nbMastPaired-1
                cellMastPaired(iCell) = cellMastPaired(iCell+1)
            end do
            nbMastPaired = nbMastPaired-1

! --------- Get informations about master element
            call aptype(elem_mast_type, &
                        elem_mast_nbnode, elem_mast_code, elem_mast_dime)

! --------- Get coordinates of master element
            call apcoor(v_mesh_connex, v_connex_lcum, jv_geom, &
                        cellMastNume, elem_mast_nbnode, elem_mast_dime, &
                        elem_mast_coor)

! --------- Cut master element in linearized sub-elements
            if (elem_mast_code .eq. "TR6") then
                elin_mast_code = "TR3"
                elin_mast_nbnode = 3
            elseif (elem_mast_code .eq. "QU8" .or. elem_mast_code .eq. "QU9") then
                elin_mast_code = "QU4"
                elin_mast_nbnode = 4
            elseif (elem_mast_code .eq. "SE3") then
                elin_mast_code = "SE2"
                elin_mast_nbnode = 2
            else
                elin_mast_code = elem_mast_code
                elin_mast_nbnode = elem_mast_nbnode
            end if

! --------- Loop on linearized slave sub-elements
            inteNeigh = 0

! --------- Projection/intersection of elements in slave parametric space
            call prjint(pair_tole, elem_slav_dime, &
                        elin_mast_nbnode, elem_mast_coor, elin_mast_code, &
                        elin_slav_nbnode, elem_slav_coor, elin_slav_code, &
                        poin_inte_sl, inteArea, nb_poin_inte, &
                        inte_neigh_=inteNeigh, ierror_=iret)
            if (iret .eq. 1) then
                vali(1) = cellSlavNume
                vali(2) = cellMastNume
                call utmess('A', 'CONTACT4_6', ni=2, vali=vali)
            end if
            if (debug) then
                write (*, *) "Intersection - Master: ", cellMastNume
                write (*, *) "Intersection - Slave : ", cellSlavNume
                write (*, *) "Intersection - Poids : ", inteArea
                write (*, *) "Intersection - Nb    : ", nb_poin_inte
                write (*, *) "Intersection - Points: ", poin_inte_sl
            end if

! --------- Non-void intersection
            if (inteArea .gt. pair_tole) then
                call aprtpm(pair_tole, elem_slav_dime, &
                            elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                            elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                            poin_inte_sl, nb_poin_inte, poin_inte_ma, &
                            poin_gauss_ma, iret)
            end if

! --------- Add pair
            if (inteArea .gt. pair_tole .and. iret .eq. 0) then
                nb_pair = nb_pair+1
                list_pair(nb_pair) = cellMastNume
                li_nb_pt_inte_sl(nb_pair) = nb_poin_inte
                ASSERT(nb_poin_inte .le. 8)
                li_pt_inte_ma(1+(nb_pair-1)*16:(nb_pair-1)*16+16) = poin_inte_ma(1:16)
                li_pt_gaus_ma(1+(nb_pair-1)*72:(nb_pair-1)*72+72) = poin_gauss_ma(1:72)
                li_pt_inte_sl(1+(nb_pair-1)*16:(nb_pair-1)*16+16) = poin_inte_sl(1:16)
                elem_mast_flag(cellMastIndx) = 1
            end if

! --------- Find neighbour of current master element
            if (inteArea .gt. pair_tole .or. l_recup) then

! ------------- Number of neighbours
                if (elem_mast_code .eq. 'SE2' .or. elem_mast_code .eq. 'SE3') then
                    nbMastNeigh = 2
                    tole_weight = 0.5
                elseif (elem_mast_code .eq. 'TR3' .or. elem_mast_code .eq. 'TR6') then
                    nbMastNeigh = 3
                    tole_weight = 0.05
                elseif (elem_mast_code .eq. 'QU4' .or. elem_mast_code .eq. 'QU8' .or. &
                        elem_mast_code .eq. 'QU9') then
                    nbMastNeigh = 4
                    tole_weight = 0.4
                else
                    ASSERT(.false.)
                end if

! ------------- Prepare next master element
                do iMastNeigh = 1, nbMastNeigh
                    cellNeighIndx = 4*(cellMastIndx-1)+iMastNeigh
                    cellNeighNume = meshMastNeigh(cellNeighIndx)
                    if (cellNeighNume .ne. 0) then
                        if (cellMastFlag(cellNeighNume+1-mastIndxMini) .ne. 1) then
                            nbMastPaired = nbMastPaired+1
                            cellMastPaired(nbMastPaired) = cellNeighNume
                            cellMastFlag(cellNeighNume+1-mastIndxMini) = 1
                        end if
                    end if
                end do

! ------------- Prepare next slave element: higher weight
                do iSlavNeigh = 1, nbSlavNeigh
                    cellNeighIndx = 4*(cellSlavIndx-1)+iSlavNeigh
                    cellNeighNume = meshSlavNeigh(cellNeighIndx)
                    if (cellNeighNume .ne. 0) then
                        if (inteNeigh(iSlavNeigh) == 1 .and. &
                            cellSlavFlag(cellNeighNume+1-slavIndxMini) .ne. 1 .and. &
                            list_slav_weight(iSlavNeigh) .lt. tole_weight) then
                            weight_test = 0.d0
                            call testvois(jv_geom, elem_slav_type, &
                                          elem_mast_coor, elem_mast_code, cellSlavNume, &
                                          pair_tole, weight_test, v_mesh_connex, &
                                          v_connex_lcum)
                            if (weight_test .gt. list_slav_weight(iSlavNeigh) .and. &
                                weight_test .gt. pair_tole) then
                                list_slav_master(iSlavNeigh) = cellMastNume
                                list_slav_weight(iSlavNeigh) = weight_test
                            end if
                        end if
                    end if
                end do
                l_recup = .false.
            end if
        end do
!
! ----- Save pairing informations (contact pair)
!
        if (nb_pair .ne. 0) then
            call apsave_pair(i_zone, cellSlavNume, &
                             nb_pair, list_pair, &
                             li_nb_pt_inte_sl, li_pt_inte_sl, li_pt_inte_ma, &
                             li_pt_gaus_ma, &
                             nb_pair_zmpi(i_proc+1), list_pair_zmpi, &
                             li_nbptsl_zmpi, li_ptintsl_zmpi, li_ptintma_zmpi, &
                             li_ptgausma_zmpi, nb_elem_slav, nb_elem_mast, &
                             nb_next_alloc)
        end if

! ----- Next elements
        if (debug) then
            write (*, *) 'Next elements - Nb: ', nbSlavNeigh
        end if
        do iSlavNeigh = 1, nbSlavNeigh
            cellNeighIndx = 4*(cellSlavIndx-1)+iSlavNeigh
            cellNeighNume = meshSlavNeigh(cellNeighIndx)
            if (cellNeighNume .ne. 0) then
                if (list_slav_master(iSlavNeigh) .ne. 0 .and. &
                    cellSlavFlag(cellNeighNume+1-slavIndxMini) .ne. 1) then
                    nbSlavStart = nbSlavStart+1
                    cellSlavStart(nbSlavStart) = cellNeighNume
                    cellSlavFlag(cellNeighNume+1-slavIndxMini) = 1
                    nbMastStart = nbMastStart+1
                    cellMastStart(nbMastStart) = list_slav_master(iSlavNeigh)
                end if
            end if
        end do
        cellMastFlag(1:mastIndxMaxi+1-mastIndxMini) = 0
    end do
    if (debug) then
        write (*, *) 'Fin boucle appariement PANG'
        write (*, *) 'maille contact trouvee', nb_pair_zmpi
    end if
    goto 120
110 continue
    if (debug) then
        write (*, *) 'Fin appariement PANG RAPIDE'
    end if
!
! - Save values for patch
!
    call apsave_patch(mesh, sdappa, i_zone, &
                      patch_weight_t, nb_proc, list_pair_zmpi, &
                      li_nbptsl_zmpi, li_ptintsl_zmpi, li_ptintma_zmpi, li_ptgausma_zmpi, &
                      nb_pair_zmpi, list_pair_zone, list_nbptit_zone, list_ptitsl_zone, &
                      list_ptitma_zone, li_ptgausma_zone, nb_pair_zone, i_proc)
    !write(*,*)"Fin APSAVE_PATCH", i_proc
    !write(*,*)"NB_pair_zone: ",nb_pair_zone ,i_proc
    !write(*,*)"list_pair_zone: ",list_pair_zone(:) ,i_proc
!
    AS_DEALLOCATE(vi=cellMastFlag)
    AS_DEALLOCATE(vi=cellSlavFlag)
    AS_DEALLOCATE(vi=elem_mast_flag)
    AS_DEALLOCATE(vi=list_pair_zmpi)
    AS_DEALLOCATE(vi=li_nbptsl_zmpi)
    AS_DEALLOCATE(vr=li_ptintsl_zmpi)
    AS_DEALLOCATE(vr=li_ptintma_zmpi)
    AS_DEALLOCATE(vr=li_ptgausma_zmpi)
    call jedetr(njv_weight_t)
    call jedetr(njv_nb_pair_zmpi)
    call jedema()
end subroutine
