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
subroutine char_crea_cart(phenom, loadType, load, mesh, valeType, &
                          nbMap, map, nbCmp, &
                          createMap_, physQuantity_, cmpName_)
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "jeveux.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisnnl.h"
#include "asterfort/nocart.h"
!
    character(len=*), intent(in) :: phenom
    character(len=16), intent(in) :: loadType
    character(len=8), intent(in) :: load, mesh
    character(len=4), intent(in) :: valeType
    integer(kind=8), intent(out) :: nbMap
    character(len=19), intent(out) :: map(LOAD_MAP_NBMAX)
    integer(kind=8), intent(out) :: nbCmp(LOAD_MAP_NBMAX)
    aster_logical, optional, intent(in) :: createMap_
    character(len=8), optional, intent(out) :: physQuantity_(LOAD_MAP_NBMAX)
    character(len=8), optional, intent(out) :: cmpName_(LOAD_MAP_NBMAX, LOAD_MAP_NBCMPMAX)
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Creation and initialization to zero of <CARTE> for Neumann loads
!
! --------------------------------------------------------------------------------------------------
!
! In  phenom       : phenomenon (MECANIQUE/THERMIQUE/ACOUSTIQUE)
! In  loadType     : type of load
! In  mesh         : name of mesh
! In  load         : name of load
! In  valeType     : affected value type (real, complex or function)
! Out nbMap        : number of <CARTE> for this Neumann load
! Out map          : <CARTE> for this Neumann load
! Out nbCmp        : number of components for each map
! In  createMap    : flag to create maps
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: cmpName(LOAD_MAP_NBMAX, LOAD_MAP_NBCMPMAX)
    character(len=13) :: obje_pref
    character(len=8) :: physQuantity(LOAD_MAP_NBMAX)
    character(len=4) :: mapType(LOAD_MAP_NBMAX)
    integer(kind=8) :: jvValv, iMap, i_cmp, iret
    aster_logical :: l_init(LOAD_MAP_NBMAX), createMap
    character(len=8), pointer :: ncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    nbMap = 0
    map = ' '
    nbCmp = 0
    cmpName = ' '
    physQuantity = ' '
    mapType = ' '
    l_init = ASTER_FALSE

! - Create or not the map ?
    createMap = ASTER_TRUE
    if (present(createMap_)) then
        createMap = createMap_
    end if

! - Prefix of <CARTE> objects
    call lisnnl(phenom, load, obje_pref)

! - Number of <CARTE> objects
    if (loadType .eq. 'EFFE_FOND') then
        nbMap = 2
    else if (loadType .eq. 'ONDE_PLANE') then
        nbMap = 2
    else if (loadType .eq. 'ROTATION') then
        nbMap = 1
    else if (loadType .eq. 'ONDE_FLUI') then
        nbMap = 1
    else if (loadType .eq. 'ECHANGE_THM') then
        nbMap = 1
    else if (loadType .eq. 'ECHANGE_THM_HR') then
        nbMap = 1
    else if (loadType .eq. 'PRES_REP') then
        nbMap = 1
    else if (loadType .eq. 'FORCE_TUYAU') then
        nbMap = 1
    else if (loadType .eq. 'PESANTEUR') then
        nbMap = 1
    else if (loadType .eq. 'VITE_FACE') then
        nbMap = 1
    else
        ASSERT(ASTER_FALSE)
    end if
    ASSERT(nbMap .le. LOAD_MAP_NBMAX)

! - Name of the <CARTE>
    if (loadType .eq. 'EFFE_FOND') then
        map(1) = obje_pref(1:13)//'.EFOND'
        map(2) = obje_pref(1:13)//'.PREFF'
    else if (loadType .eq. 'ONDE_PLANE') then
        map(1) = obje_pref(1:13)//'.ONDPL'
        map(2) = obje_pref(1:13)//'.ONDPR'
    else if (loadType .eq. 'ROTATION') then
        map(1) = obje_pref(1:13)//'.ROTAT'
    else if (loadType .eq. 'ONDE_FLUI') then
        map(1) = obje_pref(1:13)//'.ONDE'
    else if (loadType .eq. 'ECHANGE_THM') then
        map(1) = obje_pref(1:13)//'.ETHM'
    else if (loadType .eq. 'ECHANGE_THM_HR') then
        map(1) = obje_pref(1:13)//'.ETHMH'
    else if (loadType .eq. 'PRES_REP') then
        map(1) = obje_pref(1:13)//'.PRESS'
    else if (loadType .eq. 'FORCE_TUYAU') then
        map(1) = obje_pref(1:13)//'.PRESS'
    else if (loadType .eq. 'PESANTEUR') then
        map(1) = obje_pref(1:13)//'.PESAN'
    else if (loadType .eq. 'VITE_FACE') then
        map(1) = obje_pref(1:13)//'.VFACE'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Name of the <GRANDEUR>
    if (loadType .eq. 'EFFE_FOND') then
        if (valeType .eq. 'REEL') then
            physQuantity(1) = 'NEUT_R'
            physQuantity(2) = 'PRES_R'
        else if (valeType .eq. 'FONC') then
            physQuantity(1) = 'NEUT_R'
            physQuantity(2) = 'PRES_F'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ONDE_PLANE') then
        if (valeType .eq. 'FONC') then
            physQuantity(1) = 'NEUT_K8'
            physQuantity(2) = 'NEUT_R'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ROTATION') then
        if (valeType .eq. 'REEL') then
            physQuantity(1) = 'ROTA_R'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ONDE_FLUI') then
        if (valeType .eq. 'REEL') then
            physQuantity(1) = 'ONDE_R'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ECHANGE_THM') then
        if (valeType .eq. 'REEL') then
            physQuantity(1) = 'ETHM_R'
        else if (valeType .eq. 'FONC') then
            physQuantity(1) = 'ETHM_F'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ECHANGE_THM_HR') then
        if (valeType .eq. 'REEL') then
            physQuantity(1) = 'ETHMH_R'
        else if (valeType .eq. 'FONC') then
            physQuantity(1) = 'ETHMH_F'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'PRES_REP') then
        if (valeType .eq. 'REEL') then
            physQuantity(1) = 'PRES_R'
        else if (valeType .eq. 'FONC') then
            physQuantity(1) = 'PRES_F'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'FORCE_TUYAU') then
        if (valeType .eq. 'REEL') then
            physQuantity(1) = 'PRES_R'
        else if (valeType .eq. 'FONC') then
            physQuantity(1) = 'PRES_F'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'PESANTEUR') then
        if (valeType .eq. 'REEL') then
            physQuantity(1) = 'PESA_R'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'VITE_FACE') then
        if (valeType .eq. 'COMP') then
            physQuantity(1) = 'VFAC_C'
        elseif (valeType .eq. 'REEL') then
            physQuantity(1) = 'VFAC_R'
        elseif (valeType .eq. 'FONC') then
            physQuantity(1) = 'VFAC_F'
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if

! - Type of the <CARTE>
    if (loadType .eq. 'EFFE_FOND') then
        if (valeType .eq. 'REEL') then
            mapType(1) = 'R'
            mapType(2) = 'R'
        else if (valeType .eq. 'FONC') then
            mapType(1) = 'R'
            mapType(2) = 'K8'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ONDE_PLANE') then
        if (valeType .eq. 'FONC') then
            mapType(1) = 'K8'
            mapType(2) = 'R'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ROTATION') then
        if (valeType .eq. 'REEL') then
            mapType(1) = 'R'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ONDE_FLUI') then
        if (valeType .eq. 'REEL') then
            mapType(1) = 'R'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ECHANGE_THM') then
        if (valeType .eq. 'REEL') then
            mapType(1) = 'R'
        else if (valeType .eq. 'FONC') then
            mapType(1) = 'K8'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'ECHANGE_THM_HR') then
        if (valeType .eq. 'REEL') then
            mapType(1) = 'R'
        else if (valeType .eq. 'FONC') then
            mapType(1) = 'K8'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'PRES_REP') then
        if (valeType .eq. 'REEL') then
            mapType(1) = 'R'
        else if (valeType .eq. 'FONC') then
            mapType(1) = 'K8'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'FORCE_TUYAU') then
        if (valeType .eq. 'REEL') then
            mapType(1) = 'R'
        else if (valeType .eq. 'FONC') then
            mapType(1) = 'K8'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'PESANTEUR') then
        if (valeType .eq. 'REEL') then
            mapType(1) = 'R'
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (loadType .eq. 'VITE_FACE') then
        if (valeType .eq. 'COMP') then
            mapType(1) = 'C'
        elseif (valeType .eq. 'REEL') then
            mapType(1) = 'R'
        elseif (valeType .eq. 'FONC') then
            mapType(1) = 'K8'
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if

! - Components of the <CARTE>
    if (loadType .eq. 'EFFE_FOND') then
        nbCmp(1) = 1
        nbCmp(2) = 1
        cmpName(1, 1) = 'X1'
        cmpName(2, 1) = 'PRES'
    else if (loadType .eq. 'ONDE_PLANE') then
        nbCmp(1) = 2
        nbCmp(2) = 10
        cmpName(1, 1) = 'Z1'
        cmpName(1, 2) = 'Z2'
        cmpName(2, 1) = 'X1'
        cmpName(2, 2) = 'X2'
        cmpName(2, 3) = 'X3'
        cmpName(2, 4) = 'X4'
        cmpName(2, 5) = 'X5'
        cmpName(2, 6) = 'X6'
        cmpName(2, 7) = 'X7'
        cmpName(2, 8) = 'X8'
        cmpName(2, 9) = 'X9'
        cmpName(2, 10) = 'X10'
    else if (loadType .eq. 'ROTATION') then
        nbCmp(1) = 7
        cmpName(1, 1) = 'OME'
        cmpName(1, 2) = 'AR'
        cmpName(1, 3) = 'BR'
        cmpName(1, 4) = 'CR'
        cmpName(1, 5) = 'X'
        cmpName(1, 6) = 'Y'
        cmpName(1, 7) = 'Z'
    else if (loadType .eq. 'ONDE_FLUI') then
        nbCmp(1) = 1
        cmpName(1, 1) = 'PRES'
    else if (loadType .eq. 'ECHANGE_THM') then
        nbCmp(1) = 6
        cmpName(1, 1) = 'COEF1'
        cmpName(1, 2) = 'COEF2'
        cmpName(1, 3) = 'COEF3'
        cmpName(1, 4) = 'COEF4'
        cmpName(1, 5) = 'PRE1'
        cmpName(1, 6) = 'PRE2'
    else if (loadType .eq. 'ECHANGE_THM_HR') then
        nbCmp(1) = 3
        cmpName(1, 1) = 'COEF1'
        cmpName(1, 2) = 'COEF2'
        cmpName(1, 3) = 'HR1'
    else if (loadType .eq. 'PRES_REP') then
        nbCmp(1) = 2
        cmpName(1, 1) = 'PRES'
        cmpName(1, 2) = 'CISA'
    else if (loadType .eq. 'FORCE_TUYAU') then
        nbCmp(1) = 2
        cmpName(1, 1) = 'PRES'
        cmpName(1, 2) = 'CISA'
    else if (loadType .eq. 'PESANTEUR') then
        nbCmp(1) = 4
        cmpName(1, 1) = 'G'
        cmpName(1, 2) = 'AG'
        cmpName(1, 3) = 'BG'
        cmpName(1, 4) = 'CG'
    else if (loadType .eq. 'VITE_FACE') then
        nbCmp(1) = 5
        cmpName(1, 1) = 'VITE'
        cmpName(1, 2) = 'INDC'
        cmpName(1, 3) = 'DIRX'
        cmpName(1, 4) = 'DIRY'
        cmpName(1, 5) = 'DIRZ'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Creation of the <CARTE>
    if (createMap) then
        do iMap = 1, nbMap
            call exisd('CARTE', map(iMap), iret)
            if (iret .eq. 0) then
                call alcart('G', map(iMap), mesh, physQuantity(iMap))
                l_init(iMap) = ASTER_TRUE
            else
                l_init(iMap) = ASTER_FALSE
            end if
        end do
    end if
! - Initialization of the <CARTE>
    if (createMap) then
        do iMap = 1, nbMap
            if (l_init(iMap)) then
                call jeveuo(map(iMap)//'.NCMP', 'E', vk8=ncmp)
                call jeveuo(map(iMap)//'.VALV', 'E', jvValv)
                do i_cmp = 1, nbCmp(iMap)
                    ncmp(i_cmp) = cmpName(iMap, i_cmp)
                    if (mapType(iMap) .eq. 'R') then
                        zr(jvValv-1+i_cmp) = 0.d0
                    elseif (mapType(iMap) .eq. 'C') then
                        zc(jvValv-1+i_cmp) = (0.d0, 0.d0)
                    else if (mapType(iMap) .eq. 'K8') then
                        zk8(jvValv-1+i_cmp) = '&FOZERO'
                    else if (mapType(iMap) .eq. 'K16') then
                        zk16(jvValv-1+i_cmp) = ' '
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                end do
                call nocart(map(iMap), 1, nbCmp(iMap))
            end if
        end do
    end if
!
    if (present(physQuantity_)) then
        physQuantity_ = physQuantity
    end if
    if (present(cmpName_)) then
        cmpName_ = cmpName
    end if
!
    call jedema()
end subroutine
