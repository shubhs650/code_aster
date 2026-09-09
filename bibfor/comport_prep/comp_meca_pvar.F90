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
subroutine comp_meca_pvar(modelFED, comporMap, comporInfo)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_meca_exc2.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_name.h"
#include "asterfort/comp_ntvari.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lteatt.h"
#include "asterfort/wkvect.h"
!
    character(len=19), intent(in) :: modelFED, comporMap, comporInfo
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Prepare informations about internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  modelFED           : FED from model
! In  comporMap          : map for parameters of constitutive laws
! In  comporInfo         : object for information about internal state variables and behaviour
!
!    ComporInfo:
!       INFO.INFO = global parameters
!         comporInfoInfo(1) = nbCellMesh
!          => total number of elements in mesh
!         comporInfoInfo(2) = mapNbZone
!          => total number of zone in CARTE
!         comporInfoInfo(3) = nbVariMaxi
!          => maximum number of internal variables
!         comporInfoInfo(4) = ntVari
!          => total number of internal variables
!       INFO.VARI = Collection of mapNbZone (from CARTE) x Vecteur_Info
!       For each zone   : Vector_Info is list of nbVari name of internal variables (K16)
!       INFO.ZONE = list on mapNbZone (from CARTE)
!       For each zone   : number of elements with this comportement
!       INFO.RELA = list on mapNbZone (from CARTE) * 8
!       For each zone   : some information from comprotement (name of RELATION, DEFORMATION, ...)
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_excl, l_kit_meta, l_cristal, l_pmf, l_kit_thm
    aster_logical :: l_zone_read
    character(len=8) :: mesh
    integer(kind=8), pointer :: comporInfoInfo(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    integer(kind=8), pointer :: zoneRead(:) => null()
    integer(kind=8), pointer :: modelCell(:) => null()
    character(len=16), pointer :: comporInfoVari(:) => null()
    character(len=16), pointer :: comporInfoRela(:) => null()
    character(len=16), pointer :: comporVale(:) => null()
    integer(kind=8), pointer :: comporDesc(:) => null()
    integer(kind=8), pointer :: comporPtma(:) => null()
    integer(kind=8) :: nbVale, mapNbCmpMax, mapNbZone, nbVari, ntVari, nbVariMaxi, nbZoneActi
    integer(kind=8) :: mapZoneNume, iCellMesh, nbCellMesh, iret, elemTypeNume, nbVariMeca
    character(len=16) :: postIter, variExcl, reguVisc, postIncr
    character(len=16) :: relaComp, defoComp, kitComp(4), typeCpla
    character(len=16) :: adrsMGIS, elemTypeName
    integer(kind=8) :: solvBehavType
    type(BehaviourPrep_Exte), pointer :: prepExte(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    nbZoneActi = 0

! - Access to list of finite elements
    call jeveuo(modelFED//'.TYFE', 'L', vi=modelCell)

! - Access to mesh
    call dismoi('NOM_MAILLA', comporMap, 'CARTE', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCellMesh)

! - Access to map
    call jeveuo(comporMap//'.DESC', 'L', vi=comporDesc)
    call jeveuo(comporMap//'.VALE', 'L', vk16=comporVale)
    call jelira(comporMap//'.VALE', 'LONMAX', nbVale)
    mapNbZone = comporDesc(3)
    mapNbCmpMax = nbVale/comporDesc(2)
    call etenca(comporMap, modelFED, iret)
    call jeveuo(comporMap//'.PTMA', 'L', vi=comporPtma)

! - Create list of zones: for each zone (in CARTE), how many elements ?
    call wkvect(comporInfo(1:19)//'.ZONE', 'V V I', mapNbZone, vi=comporInfoZone)

! - Count number of cells by zone
    do iCellMesh = 1, nbCellMesh
        mapZoneNume = comporPtma(iCellMesh)
        if (mapZoneNume .ne. 0 .and. modelCell(iCellMesh) .ne. 0) then
            comporInfoZone(mapZoneNume) = comporInfoZone(mapZoneNume)+1
        end if
    end do

! - Prepare objects for external constitutive laws
    allocate (prepExte(mapNbZone))

! - Count total of internal state variables
    call comp_ntvari(modelCell, comporMap, comporInfo, &
                     ntVari, nbVariMaxi, prepExte)
    AS_ALLOCATE(vi=zoneRead, size=mapNbZone)

    if (ntVari .ne. 0) then

! ----- Create list of comportment information (RELATION, DEFORMATION, etc.)
        call wkvect(comporInfo(1:19)//'.RELA', 'V V K16', 5*mapNbZone, vk16=comporInfoRela)

! ----- Create list of internal variables names
        call jecrec(comporInfo(1:19)//'.VARI', 'V V K16', 'NU', 'DISPERSE', 'VARIABLE', mapNbZone)
        do mapZoneNume = 1, mapNbZone
            call jecroc(jexnum(comporInfo(1:19)//'.VARI', mapZoneNume))
        end do

        do iCellMesh = 1, nbCellMesh
! --------- Get current zone
            mapZoneNume = comporPtma(iCellMesh)
            if (mapZoneNume .eq. 0) then
                l_zone_read = ASTER_TRUE
            else
                ASSERT(mapZoneNume .ne. 0)
                l_zone_read = zoneRead(mapZoneNume) .eq. 1
            end if

            if (.not. l_zone_read) then
! ------------- Get parameters
                relaComp = comporVale(mapNbCmpMax*(mapZoneNume-1)+RELA_NAME)
                defoComp = comporVale(mapNbCmpMax*(mapZoneNume-1)+DEFO)
                typeCpla = comporVale(mapNbCmpMax*(mapZoneNume-1)+PLANESTRESS)
                kitComp(1) = comporVale(mapNbCmpMax*(mapZoneNume-1)+KIT1_NAME)
                kitComp(2) = comporVale(mapNbCmpMax*(mapZoneNume-1)+KIT2_NAME)
                kitComp(3) = comporVale(mapNbCmpMax*(mapZoneNume-1)+KIT3_NAME)
                kitComp(4) = comporVale(mapNbCmpMax*(mapZoneNume-1)+KIT4_NAME)
                postIter = comporVale(mapNbCmpMax*(mapZoneNume-1)+POSTITER)
                read (comporVale(mapNbCmpMax*(mapZoneNume-1)+NVAR), '(I16)') nbVari
                nbVariMeca = 0
                if (comporVale(mapNbCmpMax*(mapZoneNume-1)+MECA_NVAR) .ne. 'VIDE') then
                    read (comporVale(mapNbCmpMax*(mapZoneNume-1)+MECA_NVAR), '(I16)') nbVariMeca
                end if
                reguVisc = comporVale(mapNbCmpMax*(mapZoneNume-1)+REGUVISC)
                postIncr = comporVale(mapNbCmpMax*(mapZoneNume-1)+POSTINCR)

! ------------- Detection of specific cases
                call comp_meca_l(relaComp, 'KIT_THM', l_kit_thm)
                call comp_meca_l(relaComp, 'KIT_META', l_kit_meta)
                call comp_meca_l(relaComp, 'CRISTAL', l_cristal)
                l_pmf = ASTER_FALSE
                elemTypeNume = modelCell(iCellMesh)
                if (elemTypeNume .ne. 0) then
                    call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                    l_pmf = lteatt('TYPMOD2', 'PMF', typel=elemTypeName)
                end if

! ------------- Parameters for external constitutive laws
                solvBehavType = prepExte(mapZoneNume)%solvBehavType
                adrsMGIS = prepExte(mapZoneNume)%adrsMGIS

! ------------- Exception for name of internal state variables
                call comp_meca_exc2(l_cristal, l_pmf, &
                                    l_excl, variExcl)

! ------------- Save names of relation
                comporInfoRela(5*(mapZoneNume-1)+1) = relaComp
                comporInfoRela(5*(mapZoneNume-1)+2) = defoComp
                comporInfoRela(5*(mapZoneNume-1)+3) = typeCpla
                comporInfoRela(5*(mapZoneNume-1)+4) = reguVisc
                comporInfoRela(5*(mapZoneNume-1)+5) = postIncr

! ------------- Get names of internal state variables
                call jeecra(jexnum(comporInfo(1:19)//'.VARI', mapZoneNume), 'LONMAX', nbVari)
                call jeveuo(jexnum(comporInfo(1:19)//'.VARI', mapZoneNume), 'E', &
                            vk16=comporInfoVari)
                call comp_meca_name(nbVari, nbVariMeca, &
                                    l_excl, variExcl, l_kit_meta, &
                                    relaComp, defoComp, kitComp, typeCpla, postIter, &
                                    reguVisc, postIncr, &
                                    adrsMGIS, solvBehavType, comporInfoVari)

! ------------- Save current zone
                zoneRead(mapZoneNume) = 1
                nbZoneActi = nbZoneActi+1
            end if
        end do
!
    end if

! - Save general information
    call wkvect(comporInfo(1:19)//'.INFO', 'V V I', 5, vi=comporInfoInfo)
    comporInfoInfo(1) = nbCellMesh
    comporInfoInfo(2) = mapNbZone
    comporInfoInfo(3) = nbVariMaxi
    comporInfoInfo(4) = ntVari
    comporInfoInfo(5) = nbZoneActi
!
    deallocate (prepExte)
    AS_DEALLOCATE(vi=zoneRead)
!
    call jedema()
!
end subroutine
