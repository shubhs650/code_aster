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
subroutine comp_meca_pvarPt(comporList, comporInfo)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_meca_exc2.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_name.h"
#include "asterfort/comp_ntvariPt.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=16), intent(in) :: comporList(COMPOR_SIZE)
    character(len=19), intent(in) :: comporInfo
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Prepare informations about internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  comporList       : list for parameters of constitutive laws
! In  comporInfo       : object for information about internal state variables and behaviour
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
    integer(kind=8), parameter :: nbCellMesh = 1, nbZoneActi = 1, mapNbZone = 1
    aster_logical, parameter :: l_pmf = ASTER_FALSE
    aster_logical :: l_excl, l_kit_meta, l_cristal, l_kit_thm
    integer(kind=8), pointer :: comporInfoInfo(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    character(len=16), pointer :: comporInfoVari(:) => null()
    character(len=16), pointer :: comporInfoRela(:) => null()
    integer(kind=8) :: nbVari, ntVari, nbVariMaxi, nbVariMeca
    character(len=16) :: postIter, vari_excl, reguVisc, postIncr
    character(len=16) :: relaComp, defoComp, kitComp(4), typeCpla
    character(len=16) :: adrsMGIS
    integer(kind=8) :: solvBehavType
    type(BehaviourPrep_Exte) :: prepExte
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Create list of zones: for each zone (in CARTE), how many elements ?
    call wkvect(comporInfo(1:19)//'.ZONE', 'V V I', 1, vi=comporInfoZone)

! - Number of cells by zone
    comporInfoZone(1) = 1

! - Count total of internal state variables
    call comp_ntvariPt(comporList, comporInfo, &
                       ntVari, nbVariMaxi, prepExte)

    if (ntVari .ne. 0) then
! ----- Create list of comportment information (RELATION, DEFORMATION, etc.)
        call wkvect(comporInfo(1:19)//'.RELA', 'V V K16', 5, vk16=comporInfoRela)

! ----- Create list of internal variables names
        call jecrec(comporInfo(1:19)//'.VARI', 'V V K16', 'NU', 'DISPERSE', 'VARIABLE', 1)
        call jecroc(jexnum(comporInfo(1:19)//'.VARI', 1))

! ----- Get parameters
        relaComp = comporList(RELA_NAME)
        defoComp = comporList(DEFO)
        typeCpla = comporList(PLANESTRESS)
        kitComp(1) = comporList(KIT1_NAME)
        kitComp(2) = comporList(KIT2_NAME)
        kitComp(3) = comporList(KIT3_NAME)
        kitComp(4) = comporList(KIT4_NAME)
        postIter = comporList(POSTITER)
        read (comporList(NVAR), '(I16)') nbVari
        nbVariMeca = 0
        if (comporList(MECA_NVAR) .ne. 'VIDE') then
            read (comporList(MECA_NVAR), '(I16)') nbVariMeca
        end if
        reguVisc = comporList(REGUVISC)
        postIncr = comporList(POSTINCR)

! ----- Detection of specific cases
        call comp_meca_l(relaComp, 'KIT_THM', l_kit_thm)
        call comp_meca_l(relaComp, 'KIT_META', l_kit_meta)
        call comp_meca_l(relaComp, 'CRISTAL', l_cristal)

! ----- Parameters for external constitutive laws
        solvBehavType = prepExte%solvBehavType
        adrsMGIS = prepExte%adrsMGIS

! ----- Exception for name of internal state variables
        call comp_meca_exc2(l_cristal, l_pmf, &
                            l_excl, vari_excl)

! ----- Save names of relation
        comporInfoRela(1) = relaComp
        comporInfoRela(2) = defoComp
        comporInfoRela(3) = typeCpla
        comporInfoRela(4) = reguVisc
        comporInfoRela(5) = postIncr

! ----- Get names of internal state variables
        call jeecra(jexnum(comporInfo(1:19)//'.VARI', 1), 'LONMAX', nbVari)
        call jeveuo(jexnum(comporInfo(1:19)//'.VARI', 1), 'E', vk16=comporInfoVari)
        call comp_meca_name(nbVari, nbVariMeca, &
                            l_excl, vari_excl, l_kit_meta, &
                            relaComp, defoComp, kitComp, typeCpla, postIter, &
                            reguVisc, postIncr, &
                            adrsMGIS, solvBehavType, comporInfoVari)

    end if

! - Save general information
    call wkvect(comporInfo(1:19)//'.INFO', 'V V I', 5, vi=comporInfoInfo)
    comporInfoInfo(1) = nbCellMesh
    comporInfoInfo(2) = mapNbZone
    comporInfoInfo(3) = nbVariMaxi
    comporInfoInfo(4) = ntVari
    comporInfoInfo(5) = nbZoneActi
!
    call jedema()
!
end subroutine
