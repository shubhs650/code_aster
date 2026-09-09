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
subroutine comp_ntvari(modelCell, comporMap, comporInfo, &
                       ntVari, nbVariMaxi, prepExte)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/compGetMecaPart.h"
#include "asterfort/dismoi.h"
#include "asterfort/getExternalBehaviourParaFromAdr.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/teattr.h"
!
    integer(kind=8), pointer :: modelCell(:)
    character(len=19), intent(in) :: comporMap, comporInfo
    integer(kind=8), intent(out) :: ntVari, nbVariMaxi
    type(BehaviourPrep_Exte), pointer :: prepExte(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Count total of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! Ptr modelCell        : pointer to finite element on cells
! In  comporMap        : map for parameters of constitutive laws
! Out ntVari           : total number of internal variables (on all <CARTE> COMPOR)
! Out nbVariMaxi       : maximum number of internal variables on all comportments"
! Ptr prepExte         : pointer to external behaviours parameters
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), pointer :: comporVale(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    integer(kind=8), pointer :: comporDesc(:) => null()
    integer(kind=8), pointer :: comporLima(:) => null()
    integer(kind=8), pointer :: comporLimaCumu(:) => null()
    integer(kind=8) :: nbVale, mapNbCmpMax, nbVari, nbCell, nbCellMesh, mapNbZone
    integer(kind=8) :: iMapZone, iret, iCell, posit
    integer(kind=8) :: affeZoneType, affeZoneNume, elemTypeNume, cellNume
    character(len=16) :: elemTypeName
    character(len=16) :: relaComp, defoComp, kitComp(4), typeCpla, relaMeca
    character(len=16) :: principal, adrsMGIS
    character(len=8) :: mesh
    aster_logical :: l_mfront_cp
!
! --------------------------------------------------------------------------------------------------
!
    ntVari = 0
    nbVariMaxi = 0

! - Access to map
    call jeveuo(comporMap//'.DESC', 'L', vi=comporDesc)
    call jeveuo(comporMap//'.VALE', 'L', vk16=comporVale)
    call jelira(comporMap//'.VALE', 'LONMAX', nbVale)
    call jeveuo(jexnum(comporMap//'.LIMA', 1), 'L', vi=comporLima)
    call jeveuo(jexatr(comporMap//'.LIMA', 'LONCUM'), 'L', vi=comporLimaCumu)
    mapNbZone = comporDesc(3)
    mapNbCmpMax = nbVale/comporDesc(2)
    call dismoi('NOM_MAILLA', comporMap, 'CARTE', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCellMesh)

! - Create list of zones: for each zone (in CARTE), how many elements
    call jeveuo(comporInfo(1:19)//'.ZONE', 'L', vi=comporInfoZone)

! - Count internal variables by comportment
    do iMapZone = 1, mapNbZone

! ----- Get parameters
        relaComp = comporVale(mapNbCmpMax*(iMapZone-1)+RELA_NAME)
        adrsMGIS = comporVale(mapNbCmpMax*(iMapZone-1)+MGIS_ADDR)
        defoComp = comporVale(mapNbCmpMax*(iMapZone-1)+DEFO)
        typeCpla = comporVale(mapNbCmpMax*(iMapZone-1)+PLANESTRESS)
        kitComp(1) = comporVale(mapNbCmpMax*(iMapZone-1)+KIT1_NAME)
        kitComp(2) = comporVale(mapNbCmpMax*(iMapZone-1)+KIT2_NAME)
        kitComp(3) = comporVale(mapNbCmpMax*(iMapZone-1)+KIT3_NAME)
        kitComp(4) = comporVale(mapNbCmpMax*(iMapZone-1)+KIT4_NAME)

! ----- Get mechanical part of behaviour
        relaMeca = 'VIDE'
        call compGetMecaPart(relaComp, kitComp, relaMeca)

! ----- Find right TYPELEM
        affeZoneType = comporDesc(1+3+(iMapZone-1)*2)
        affeZoneNume = comporDesc(1+4+(iMapZone-1)*2)
        if (affeZoneType .eq. 3) then
            nbCell = comporLimaCumu(1+affeZoneNume)-comporLimaCumu(affeZoneNume)
            posit = comporLimaCumu(affeZoneNume)
        elseif (affeZoneType .eq. 1) then
            nbCell = nbCellMesh
            posit = 0
        else
            ASSERT(ASTER_FALSE)
        end if

        do iCell = 1, nbCell

! --------- Get current cell
            if (affeZoneType .eq. 3) then
                cellNume = comporLima(posit+iCell-1)
            elseif (affeZoneType .eq. 1) then
                cellNume = iCell
            elseif (affeZoneType .eq. 0) then
                cellNume = 1
            else
                ASSERT(.false.)
            end if

! --------- Get type of finite element
            if (cellNume .ne. 0 .and. affeZoneType .gt. 0) then
                elemTypeNume = modelCell(cellNume)
                if (elemTypeNume .ne. 0) then
                    call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                    call teattr('C', 'PRINCIPAL', principal, iret, typel=elemTypeName)
                    if (principal .eq. 'OUI') then
                        exit
                    end if
                end if
            end if
        end do

! ----- Get parameters for external programs (MFRONT/UMAT)
        l_mfront_cp = typeCpla .eq. 'ANALYTIQUE'
        call getExternalBehaviourParaFromAdr(elemTypeNume, l_mfront_cp, &
                                             adrsMGIS, relaMeca, defoComp, &
                                             prepExte(iMapZone))

! ----- Get number of internal variables
        read (comporVale(mapNbCmpMax*(iMapZone-1)+NVAR), '(I16)') nbVari

        ntVari = ntVari+nbVari
        nbVariMaxi = max(nbVariMaxi, nbVari)
    end do
!
end subroutine
