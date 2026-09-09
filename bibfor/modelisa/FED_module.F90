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
! ==================================================================================================
!
! Module for the management of Finite Element Descriptor (FED)
!
! ==================================================================================================
!
module FED_module
! ==================================================================================================
    use mesh_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: createFEDFromList
    public :: getAllCellsWithOption, getCellsAffectedByFE
    private :: getElemType
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/adalig.h"
#include "asterfort/alchml.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cormgi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getelem.h"
#include "asterfort/initel.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedupo.h"
#include "asterfort/jedup1.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! createFEDFromList
!
! Create FED (ligrel) from list of cells
!
! In  model             : model
! In  jvBase            : jeveux base to create ligrel
! In  ligrel            : finite element descriptor (ligrel)
! In  nbCell            : number of cells
! Ptr listCell          : list of cells
!
! --------------------------------------------------------------------------------------------------
    subroutine createFEDFromList(modelZ, jvBaseZ, ligrelZ, &
                                 nbCell, listCell)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelZ, ligrelZ, jvBaseZ
        integer(kind=8), intent(in) :: nbCell
        integer(kind=8), pointer :: listCell(:)
! ----- Locals
        character(len=1) :: jvBase
        character(len=8) :: model
        character(len=19) :: modelLigrel, ligrel
        integer(kind=8) :: jvDummy
        integer(kind=8) :: cellNume
        integer(kind=8) :: nbElemAffe, numeElemType, nbElemType, nbElem
        integer(kind=8) :: elemShift, elemTypeMaxi, iElem, iElemType
        integer(kind=8), pointer :: liel(:) => null()
        integer(kind=8), pointer :: listElemAffe(:) => null()
        integer(kind=8), pointer :: listElemType(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ
        ligrel = ligrelZ
        jvBase = jvBaseZ
        ASSERT(nbCell .ge. 1)

! ----- Access to model FED
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! ----- Create .NBNO object
        call wkvect(ligrel(1:19)//'.NBNO', jvBase//' V I', 1, jvDummy)

! ----- Create .LGRF object
        call jedupo(modelLigrel(1:19)//'.LGRF', jvBase, ligrel(1:19)//'.LGRF', .false._1)

! ----- Create .TYFE object
        call jedup1(modelLigrel//'.TYFE', jvBase, ligrel//'.TYFE')

! ----- Get cells affected by a finite element
        call getCellsAffectedByFE(modelLigrel, nbCell, listCell, &
                                  nbElemAffe, listElemAffe)
        if (nbElemAffe .eq. 0) then
            call utmess('F', 'MODELISA4_51')
        end if

! ----- Count number of elements by their type
        call getElemType(modelLigrel, nbElemAffe, listElemAffe, &
                         nbElemType, listElemType, elemTypeMaxi)
        ASSERT(nbElemType .ge. 1)

! ----- Create .LIEL object
        call jecrec(ligrel(1:19)//'.LIEL', jvBase//' V I', 'NU', 'CONTIG', 'VARIABLE', nbElemType)
        call jeecra(ligrel(1:19)//'.LIEL', 'LONT', nbElemType*(elemTypeMaxi+1))
        call jeveuo(ligrel(1:19)//'.LIEL', 'E', vi=liel)

! ----- Fill .LIEL object
        elemShift = 0
        cellNume = 0
        do iElemType = 1, nbElemType
            nbElem = listElemType(nbElemAffe+iElemType)
            numeElemType = listElemType(iElemType)
            call jecroc(jexnum(ligrel(1:19)//'.LIEL', iElemType))
            call jeecra(jexnum(ligrel(1:19)//'.LIEL', iElemType), 'LONMAX', nbElem+1)
            do iElem = 1, nbElem
                elemShift = elemShift+1
                cellNume = cellNume+1
                liel(elemShift) = listElemAffe(cellNume)
            end do
            elemShift = elemShift+1
            liel(elemShift) = numeElemType
        end do

! ----- Finish ligrel
        call adalig(ligrel)
        call cormgi(jvBase, ligrel)
        call initel(ligrel)

! ----- Clean
        AS_DEALLOCATE(vi=listElemAffe)
        AS_DEALLOCATE(vi=listElemType)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCellsAffectedByFE
!
! Get cells affected by a finite element
!
! In  ligrel            : finite element descriptor (ligrel)
! In  nbCell            : number of cells
! Ptr listCell          : list of cells
! Out nbElem            : number of elements (cell with model)
! Ptr listElem          : list of elements (cell with model))
!
! --------------------------------------------------------------------------------------------------
    subroutine getCellsAffectedByFE(ligrelZ, nbCell, listCell, &
                                    nbElem, listElem)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: ligrelZ
        integer(kind=8), intent(in) :: nbCell
        integer(kind=8), pointer :: listCell(:)
        integer(kind=8), intent(out) :: nbElem
        integer(kind=8), pointer :: listElem(:)
! ----- Locals
        character(len=24) :: ligrel
        integer(kind=8) :: grelNume, iCell, cellNume
        integer(kind=8), pointer :: repe(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        ligrel = ligrelZ
        ASSERT(nbCell .ge. 1)

! ----- Access to FED
        call jeveuo(ligrel(1:19)//'.REPE', 'L', vi=repe)

! ----- Create list of elements
        AS_ALLOCATE(vi=listElem, size=nbCell)
        nbElem = 0
        do iCell = 1, nbCell
            cellNume = listCell(iCell)
            if (cellNume .gt. 0) then
                grelNume = repe(1+2*(cellNume-1))
                if (grelNume .ne. 0) then
                    nbElem = nbElem+1
                    listElem(nbElem) = cellNume
                end if
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getElemType
!
! Get type and number of element by type
!
! In  ligrel            : finite element descriptor (ligrel)
! In  nbElem            : number of elements (cell with model)
! Ptr listElem          : list of elements (cell with model)
! Out nbElemType        : number of different element type in list of elements
! Ptr listElemType      : list of different element type
! Out elemTypeMaxi      : maximum number of elements in all types
!
! --------------------------------------------------------------------------------------------------
    subroutine getElemType(ligrelZ, nbElem, listElem, &
                           nbElemType, listElemType, elemTypeMaxi)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: ligrelZ
        integer(kind=8), intent(in) :: nbElem
        integer(kind=8), pointer :: listElem(:), listElemType(:)
        integer(kind=8), intent(out) :: elemTypeMaxi, nbElemType
! ----- Locals
        character(len=24) :: ligrel
        integer(kind=8) :: elemTypeRefe, elemType
        integer(kind=8) :: iElem, iElemType, cellNume, grelNume, numeElemType
        integer(kind=8), pointer :: repe(:) => null()
        integer(kind=8), pointer :: liel(:) => null(), lielCumu(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        ligrel = ligrelZ
        ASSERT(nbElem .ge. 1)
        elemTypeMaxi = 0
        nbElemType = 0

! ----- Access to FED
        call jeveuo(ligrel(1:19)//'.REPE', 'L', vi=repe)
        call jeveuo(jexatr(ligrel(1:19)//'.LIEL', 'LONCUM'), 'L', vi=lielCumu)
        call jeveuo(ligrel(1:19)//'.LIEL', 'L', vi=liel)

! ----- Create list
        AS_ALLOCATE(vi=listElemType, size=2*nbElem)
        elemTypeRefe = 0
        nbElemType = 0
        numeElemType = 0
        do iElem = 1, nbElem
            cellNume = listElem(iElem)
            grelNume = repe(1+2*(cellNume-1))
            ASSERT(grelNume .gt. 0)
            elemType = liel(lielCumu(1+grelNume)-1)
            if (elemType .eq. elemTypeRefe) then
! ------------- Another cell with this type
                listElemType(nbElem+numeElemType) = listElemType(nbElem+numeElemType)+1
            else
! ------------- A new type of element
                nbElemType = nbElemType+1
                numeElemType = nbElemType
                elemTypeRefe = elemType
                listElemType(nbElem+numeElemType) = 1
                listElemType(nbElemType) = elemType
            end if
        end do

! ----- Maximum number of elements for all types
        elemTypeMaxi = 0
        do iElemType = 1, nbElemType
            elemTypeMaxi = max(elemTypeMaxi, listElemType(nbElem+iElemType))
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getAllCellsWithOption
!
! Get all cells with option
!
! In  ligrel            : finite element descriptor (ligrel)
! In  option            : name of option
! In  paramIn           : name of a "representative" input parameter of option
! Out nbCellAffe        : number of cells with this option
! Ptr listCellAffe      : list of cells with this option
!
! --------------------------------------------------------------------------------------------------
    subroutine getAllCellsWithOption(ligrelZ, optionZ, paramInZ, nbCellAffe, listCellAffe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: ligrelZ, optionZ, paramInZ
        integer(kind=8), intent(out) :: nbCellAffe
        integer(kind=8), pointer :: listCellAffe(:)
! ----- Locals
        character(len=24), parameter :: ces = '&&SRLIMA.CES_ELNO', cel = '&&SRLIMA.CEL_ELNO'
        character(len=24) :: ligrel
        integer(kind=8) :: iret, iElem, nbElem
        integer(kind=8) :: jvCesd, jvCesl, iad
!   ------------------------------------------------------------------------------------------------
!
        ligrel = ligrelZ
        nbCellAffe = 0

! ----- Create map
        call alchml(ligrel, optionZ, paramInZ, 'V', cel, iret, ' ')
        ASSERT(iret .eq. 0)
        call celces(cel, 'V', ces)
        call detrsd('CHAMP', cel)
        call jeveuo(ces(1:19)//'.CESD', 'L', jvCesd)
        call jeveuo(ces(1:19)//'.CESL', 'L', jvCesl)

! ----- Create list
        nbElem = zi(jvCesd-1+1)
        AS_ALLOCATE(vi=listCellAffe, size=nbElem)
        do iElem = 1, nbElem
            call cesexi('C', jvCesd, jvCesl, iElem, 1, 1, 1, iad)
            if (iad .gt. 0) then
                nbCellAffe = nbCellAffe+1
                listCellAffe(nbCellAffe) = iElem
            end if
        end do
        call detrsd('CHAM_ELEM_S', ces)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module FED_module
