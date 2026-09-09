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
subroutine peeint(tableOut, model, nbocc)
!
    use MGIS_module
    implicit none
!
#include "asterc/indik8.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/convertFieldNodeToNeutElem.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismlg.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim1.h"
#include "asterfort/getelem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/peecal.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsGetOneBehaviourFromResult.h"
#include "asterfort/rsSelectStoringIndex.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/umalma.h"
#include "asterfort/utflmd.h"
#include "asterfort/utmess.h"
#include "asterfort/varinonu.h"
#include "jeveux.h"
!
    integer(kind=8) :: nbocc
    character(len=8) :: model
    character(len=19) :: tableOut
!
! --------------------------------------------------------------------------------------------------
!
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "INTEGRALE"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbParaResult = 4, nbParaField = 2, nbCmpOk = 6
    character(len=8), parameter :: paraTypeResult(nbParaResult) = (/'K16', 'I  ', 'R  ', 'R  '/)
    character(len=8), parameter :: paraTypeField(nbParaField) = (/'K16', 'R  '/)
    character(len=16), parameter :: paraNameResult(nbParaResult) = (/'NOM_CHAM  ', 'NUME_ORDRE', &
                                                                     'INST      ', 'VOL       '/)
    character(len=16), parameter :: paraNameField(nbParaField) = (/'CHAM_GD ', 'VOL     '/)
    character(len=3), parameter :: cmpNameOk(nbCmpOk) = (/'N  ', 'VY ', 'VZ ', 'MT ', 'MFY', 'MFZ'/)
    integer(kind=8) :: iret, ibid, iocc, nbret
    integer(kind=8) :: cmpNume, numeStore
    integer(kind=8) :: iCmp, iCellCompute, iGroup, iStore, iCell, iCmpOk
    integer(kind=8) :: nbCellMesh, nbCellUser, nbCellFilter, nbCellCompute
    integer(kind=8) :: nbCmp, nbStore, nbCmpField, nbGroup, nbVari
    integer(kind=8) :: pdtElemType
    real(kind=8) :: inst
    character(len=8) :: mesh, resultIn, physName
    character(len=4) :: fieldSupp, lStructElem
    character(len=8), parameter :: locaNameAll = 'TOUT', locaNameGroup = 'GROUP_MA'
    character(len=24), parameter :: locaNameUnion = 'UNION_GROUP_MA'
    character(len=24), parameter :: factorKeyword = 'INTEGRALE'
    character(len=24), parameter :: listCellUser = '&&PEEINT.TYFE_USER'
    character(len=24) :: listCellFilter
    character(len=24) :: numeStoreJv, timeStoreJv, compor
    character(len=19) :: field, fieldFromUser
    character(len=19) :: ligrel, modelligrel
    character(len=19), parameter :: cespoi = '&&PEEINT.CESPOI'
    character(len=19) :: fieldInput
    character(len=24) :: fieldName, groupName
    aster_logical :: convToNeut, lFromField, lFromResult, lVariName, lCmpOk, l_pmesh, l_empty
    aster_logical :: l_newlig
    integer(kind=8) :: filterTypeNume
    character(len=8) :: filterTypeName
    character(len=8), pointer :: cmpNameNode(:) => null(), cmpNameNeut(:) => null()
    character(len=8), pointer :: cmpNameInit(:) => null()
    integer(kind=8), pointer :: cellCompute(:) => null()
    integer(kind=8), pointer :: cellFilter(:) => null()
    integer(kind=8), pointer :: listNumeStore(:) => null()
    integer(kind=8), pointer :: listElemType(:) => null()
    real(kind=8), pointer :: listTimeStore(:) => null()
    character(len=16), pointer :: variName(:) => null()
    character(len=8), pointer :: cmpName(:) => null()
    character(len=8), pointer :: cmpNameAll(:) => null()
    character(len=24), pointer :: groupCell(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Main parameters
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCellMesh)
    l_pmesh = isParallelMesh(mesh)

! - Origin of fields
    call getvid(' ', 'RESULTAT', scal=resultIn, nbret=nbret)
    lFromResult = nbret .ne. 0
    call getvid(' ', 'CHAM_GD', scal=fieldFromUser, nbret=nbret)
    lFromField = nbret .ne. 0
    if (lFromField) then
        ASSERT(.not. lFromResult)
        fieldInput = 'TMP_CHAMP_GD'
        call copisd('CHAMP', 'V', fieldFromUser, fieldInput)
    end if

! - Select storing index from user
    call rsSelectStoringIndex(resultIn, lFromField, &
                              nbStore, numeStoreJv, timeStoreJv)
    call jeveuo(numeStoreJv, 'L', vi=listNumeStore)
    call jeveuo(timeStoreJv, 'L', vr=listTimeStore)

! - Create table
    call tbcrsd(tableOut, 'G')
    if (lFromResult) then
        call tbajpa(tableOut, nbParaResult, paraNameResult, paraTypeResult)
    else
        call tbajpa(tableOut, nbParaField, paraNameField, paraTypeField)
    end if
!
    do iocc = 1, nbocc
! ----- Get list of cells from user to create reduced domain
        call getelem(mesh, factorKeyword, iocc, 'F', &
                     listCellUser, nbCellUser, l_keep_propz=ASTER_TRUE, model=model)
        call jeexin(listCellUser, iret)
        l_empty = ASTER_FALSE
        if (iret .eq. 0) then
            ASSERT(l_pmesh)
            l_empty = ASTER_TRUE
        end if

! ----- Sort with topological dimension of cells
        call getvtx(factorKeyword, 'TYPE_MAILLE', iocc=iocc, scal=filterTypeName, nbret=iret)
        if (iret .eq. 0) then
            listCellFilter = listCellUser
            nbCellFilter = nbCellUser
        else
            if (filterTypeName .eq. '1D') then
                filterTypeNume = 1
            else if (filterTypeName .eq. '2D') then
                filterTypeNume = 2
            else if (filterTypeName .eq. '3D') then
                filterTypeNume = 3
            else
                ASSERT(ASTER_FALSE)
            end if
            listCellFilter = '&&PEEINT.MAILLES_FILTRE'
            nbCellFilter = 0
            if (.not. l_empty) then
                call utflmd(mesh, listCellUser, nbCellUser, filterTypeNume, ' ', &
                            nbCellFilter, listCellFilter)
            end if
            if (nbCellFilter .eq. 0) then
                if (.not. l_pmesh) call utmess('F', 'PREPOST2_8')
            else
                call utmess('I', 'PREPOST2_7', si=(nbCellUser-nbCellFilter))
            end if
        end if

! ----- Create LIGREL
        l_newlig = ASTER_FALSE
        if (.not. l_empty) then
            ligrel = '&&PEEINT.LIGREL'
            call jeveuo(listCellFilter, 'L', vi=cellFilter)
            call exlim1(cellFilter, nbCellFilter, model, 'V', ligrel)
            l_newlig = ASTER_TRUE
        else
            call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
        end if

! ----- Get name of components
        call getvtx(factorKeyword, 'NOM_CMP', iocc=iocc, nbval=0, nbret=nbCmp)
        nbCmp = -nbCmp
        lVariName = ASTER_FALSE

        if (nbCmp .eq. 0) then
            ASSERT(.not. l_pmesh)
! --------- Get list for internal state variables
            lVariName = ASTER_TRUE
            call getvtx(factorKeyword, 'NOM_VARI', iocc=iocc, nbval=0, nbret=nbVari)
            nbVari = -nbVari
            ASSERT(nbVari .gt. 0)
            if (nbVari .gt. 0) then
                if (.not. lFromResult) then
                    call utmess('F', 'POSTELEM_6')
                end if
            end if
            AS_ALLOCATE(vk16=variName, size=nbVari)
            call getvtx(factorKeyword, 'NOM_VARI', iocc=iocc, nbval=nbVari, vect=variName)

! --------- Get behaviour (only one !)
            if (lFromResult) then
                call rsGetOneBehaviourFromResult(resultIn, nbStore, listNumeStore, compor)
                if (compor .eq. '#SANS') then
                    call utmess('F', 'RESULT1_5')
                end if
                if (compor .eq. '#PLUSIEURS') then
                    call utmess('F', 'RESULT1_6')
                end if
            else
                call utmess('F', 'RESULT1_7')
            end if

! --------- Get name of internal state variables
            nbCmp = nbVari
            AS_ALLOCATE(vk8=cmpName, size=nbCmp)
            AS_ALLOCATE(vk8=cmpNameAll, size=nbCellFilter*nbCmp)
            if (hasMFront(compor)) then
                call utmess('F', "COMPOR6_6")
            end if
            call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelligrel)
            call varinonu(modelligrel, compor, &
                          nbCellFilter, cellFilter, &
                          nbVari, variName, cmpNameAll)
        else
            AS_ALLOCATE(vk8=cmpName, size=nbCmp)
            call getvtx(factorKeyword, 'NOM_CMP', iocc=iocc, nbval=nbCmp, vect=cmpName, nbret=iret)
        end if

! ----- Copy name of components
        AS_ALLOCATE(vk8=cmpNameInit, size=nbCmp)
        do iCmp = 1, nbCmp
            if (lVariName) then
                cmpNameInit(iCmp) = variName(iCmp) (1:8)
                cmpName(icmp) = cmpNameAll((icmp-1)*nbCellFilter+1)
            else
                cmpNameInit(iCmp) = cmpName(iCmp)
            end if
        end do

! ----- No structural elements !
!       Except POU_D_T, components N, VY, VZ, MFT, MFY, MFZ
        if (.not. l_empty) then
            call dismlg('EXI_RDM', ligrel, ibid, lStructElem, iret)
        else
            lStructElem = "NON"
        end if
        if (lStructElem .eq. 'OUI') then
            call jenonu(jexnom('&CATA.TE.NOMTE', 'MECA_POU_D_T'), pdtElemType)
            call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelligrel)
            call jeveuo(modelligrel//'.TYFE', 'L', vi=listElemType)
!           Check components
            do iCmp = 1, nbCmp
                lCmpOk = .false.
                do iCmpOk = 1, nbCmpOk
                    if (cmpNameInit(iCmp) .eq. cmpNameOk(iCmpOk)) then
                        lCmpOk = .true.
                        exit
                    end if
                end do
                if (.not. lCmpOk) call utmess('F', 'UTILITAI8_60')
            end do
!           Check elements types
            do iCell = 1, nbCellFilter
                if (listElemType(cellFilter(iCell)) .ne. pdtElemType) then
                    call utmess('F', 'UTILITAI8_60')
                end if
            end do
        end if

! ----- Loop on storing index
        do iStore = 1, nbStore
            if (lFromResult) then
                numeStore = listNumeStore(iStore)
                inst = listTimeStore(iStore)
                call getvtx(factorKeyword, 'NOM_CHAM', iocc=iocc, scal=fieldName, nbret=iret)
                if (iret .eq. 0) then
                    call utmess('F', 'POSTELEM_4')
                end if
                if (fieldName .eq. 'FORC_NODA' .or. fieldName .eq. 'REAC_NODA') then
                    call utmess('F', 'POSTELEM_5')
                end if
                call rsexch('F', resultIn, fieldName, numeStore, fieldInput, iret)
            else
                numeStore = nbStore
                fieldName = fieldFromUser
            end if
!
            call dismoi('TYPE_CHAMP', fieldInput, 'CHAMP', repk=fieldSupp, arret='C', ier=iret)
            call dismoi('NOM_GD', fieldInput, 'CHAMP', repk=physName, arret='C', ier=iret)
!
            if (physName(6:6) .eq. 'C') then
                exit
            end if
!
! --------- Prepare field
            convToNeut = ASTER_FALSE
            if (fieldSupp(1:2) .eq. 'EL') then
                field = fieldInput
            else
                convToNeut = ASTER_TRUE
                field = '&&PEEINT.FIELD'
                call convertFieldNodeToNeutElem(model, &
                                                ligrel, fieldInput, field, &
                                                nbCmpField, cmpNameNode, cmpNameNeut)
            end if
!
            if (.not. l_empty) then
                call dismoi('TYPE_CHAMP', field, 'CHAMP', repk=fieldSupp, arret='C', ier=iret)
                ASSERT(fieldSupp(1:2) .eq. 'EL')
            end if
!
            if (convToNeut .and. .not. l_empty) then
                do iCmp = 1, nbCmp
                    cmpNume = indik8(cmpNameNode, cmpNameInit(iCmp), 1, nbCmpField)
                    cmpName(iCmp) = cmpNameNeut(cmpNume)
                end do
            end if
            AS_DEALLOCATE(vk8=cmpNameNode)
            AS_DEALLOCATE(vk8=cmpNameNeut)
!
! --------- CALCUL ET STOCKAGE DES MOYENNES : MOT-CLE 'TOUT'
            call getvtx(factorKeyword, 'TOUT', iocc=iocc, nbval=0, nbret=iret)
            if (iret .ne. 0) then
                nbCellCompute = nbCellMesh
                AS_ALLOCATE(vi=cellCompute, size=nbCellCompute)
                do iCellCompute = 1, nbCellCompute
                    cellCompute(iCellCompute) = iCellCompute
                end do
                call peecal(fieldSupp, tableOut, fieldName, &
                            locaNameAll, locaNameAll, &
                            cellCompute, nbCellCompute, &
                            model, lFromResult, field, &
                            nbCmp, cmpName, cmpNameInit, &
                            numeStore, inst, iocc, &
                            ligrel, cespoi)
                AS_DEALLOCATE(vi=cellCompute)
            end if
!
! --------- CALCUL ET STOCKAGE DES MOYENNES : MOT-CLE 'GROUP_MA'
            call getvtx(factorKeyword, 'GROUP_MA', iocc=iocc, nbval=0, nbret=nbret)
            if (nbret .ne. 0) then
                nbGroup = -nbret
                AS_ALLOCATE(vk24=groupCell, size=nbGroup)
                call getvtx(factorKeyword, 'GROUP_MA', iocc=iocc, nbval=nbGroup, vect=groupCell)
                do iGroup = 1, nbGroup
                    groupName = groupCell(iGroup)
                    call jeexin(jexnom(mesh//'.GROUPEMA', groupName), iret)
                    if (iret .eq. 0 .and. .not. l_empty) then
                        call utmess('A', 'UTILITAI3_46', sk=groupName)
                        cycle
                    end if
                    nbCellCompute = 0
                    if (.not. l_empty) then
                        call jelira(jexnom(mesh//'.GROUPEMA', groupName), 'LONUTI', nbCellCompute)
                        if (nbCellCompute .eq. 0) then
                            call utmess('A', 'UTILITAI3_47', sk=groupName)
                            cycle
                        end if
                        call jeveuo(jexnom(mesh//'.GROUPEMA', groupName), 'L', vi=cellCompute)
                    end if
                    call peecal(fieldSupp, tableOut, fieldName, &
                                locaNameGroup, groupName, &
                                cellCompute, nbCellCompute, &
                                model, lFromResult, field, &
                                nbCmp, cmpName, cmpNameInit, &
                                numeStore, inst, iocc, &
                                ligrel, cespoi)
                end do
! --- UNION
                if (nbGroup > 1) then
                    call umalma(mesh, groupCell, nbGroup, cellCompute, nbCellCompute)
                    ASSERT(nbCellCompute > 0)
                    call peecal(fieldSupp, tableOut, fieldName, &
                                locaNameGroup, locaNameUnion, &
                                cellCompute, nbCellCompute, &
                                model, lFromResult, field, nbCmp, cmpName, &
                                cmpNameInit, numeStore, inst, iocc, ligrel, cespoi)
                    AS_DEALLOCATE(vi=cellCompute)
                end if
                AS_DEALLOCATE(vk24=groupCell)

            end if
        end do
!
        call jedetr(listCellUser)
        call jedetr(listCellFilter)
        if (l_newlig) call detrsd('LIGREL', ligrel)
        call detrsd('CHAM_ELEM_S', cespoi)
        call jedetr(cespoi//'.PDSM')
        AS_DEALLOCATE(vk16=variName)
        AS_DEALLOCATE(vk8=cmpName)
        AS_DEALLOCATE(vk8=cmpNameAll)
        AS_DEALLOCATE(vk8=cmpNameInit)
    end do
!
! - Clean
!
    call jedetr(numeStoreJv)
    call jedetr(timeStoreJv)
    if (lFromField) then
        call detrsd('CHAMP', fieldInput)
    end if
!
    call jedema()
!
end subroutine
