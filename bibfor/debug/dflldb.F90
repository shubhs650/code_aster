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
subroutine dflldb(sdlist)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/assert.h"
#include "asterfort/dflld2.h"
#include "asterfort/dflld3.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: sdlist
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_LIST_INST
!
! Debug
!
! --------------------------------------------------------------------------------------------------
!
! In  sdlist           : name of DEFI_LIST_INST datastructure
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iFail, nbFail, nbInst, iAdap, nbAdap
    integer(kind=8) :: nbStepMaxi
    real(kind=8) :: dtmin, timeStepMini, timeStepMaxi
    character(len=24) :: sdlistEEvenrName
    real(kind=8), pointer :: sdlistEEvenr(:) => null()
    character(len=24) :: sdlistEEvenkName
    character(len=16), pointer :: sdlistEEvenk(:) => null()
    character(len=24) :: sdlistESubdrName
    real(kind=8), pointer :: sdlistESubdr(:) => null()
    character(len=24) :: sdlistInforName
    real(kind=8), pointer :: sdlistInfor(:) => null()
    character(len=24) :: sdlistAEvenrName
    real(kind=8), pointer :: sdlistAEvenr(:) => null()
    real(kind=8) :: valeRefe, peneMaxi, resiGlobMaxi
    character(len=16):: fieldName, cmpName, cmpCrit
    integer(kind=8) :: nb_incr_seuil, nb_iter_newt, crit_compi
    integer(kind=8) :: eventType, action_type
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to datastructures
    sdlistInforName = sdlist(1:8)//'.LIST.INFOR'
    call jeveuo(sdlistInforName, 'L', vr=sdlistInfor)
    sdlistEEvenrName = sdlist(1:8)//'.ECHE.EVENR'
    sdlistEEvenkName = sdlist(1:8)//'.ECHE.EVENK'
    sdlistESubdrName = sdlist(1:8)//'.ECHE.SUBDR'
    call jeveuo(sdlistEEvenrName, 'L', vr=sdlistEEvenr)
    call jeveuo(sdlistEEvenkName, 'L', vk16=sdlistEEvenk)
    call jeveuo(sdlistESubdrName, 'L', vr=sdlistESubdr)

! - Get main parameters
    timeStepMini = sdlistInfor(2)
    timeStepMaxi = sdlistInfor(3)
    nbStepMaxi = nint(sdlistInfor(4))
    dtmin = sdlistInfor(5)
    nbFail = nint(sdlistInfor(9))
    nbInst = nint(sdlistInfor(8))
    nbAdap = nint(sdlistInfor(10))

! - Time list management
    if (nint(sdlistInfor(1)) .eq. 1) then
        call utmess('I', 'DISCRETISATION3_1')
    else if (nint(sdlistInfor(1)) .eq. 2) then
        call utmess('I', 'DISCRETISATION3_2')
        call utmess('I', 'DISCRETISATION3_3', &
                    nr=2, valr=[timeStepMini, timeStepMaxi], &
                    si=nbStepMaxi)
    else
        ASSERT(.false.)
    end if
    call utmess('I', 'DISCRETISATION3_4', si=nbInst, sr=dtmin)

! - Failures
    if (nbFail .gt. 0) then
        call utmess('I', 'DISCRETISATION3_5', si=nbFail)
        do iFail = 1, nbFail
            valeRefe = sdlistEEvenr(SIZE_LEEVR*(iFail-1)+5)
            peneMaxi = sdlistEEvenr(SIZE_LEEVR*(iFail-1)+6)
            resiGlobMaxi = sdlistEEvenr(SIZE_LEEVR*(iFail-1)+7)
            fieldName = sdlistEEvenk(SIZE_LEEVK*(iFail-1)+1)
            cmpName = sdlistEEvenk(SIZE_LEEVK*(iFail-1)+2)
            cmpCrit = sdlistEEvenk(SIZE_LEEVK*(iFail-1)+3)
            eventType = nint(sdlistEEvenr(SIZE_LEEVR*(iFail-1)+1))
            if (eventType .eq. FAIL_EVT_ERROR) then
                call utmess('I', 'DISCRETISATION3_10', si=iFail)
            else if (eventType .eq. FAIL_EVT_INCR_QUANT) then
                call utmess('I', 'DISCRETISATION3_11', si=iFail)
                call utmess('I', 'DISCRETISATION3_21', &
                            nk=3, valk=[fieldName, cmpName, cmpCrit], &
                            sr=valeRefe)
            else if (eventType .eq. FAIL_EVT_INTERPENE) then
                call utmess('I', 'DISCRETISATION3_13', si=iFail)
                call utmess('I', 'DISCRETISATION3_22', sr=peneMaxi)
            else if (eventType .eq. FAIL_EVT_DIVE_RESI) then
                call utmess('I', 'DISCRETISATION3_14', si=iFail)
            else if (eventType .eq. FAIL_EVT_INSTABILITY) then
                call utmess('I', 'DISCRETISATION3_15', si=iFail)
            else if (eventType .eq. FAIL_EVT_RESI_MAXI) then
                call utmess('I', 'DISCRETISATION3_16', si=iFail)
                call utmess('I', 'DISCRETISATION3_23', sr=resiGlobMaxi)
            else
                ASSERT(.false.)
            end if

! --------- Action
            action_type = nint(sdlistEEvenr(SIZE_LEEVR*(iFail-1)+2))
            if (action_type .eq. FAIL_ACT_STOP) then
                call utmess('I', 'DISCRETISATION3_30')
            else if (action_type .eq. FAIL_ACT_CUT) then
                call utmess('I', 'DISCRETISATION3_31')
                call dflld2(sdlist, iFail)
            else if (action_type .eq. FAIL_ACT_CONTINUE) then
                call utmess('I', 'DISCRETISATION3_36')
            else if (action_type .eq. FAIL_ACT_ARCHI) then
                call utmess('I', 'DISCRETISATION3_37')
            else
                ASSERT(.false.)
            end if
        end do
    end if

! - Adaptation
    if (nbAdap .gt. 0) then
        sdlistAEvenrName = sdlist(1:8)//'.ADAP.EVENR'
        call jeveuo(sdlistAEvenrName, 'L', vr=sdlistAEvenr)
        call utmess('I', 'DISCRETISATION3_6', si=nbAdap)
        do iAdap = 1, nbAdap
            nb_incr_seuil = nint(sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+2))
            nb_iter_newt = nint(sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+5))
            crit_compi = nint(sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+4))
            eventType = nint(sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+1))
            if (eventType .eq. ADAP_EVT_NONE) then
                call utmess('I', 'DISCRETISATION3_50', si=iAdap)
            else if (eventType .eq. ADAP_EVT_ALLSTEPS) then
                call utmess('I', 'DISCRETISATION3_51', si=iAdap)
                call dflld3(sdlist, iAdap)
            else if (eventType .eq. ADAP_EVT_TRIGGER) then
                call utmess('I', 'DISCRETISATION3_52', si=iAdap)
                if (crit_compi .eq. 1) then
                    call utmess('I', 'DISCRETISATION3_64', ni=2, &
                                vali=[nb_incr_seuil, nb_iter_newt])
                elseif (crit_compi .eq. 2) then
                    call utmess('I', 'DISCRETISATION3_66', ni=2, &
                                vali=[nb_incr_seuil, nb_iter_newt])
                elseif (crit_compi .eq. 3) then
                    call utmess('I', 'DISCRETISATION3_63', ni=2, &
                                vali=[nb_incr_seuil, nb_iter_newt])
                elseif (crit_compi .eq. 4) then
                    call utmess('I', 'DISCRETISATION3_65', ni=2, &
                                vali=[nb_incr_seuil, nb_iter_newt])
                end if
                call dflld3(sdlist, iAdap)
            end if
        end do
    end if
!
    call jedema()
end subroutine
