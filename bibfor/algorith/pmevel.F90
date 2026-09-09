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
subroutine pmevel(sddisc, tablType, tablIncr, &
                  lerror, conver)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/eneven.h"
#include "asterfort/getFailEvent.h"
#include "asterfort/pmevdg.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "event_def.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: tablType
    character(len=24), intent(in) :: tablIncr
    aster_logical, intent(in) :: lerror, conver
!
! --------------------------------------------------------------------------------------------------
!
! SIMU_POINT_MAT
!
! Detect first event
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbEvent, iEvent, iEventActi
    integer(kind=8) :: eventType
!
! --------------------------------------------------------------------------------------------------
!
    iEventActi = 0

! - Number of events
    call utdidt('L', sddisc, 'LIST', 'NECHEC', vali_=nbEvent)

    do iEvent = 1, nbEvent
! ----- Get event type
        call getFailEvent(sddisc, iEvent, eventType)

! ----- Set event to inactive
        call eneven(sddisc, iEvent, ASTER_FALSE)

! ----- Detect event
        if (eventType .eq. FAIL_EVT_ERROR) then
            if (lerror) then
                iEventActi = iEvent
                goto 99
            end if

        else if (eventType .eq. FAIL_EVT_INCR_QUANT) then
            if (conver) then
                call pmevdg(sddisc, tablType, tablIncr, iEvent, iEventActi)
                if (iEventActi .ne. 0) then
                    goto 99
                end if
            end if

        else
            call utmess('F', 'COMPOR2_9')

        end if
    end do
!
99  continue

! - Set event to active
    if (iEventActi .ne. 0) then
        call eneven(sddisc, iEventActi, ASTER_TRUE)
    end if
!
end subroutine
