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
subroutine nmfinp(sddisc, numeInst, lastTimeStep)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/didern.h"
#include "asterfort/diinst.h"
#include "asterfort/nmjalo.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "event_def.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeInst
    aster_logical, intent(out) :: lastTimeStep
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE
!
! Detect last time step
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  numeInst         : index of current time step
! Out lastTimeStep     : flag for last time step
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: prec, jalon, timeCurr
    character(len=16) :: metlis
    integer(kind=8) :: nbStepMaxi
    aster_logical :: lOutOFTime
!
! --------------------------------------------------------------------------------------------------
!
    lastTimeStep = ASTER_FALSE

! - Current time step
    timeCurr = diinst(sddisc, numeInst)

! - PRECISION SUR LES INSTANTS
! - (LIEE A CELLE DE VAL_MIN DE PAS_MINI DANS DEFI_LIST_INST.CAPY)
    prec = 1.d-12

! - METHODE DE GESTION DE LA LISTE D'INSTANTS
    call utdidt('L', sddisc, 'LIST', 'METHODE', valk_=metlis)

! - Detect last time step by "manual" time stepping
    if (didern(sddisc, numeInst)) then
        lastTimeStep = ASTER_TRUE
    end if

! - Detect last time step by NB_PAS_MAXI
    call utdidt('L', sddisc, 'LIST', 'NB_PAS_MAXI', vali_=nbStepMaxi)
    if (numeInst .ne. 0) then
        if (numeInst .eq. nbStepMaxi) then
            call utmess('I', 'ADAPTATION_13')
            lastTimeStep = ASTER_TRUE
        end if
    end if

! -  Detect last time step by automatic time stepping
    if (metlis .eq. 'AUTO') then
        call nmjalo(sddisc, timeCurr, prec, jalon, lOutOFTime)
        if (lOutOFTime) then
            lastTimeStep = ASTER_TRUE
        end if
    end if
!
end subroutine
