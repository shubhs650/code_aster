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
subroutine nmevr0(sddisc)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dieven.h"
#include "asterfort/utdidt.h"
#include "event_def.h"
!
    character(len=19), intent(in) :: sddisc
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - EVENEMENTS)
!
! REINITIALISATIONS DES EVENEMENTS
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc          : datastructure for time discretization
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lActivate
    integer(kind=8) :: iFail, nbFail, actionType
!
! --------------------------------------------------------------------------------------------------
!
    call utdidt('L', sddisc, 'LIST', 'NECHEC', vali_=nbFail)
    do iFail = 1, nbFail
        call dieven(sddisc, iFail, lActivate)
    end do
!
end subroutine
