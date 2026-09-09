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
subroutine comp_meca_incr(lInitialState, relaComp, defoComp, typeComp)
!
    implicit none
!
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterc/lctest.h"
#include "asterf_types.h"
!
    aster_logical, intent(in) :: lInitialState
    character(len=16), intent(in) :: relaComp, defoComp
    character(len=16), intent(out) :: typeComp
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Select type of comportment (incremental or total)
!
! --------------------------------------------------------------------------------------------------
!
! In  lInitialState    : .true. if initial state is defined
! In  relaComp         : behaviour (RELATION keyword)
! In  defoComp         : model of strain (DEFORMATION keyword)
! Out typeComp         : type of behaviour (incremental or total)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=16) :: relaCompPY
!
! --------------------------------------------------------------------------------------------------
!
    call lccree(1, relaComp, relaCompPY)
    call lctest(relaCompPY, 'PROPRIETES', 'COMP_ELAS', iret)
    call lcdiscard(relaCompPY)
    if (iret .eq. 0) then
        typeComp = 'COMP_INCR'
    else
        typeComp = 'COMP_ELAS'
        if (lInitialState) then
            typeComp = 'COMP_INCR'
        end if
        if (defoComp .eq. 'PETIT_REAC') then
            typeComp = 'COMP_INCR'
        end if
        if (relaComp .eq. 'ELAS' .and. defoComp .eq. 'GROT_GDEP') then
            typeComp = 'COMP_INCR'
        end if
        if (relaComp .eq. 'ELAS' .and. defoComp .eq. 'GREEN_LAGRANGE') then
            typeComp = 'COMP_INCR'
        end if
    end if
!
end subroutine
