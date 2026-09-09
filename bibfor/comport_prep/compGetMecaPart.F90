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
subroutine compGetMecaPart(relaComp, kitComp, relaMeca)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/comp_meca_l.h"
!
    character(len=16), intent(in) :: relaComp, kitComp(4)
    character(len=16), intent(out) :: relaMeca
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics) - Utility
!
! Get relation for mechanical part of behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  relaComp         : RELATION comportment
! In  kitComp          : KIT comportment
! Out relaMeca         : mecanical part of behaviour
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_kit_thm, l_kit_ddi
!
! --------------------------------------------------------------------------------------------------
!
    relaMeca = relaComp
    call comp_meca_l(relaComp, 'KIT_THM', l_kit_thm)
    call comp_meca_l(relaComp, 'KIT_DDI', l_kit_ddi)
    if (l_kit_thm) then
        relaMeca = kitComp(1)

    elseif (l_kit_ddi) then
        relaMeca = kitComp(1)

    end if
!
end subroutine
