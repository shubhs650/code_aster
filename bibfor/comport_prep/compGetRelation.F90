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
subroutine compGetRelation(factorKeyword, iFactorKeyword, relaComp)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/deprecated_behavior.h"
#include "asterfort/getvtx.h"
!
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=16), intent(out) :: relaComp
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Get type of relation
!
! --------------------------------------------------------------------------------------------------
!
! In  iFactorKeyword   : factor keyword index
! Out relaComp         : behaviour (RELATION keyword)
!
! --------------------------------------------------------------------------------------------------
!
    relaComp = 'VIDE'
    call getvtx(factorKeyword, 'RELATION', iocc=iFactorKeyword, scal=relaComp)
    call deprecated_behavior(relaComp)
    if ((relaComp(1:4) .eq. 'META') .and. (relaComp .ne. 'META_LEMA_ANI')) then
        relaComp = 'KIT_META'
    end if
!
end subroutine
