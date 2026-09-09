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
#include "asterf_types.h"
!
interface
    subroutine pmconv(resi, resiInit, resiEval, &
                      ds_conv, &
                      timeCurr, iterNewt, &
                      coefAdim, sigmCurr, &
                      conver, lIterNewtMaxi)
        use NonLin_Datastructure_type
        real(kind=8), intent(in) :: resi(12), resiInit(12)
        real(kind=8), intent(inout) :: resiEval(12)
        type(NL_DS_Conv), intent(in) :: ds_conv
        real(kind=8), intent(in) :: timeCurr
        integer(kind=8), intent(in) :: iterNewt
        real(kind=8), intent(in) :: coefAdim, sigmCurr(6)
        aster_logical, intent(out) :: conver, lIterNewtMaxi
    end subroutine pmconv
end interface
