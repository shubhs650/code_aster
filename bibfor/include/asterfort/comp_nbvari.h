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
    subroutine comp_nbvari(relaComp, defoComp, typeCpla, kitComp, &
                           postIter, multComp, reguVisc, postIncr, &
                           solvBehavType, adrsMGIS, &
                           nbVariUMAT, &
                           nbVari, numeLaw, nbVariKit, numeLawKit)
        character(len=16), intent(in) :: relaComp, defoComp, typeCpla, kitComp(4)
        character(len=16), intent(in) :: postIter, multComp, reguVisc, postIncr
        integer(kind=8), intent(in) :: solvBehavType
        character(len=16), intent(in) :: adrsMGIS
        integer(kind=8), intent(in) :: nbVariUMAT
        integer(kind=8), intent(out) :: nbVari, numeLaw, nbVariKit(4), numeLawKit(4)
    end subroutine comp_nbvari
end interface
