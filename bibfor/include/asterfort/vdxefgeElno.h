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
    subroutine vdxefgeElno(plateCara, plateOrie, &
                           nomte, nodeCoor, &
                           nbLayer, efgeElno, &
                           matevn)
        use plate_type
        type(plateCara_Para), intent(in) :: plateCara
        type(plateOrie_Para), intent(inout) :: plateOrie
        character(len=16), intent(in) :: nomte
        real(kind=8), intent(in) :: nodeCoor(3, 9)
        integer(kind=8), intent(in) :: nbLayer
        real(kind=8), intent(out) :: efgeElno(8, 9), matevn(2, 2, 10)
    end subroutine vdxefgeElno
end interface
