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
    subroutine vdesga(plateCara, plateOrie, &
                      kwgt, nb1, &
                      disp, btild, vectBaseKpg, &
                      hasTemp_, alpha_, tempKpg_, &
                      siefKpg_, epsiKpg_)
        use plate_type
        type(plateCara_Para), intent(in) :: plateCara
        type(plateOrie_Para), intent(in) :: plateOrie
        integer(kind=8), intent(in) :: kwgt, nb1
        real(kind=8), intent(in) :: disp(42), btild(5, 42), vectBaseKpg(3, 3)
        aster_logical, optional, intent(in) :: hasTemp_
        real(kind=8), optional, intent(in) :: alpha_, tempKpg_
        real(kind=8), optional, intent(out) :: siefKpg_(6, *), epsiKpg_(6, *)
    end subroutine vdesga
end interface
