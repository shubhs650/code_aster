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
    subroutine merimo(jvBase, &
                      lXFEM, lMacrElem, &
                      model, caraElem, iterNewt, &
                      ds_constitutive, ds_material, ds_system, &
                      hval_incr, hval_algo, &
                      optionZ, ldccvg, sddynaZ)
        use NonLin_Datastructure_type
        character(len=1), intent(in) :: jvBase
        aster_logical, intent(in) :: lXFEM, lMacrElem
        character(len=24), intent(in) :: model, caraElem
        integer(kind=8), intent(in) :: iterNewt
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_System), intent(in) :: ds_system
        character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
        character(len=*), intent(in) :: optionZ
        integer(kind=8), intent(out) :: ldccvg
        character(len=*), optional, intent(in) :: sddynaZ
    end subroutine merimo
end interface
