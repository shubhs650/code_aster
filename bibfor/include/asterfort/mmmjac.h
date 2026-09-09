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
#include "asterf_types.h"
#include "MeshTypes_type.h"
!
interface
    subroutine mmmjac(l_axis   , nb_node     , elem_dime,&
                      elem_code, elem_coor   ,&
                      ff       , dff         ,&
                      jacobi   , l_axis_warn_)
        aster_logical, intent(in) :: l_axis
        character(len=8), intent(in) :: elem_code
        integer(kind=8), intent(in) :: elem_dime, nb_node
        real(kind=8), intent(in) :: elem_coor(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: ff(MT_NNOMAX2D), dff(2, MT_NNOMAX2D)
        real(kind=8), intent(out) :: jacobi
        aster_logical, intent(out), optional :: l_axis_warn_
    end subroutine mmmjac
end interface
