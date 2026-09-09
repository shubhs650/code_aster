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
!
interface
#include "MeshTypes_type.h"
    subroutine apcoma(mesh, newgeo, elem_nume, elem_nbnode, elem_coor)
        character(len=8), intent(in) :: mesh
        character(len=19), intent(in) :: newgeo
        integer(kind=8), intent(in) :: elem_nume
        integer(kind=8), intent(in) :: elem_nbnode
        real(kind=8), intent(out) :: elem_coor(3,MT_NNOMAX2D)
    end subroutine apcoma
end interface
