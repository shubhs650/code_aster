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
interface
    subroutine lcgeog(elem_dime     , i_reso_geom   ,&
                      nb_lagr       , indi_lagc     ,&
                      nb_node_slav  , nb_node_mast  ,&
                      elem_mast_init, elem_slav_init,&
                      elem_mast_coor, elem_slav_coor)
#include "MeshTypes_type.h"
        integer(kind=8), intent(in) :: elem_dime, i_reso_geom
        integer(kind=8), intent(in) :: nb_lagr, indi_lagc(10)
        integer(kind=8), intent(in) :: nb_node_slav, nb_node_mast
        real(kind=8), intent(in) :: elem_slav_init(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: elem_mast_init(3, MT_NNOMAX2D)
        real(kind=8), intent(inout) :: elem_slav_coor(3, MT_NNOMAX2D)
        real(kind=8), intent(inout) :: elem_mast_coor(3, MT_NNOMAX2D)
    end subroutine lcgeog
end interface
