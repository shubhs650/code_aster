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
#include "MeshTypes_type.h"
!
interface
    subroutine lcgeom_prep(elem_dime    , l_axis        , l_upda_jaco,&
                           nb_node_slav , elem_slav_code, elem_slav_init, elem_slav_coor,&
                           nb_node_mast , elem_mast_code, elem_mast_init, elem_mast_coor,&
                           gauss_coot_sl, gauss_coot_ma ,&
                           shape_func_sl, shape_func_ma ,&
                           jacobian_sl  , jacobian_ma   ,&
                           dist_vect_sl , dist_vect_ma)
        integer(kind=8), intent(in) :: elem_dime
        aster_logical, intent(in) :: l_axis, l_upda_jaco
        integer(kind=8), intent(in) :: nb_node_slav, nb_node_mast
        character(len=8), intent(in) :: elem_slav_code, elem_mast_code
        real(kind=8), intent(in) :: elem_mast_init(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: elem_slav_init(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: elem_mast_coor(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: elem_slav_coor(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: gauss_coot_sl(2), gauss_coot_ma(2)
        real(kind=8), intent(out) :: shape_func_sl(MT_NNOMAX2D), shape_func_ma(MT_NNOMAX2D)
        real(kind=8), intent(out) :: jacobian_sl, jacobian_ma
        real(kind=8), intent(out) :: dist_vect_sl(3), dist_vect_ma(3)
    end subroutine lcgeom_prep
end interface
