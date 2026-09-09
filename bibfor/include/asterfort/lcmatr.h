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
!
interface
    subroutine lcmatr(elem_dime   ,&
                      l_axis      , l_upda_jaco   , l_norm_smooth ,&
                      nb_lagr     , indi_lagc     , elga_fami     ,&
                      nb_node_slav, elem_slav_code, elem_slav_init, elem_slav_coor,&
                      nb_node_mast, elem_mast_code, elem_mast_init, elem_mast_coor,&
                      nb_poin_inte, poin_inte_sl  , poin_inte_ma  ,&
                      matr)
#include "MeshTypes_type.h"
        integer(kind=8), intent(in) :: elem_dime
        aster_logical, intent(in) :: l_axis, l_upda_jaco, l_norm_smooth
        integer(kind=8), intent(in) :: nb_lagr, indi_lagc(10)
        character(len=8), intent(in) :: elem_slav_code, elem_mast_code
        integer(kind=8), intent(in) :: nb_node_slav, nb_node_mast
        integer(kind=8), intent(in) :: nb_poin_inte
        real(kind=8), intent(in) :: poin_inte_sl(2,8), poin_inte_ma(2,8)
        real(kind=8), intent(in) :: elem_mast_init(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: elem_slav_init(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: elem_mast_coor(3, MT_NNOMAX2D)
        real(kind=8), intent(in) :: elem_slav_coor(3, MT_NNOMAX2D)
        character(len=8), intent(in) :: elga_fami
        real(kind=8), intent(inout) :: matr(55, 55)
    end subroutine lcmatr
end interface
