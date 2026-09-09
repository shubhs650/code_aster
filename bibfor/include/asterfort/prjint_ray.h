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
    subroutine prjint_ray(proj_tole, dist_ratio       , elem_dime     ,&
                      elem_mast_nbnode, elem_mast_coor, elem_mast_code,&
                      elem_slav_nbnode, elem_slav_coor, elem_slav_code,&
                      poin_inte_ma, poin_inte_es, inte_weight   , nb_poin_inte  ,&
                      inte_neigh_     , ierror_)
#include "MeshTypes_type.h"
        real(kind=8), intent(in) :: proj_tole, dist_ratio
        integer(kind=8), intent(in) :: elem_dime
        real(kind=8), intent(in) :: elem_slav_coor(3,MT_NNOMAX2D)
        integer(kind=8), intent(in) :: elem_slav_nbnode
        character(len=8), intent(in) :: elem_slav_code
        real(kind=8), intent(in) :: elem_mast_coor(3,MT_NNOMAX2D)
        integer(kind=8), intent(in) :: elem_mast_nbnode
        character(len=8), intent(in) :: elem_mast_code
        real(kind=8), intent(out) :: poin_inte_ma(2,8)
        real(kind=8), intent(out) :: poin_inte_es(2,8)
        real(kind=8), intent(out) :: inte_weight
        integer(kind=8), intent(out) :: nb_poin_inte
        integer(kind=8), optional, intent(inout) :: inte_neigh_(4)
        integer(kind=8), optional, intent(inout) :: ierror_
    end subroutine prjint_ray
end interface
