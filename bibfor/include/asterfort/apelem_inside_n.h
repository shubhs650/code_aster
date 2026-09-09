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
#include "asterf_types.h"
!
interface
    subroutine apelem_inside_n(pair_tole   , elem_dime, elem_code, &
                               nb_poin_coor, poin_coor,&
                               nb_poin_inte, poin_inte, inte_neigh,&
                               poin_coor_ori, poin_inte_ori)
        real(kind=8), intent(in) :: pair_tole
        integer(kind=8), intent(in) :: elem_dime
        character(len=8), intent(in) :: elem_code
        integer(kind=8), intent(in) :: nb_poin_coor
        real(kind=8), intent(in) :: poin_coor(2,4)
        real(kind=8), intent(in) :: poin_coor_ori(2,4)
        integer(kind=8), intent(inout) :: nb_poin_inte
        real(kind=8), intent(inout) :: poin_inte(2,16)
        real(kind=8), intent(inout) :: poin_inte_ori(2,16)
        integer(kind=8), intent(inout) :: inte_neigh(4)
    end subroutine apelem_inside_n
end interface
