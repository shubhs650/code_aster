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
    subroutine interNodesEdge(proj_tole       , elem_dime     , &
                              elem_mast_code, elem_slave_code, &
                              proj_coor       , nb_node_proj, &
                              nb_poin_inte,  poin_inte, inte_neigh,&
                              poin_inte_ori)
        real(kind=8), intent(in) :: proj_tole
        integer(kind=8), intent(in) :: elem_dime
        character(len=8), intent(in) :: elem_mast_code, elem_slave_code
        real(kind=8), intent(in) :: proj_coor(2,9)
        integer(kind=8), intent(in) :: nb_node_proj
        integer(kind=8), intent(inout) :: inte_neigh(4), nb_poin_inte
        real(kind=8), intent(out) :: poin_inte(2,16)
        real(kind=8), intent(out) :: poin_inte_ori(2,16)
    end subroutine interNodesEdge
end interface
