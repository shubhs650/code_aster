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
    subroutine ptinma(elem_nbnode, elem_dime , elem_code, elem_coor, pair_tole,&
                      poin_coorx , poin_coory, test, cor_inte_ori)
        integer(kind=8), intent(in) :: elem_nbnode
        integer(kind=8), intent(in) :: elem_dime
        character(len=8), intent(in) :: elem_code
        real(kind=8), intent(in) :: elem_coor(2,elem_nbnode)
        real(kind=8), intent(in) :: pair_tole
        real(kind=8), intent(in) :: poin_coorx
        real(kind=8), intent(in) :: poin_coory
        integer(kind=8), intent(out) :: test
        real(kind=8),optional ,intent(out) :: cor_inte_ori(2)
    end subroutine ptinma
end interface
