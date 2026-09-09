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
    subroutine vectan(nb1, nb2, &
                      nodeCoor, desr, &
                      vectNorm, vectTang)
        integer(kind=8), intent(in) :: nb1, nb2
        real(kind=8), intent(in) :: nodeCoor(3, *)
        real(kind=8), intent(inout) :: desr(*)
        real(kind=8), intent(out) :: vectNorm(9, 3), vectTang(9, 2, 3)
    end subroutine vectan
end interface
