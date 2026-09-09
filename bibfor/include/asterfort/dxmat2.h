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
    subroutine dxmat2(plateCara, plateOrie, &
                      iLayer, npg, &
                      ordi, epi, &
                      epais, dm)
        use plate_type
        type(plateCara_Para), intent(in) :: plateCara
        type(plateOrie_Para), intent(in) :: plateOrie
        integer(kind=8), intent(in) :: iLayer, npg
        real(kind=8), intent(out) :: ordi, epi
        real(kind=8), intent(out) :: epais, dm(3, 3)
    end subroutine dxmat2
end interface
