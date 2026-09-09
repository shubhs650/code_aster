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
    subroutine affori(typ, nomt, cara, val, jad, jin,&
                      jdno, jdco, nutyma, ntseg, nutyel, ntyele, &
                      lseuil, nbseuil, alphayz)
        character(len=*) :: typ
        character(len=*) :: nomt
        character(len=*) :: cara
        real(kind=8) :: val(6)
        integer(kind=8) :: jad
        integer(kind=8) :: jin
        integer(kind=8) :: jdno
        integer(kind=8) :: jdco
        integer(kind=8) :: nutyma
        integer(kind=8) :: ntseg
        integer(kind=8) :: nutyel
        integer(kind=8) :: ntyele(*)
        real(kind=8), intent(in), optional :: lseuil
        integer(kind=8), intent(inout), optional :: nbseuil
        real(kind=8), intent(in), optional :: alphayz(2)
    end subroutine affori
end interface
