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
    subroutine dilcar(option, compor, icontm, ivarim, ideplm, ideplp, &
                      igeom, imatuu, ivectu, icontp, &
                      ivarip, ichg, ichn, jcret, jvCarcri, iinstm, iinstp)
        character(len=16) :: option
        character(len=16), pointer :: compor(:)
        integer(kind=8) :: icontm
        integer(kind=8) :: ivarim
        integer(kind=8) :: ideplm
        integer(kind=8) :: ideplp
        integer(kind=8) :: igeom
        integer(kind=8) :: imatuu
        integer(kind=8) :: ivectu
        integer(kind=8) :: icontp
        integer(kind=8) :: ivarip
        integer(kind=8) :: ichg
        integer(kind=8) :: ichn
        integer(kind=8) :: jcret
        integer(kind=8) :: jvCarcri
        integer(kind=8) :: iinstm
        integer(kind=8) :: iinstp
    end subroutine dilcar
end interface
