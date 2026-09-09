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
    subroutine nmasym(materPara, option, &
                      xlong0, a, dlong0, &
                      effnom, vim, effnop, vip, klv, &
                      fono)
        use MaterialPara_type
        integer(kind=8), parameter :: neq = 6, nbt = 21, nvar = 4
        type(Material_Para), intent(in) :: materPara
        character(len=*) :: option
        real(kind=8) :: xlong0, a, syc, syt, etc, ett, cr
        real(kind=8) :: e, dlong0
        real(kind=8) :: effnom, vim(nvar)
        real(kind=8) :: effnop, vip(nvar), fono(neq), klv(nbt)
    end subroutine nmasym
end interface
