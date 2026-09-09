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
    subroutine lcumfp(fami, kpg, ksp, ndim, typmod, &
                      imate, tinstm, tinstp, epsm, &
                      deps, sigm, vim, option, relaPlas, &
                      sigp, vip, dsidep)
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: imate
        integer(kind=8), intent(in) :: kpg
        integer(kind=8), intent(in) :: ksp
        character(len=8), intent(in) :: typmod(*)
        character(len=16), intent(in) :: relaPlas
        character(len=16), intent(in) :: option
        character(len=*), intent(in) :: fami
        real(kind=8) :: tinstm, tinstp
        real(kind=8) :: epsm(*), deps(*), sigm(*), sigp(*), vim(*), vip(*)
        real(kind=8) :: dsidep(6, 6)
    end subroutine lcumfp
end interface
