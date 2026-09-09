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
    subroutine relax_acier_cable(fami, kpg, ksp, imate, sigm, epsm, deps, vim, &
                                 sigp, vip, dsde)
        character(len=*)  :: fami
        integer(kind=8)           :: kpg
        integer(kind=8)           :: ksp
        integer(kind=8)           :: imate
        real(kind=8)      :: sigm
        real(kind=8)      :: epsm
        real(kind=8)      :: deps
        real(kind=8)      :: vim(*)
        real(kind=8)      :: sigp
        real(kind=8)      :: vip(*)
        real(kind=8)      :: dsde
    end subroutine relax_acier_cable
end interface
