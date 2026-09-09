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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine lc0009(fami, kpg, ksp, ndim, imate, &
                      instam, instap, epsm, &
                      deps, sigm, vim, option, &
                      sigp, vip, typmod, &
                      nvi, dsidep, codret)
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        integer(kind=8) :: ndim
        integer(kind=8) :: imate
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: epsm(*)
        real(kind=8) :: deps(*)
        real(kind=8) :: sigm(*)
        real(kind=8) :: vim(nvi)
        character(len=16) :: option
        real(kind=8) :: sigp(*)
        real(kind=8) :: vip(nvi)
        character(len=8) :: typmod(*)
        integer(kind=8) :: nvi
        real(kind=8) :: dsidep(*)
        integer(kind=8) :: codret
    end subroutine lc0009
end interface
