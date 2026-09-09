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
    subroutine lc0054(fami, kpg, ksp, ndim, imate, &
                      compor, carcri, instam, instap, epsm, &
                      deps, sigm, nvi, vim, option, &
                      sigp, vip, typmod, &
                      dsidep, codret)
        integer(kind=8) :: imate, ndim, kpg, ksp, codret
        character(len=16), intent(in) :: compor(COMPOR_SIZE), option
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8) :: instam, instap
        integer(kind=8) :: nvi
        real(kind=8) :: epsm(6), deps(6)
        real(kind=8) :: sigm(6), sigp(6)
        real(kind=8) :: vim(nvi), vip(nvi)
        real(kind=8) :: dsidep(6, 6)
        character(len=8) :: typmod(*)
        character(len=*) :: fami
    end subroutine lc0054
end interface
