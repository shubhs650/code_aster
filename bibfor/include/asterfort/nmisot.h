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
#include "asterf_types.h"
!
interface
    subroutine nmisot(fami, kpg, ksp, ndim, typmod, l_epsi_varc, &
                      imate, relaComp, carcri, deps, sigm, &
                      vim, option, sigp, vip, dsidep, &
                      iret)
        aster_logical, intent(in) :: l_epsi_varc
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        integer(kind=8) :: ndim
        character(len=8) :: typmod(*)
        integer(kind=8) :: imate
        character(len=16), intent(in) :: relaComp, option
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8) :: deps(6)
        real(kind=8) :: sigm(6)
        real(kind=8) :: vim(*)
        real(kind=8) :: sigp(6)
        real(kind=8) :: vip(*)
        real(kind=8) :: dsidep(6, 6)
        integer(kind=8) :: iret
    end subroutine nmisot
end interface
