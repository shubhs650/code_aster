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
! aslint: disable=W1504,W0104
!
subroutine lc0040(fami, kpg, ksp, ndim, imate, &
                  carcri, instam, instap, neps, &
                  epsm, deps, nsig, sigm, nvi, vim, &
                  option, sigp, vip, typmod, &
                  ndsde, dsidep, codret)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/lcdp_wrap.h"
#include "asterfort/Behaviour_type.h"

    integer(kind=8) :: imate, ndim, kpg, ksp, codret
    integer(kind=8) :: nvi, neps, nsig, ndsde
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(*), deps(*)
    real(kind=8) :: sigm(*), sigp(*)
    real(kind=8) :: vim(nvi), vip(nvi)
    real(kind=8) :: dsidep(*)
    character(len=16) :: option
    character(len=8) :: typmod(*)
    character(len=*) :: fami
! ----------------------------------------------------------------------
!  Loi de comportement DRUCK_PRAG_N_A
! ----------------------------------------------------------------------
    ASSERT(neps .eq. nint(sqrt(float(ndsde))))
    ASSERT(neps .eq. nsig)

    call lcdp_wrap(fami, kpg, ksp, ndim, imate, &
                   carcri, neps, epsm, &
                   deps, vim, option, sigp, vip, &
                   dsidep, codret)

end subroutine
