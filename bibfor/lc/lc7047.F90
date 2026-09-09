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

subroutine lc7047(BEHInteg, &
                  fami, kpg, ksp, ndim, imate, &
                  instam, instap, epsm, &
                  deps, nvi, vim, option, &
                  sigp, vip, typmod, &
                  dsidep, codret)
!
    use Behaviour_type
    implicit none
!
#include "asterfort/lcejdm.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate, nvi
    real(kind=8), intent(in) :: instam, instap
    real(kind=8), intent(in) :: epsm(6)
    real(kind=8), intent(in) :: deps(6)
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(6)
    real(kind=8), intent(out) :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! JOINT_MECA_ENDO
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    call lcejdm(BEHInteg, &
                fami, kpg, ksp, ndim, imate, &
                option, epsm, deps, sigp, &
                dsidep, vim, vip, typmod, &
                instam, instap)
end subroutine
