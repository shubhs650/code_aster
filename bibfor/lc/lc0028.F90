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
subroutine lc0028(BEHInteg, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, instam, instap, epsm, &
                  deps, sigm, nvi, vim, option, &
                  sigp, vip, typmod, &
                  dsidep, codret)
!
    use Behaviour_type
    implicit none
!
#include "asterfort/nmvpir.h"
#include "asterfort/Behaviour_type.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    integer(kind=8) :: jvMaterCode, ndim, kpg, ksp, codret, nvi
    character(len=16), intent(in) :: compor(COMPOR_SIZE), option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: sigm(6), sigp(6)
    real(kind=8) :: vim(nvi), vip(nvi)
    real(kind=8) :: dsidep(6, 6)
    character(len=8) :: typmod(2)
    character(len=*) :: fami
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: relaComp
!
! --------------------------------------------------------------------------------------------------
!
    relaComp = compor(RELA_NAME)
    call nmvpir(BEHInteg, &
                fami, kpg, ksp, ndim, typmod, &
                jvMaterCode, relaComp, carcri, instam, instap, &
                deps, sigm, nvi, vim, option, &
                sigp, vip, dsidep, codret)
!
end subroutine
