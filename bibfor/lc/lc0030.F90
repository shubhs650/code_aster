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
subroutine lc0030(BEHInteg, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, instam, instap, epsm, &
                  deps, sigm, nvi, vim, option, &
                  sigp, vip, &
                  typmod, dsidep, codret)
!
    use Behaviour_type
    implicit none
!
#include "asterfort/plasti.h"
#include "asterfort/Behaviour_type.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    real(kind=8), intent(in) :: epsm(*)
    real(kind=8), intent(in) :: deps(*)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: sigm(6)
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(6)
    real(kind=8), intent(out) :: vip(nvi)
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! beton_burger, irrad3m, rouss_pr, rouss_visc
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHInteg       : parameters for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
    call plasti(BEHInteg, &
                option, typmod, &
                fami, kpg, ksp, jvMaterCode, &
                compor, carcri, instam, instap, &
                epsm, deps, &
                sigm, &
                nvi, vim, &
                sigp, vip, &
                dsidep, codret)
!
end subroutine
