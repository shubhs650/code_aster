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
subroutine lc8029(BEHInteg, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, instam, instap, neps, &
                  epsm, deps, nsig, sigm, nvi, vim, &
                  option, sigp, vip, &
                  typmod, ndsde, dsidep, codret)
!
    use Behaviour_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmcomp.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam, instap
    integer(kind=8), intent(in) :: neps, nvi
    real(kind=8), intent(in) :: epsm(*)
    real(kind=8), intent(in) :: deps(*)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(*)
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(*)
    real(kind=8), intent(out) :: vip(nvi)
    character(len=8), intent(in) :: typmod(2)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(out) :: dsidep(*)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! KIT_DDI: BETON_UMLV / MAZARS
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    character(len=16) :: relaFlua, relaPlas
    character(len=16) :: comporFlua(COMPOR_SIZE)
    integer(kind=8) :: nume_flua, nvi_flua
!
! --------------------------------------------------------------------------------------------------
!
    comporFlua = 'VIDE'
    relaFlua = compor(CREEP_NAME)
    relaPlas = compor(PLAS_NAME)
    read (compor(CREEP_NVAR), '(I16)') nvi_flua
    read (compor(CREEP_NUME), '(I16)') nume_flua
    comporFlua(RELA_NAME) = relaFlua
    write (comporFlua(NVAR), '(I16)') nvi_flua
    comporFlua(DEFO) = compor(DEFO)
    write (comporFlua(NUME), '(I16)') nume_flua
    comporFlua(CREEP_NAME) = relaFlua
    comporFlua(PLAS_NAME) = relaPlas
    BEHInteg%behavPara%nvi = nvi_flua
    BEHInteg%behavPara%numlc = nume_flua
!
    call nmcomp(BEHInteg, &
                ndim, option, typmod, &
                instam, instap, &
                comporFlua, carcri, multComp, &
                neps, epsm, deps, &
                nsig, sigm, &
                vim, &
                sigp, vip, &
                ndsde, dsidep, &
                codret)
!
end subroutine
