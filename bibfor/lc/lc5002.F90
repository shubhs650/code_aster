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
! aslint: disable=W1504,W0104,C1505,W1306

subroutine lc5002(BEHInteg, &
                  fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, typmod, ndsde, &
                  dsidep, codret)

    use Behaviour_type
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/nm1dis.h"
#include "asterfort/rcvalb.h"
#include "asterfort/verift.h"
#include "asterfort/Behaviour_type.h"
! --------------------------------------------------------------------------------------------------
    type(Behaviour_Integ)        :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
    character(len=16), intent(in) :: compor(COMPOR_SIZE), option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    real(kind=8)                 :: sigp(nsig)
    real(kind=8)                 :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out):: codret
! --------------------------------------------------------------------------------------------------
!   Relations VMIS_ISOT_LINE et VMIS_ISOT_TRAC
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: materPoin = " "
    aster_logical :: lMatr, lSigm, lVari
    integer(kind=8) :: ndimsi, propCode(1)
    real(kind=8) :: vi(nvi), propVale(1)
    real(kind=8) :: sig, deps_ther, deps_meca, em, ep, dsde
! --------------------------------------------------------------------------------------------------
    ndimsi = BEHInteg%behavPara%ndimsi
    ASSERT(ndimsi .eq. 1)

    sig = 0
    vi = 0
    dsde = 0

    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)

    if (lVari) vip = 0

    call rcvalb(BEHInteg%materPara%schemePara%fami, &
                BEHInteg%materPara%schemePara%kpg, &
                BEHInteg%materPara%schemePara%ksp, &
                '-', &
                BEHInteg%materPara%jvMaterCode, &
                materPoin, 'ELAS', &
                0, ' ', [0.d0], &
                1, 'E', propVale, propCode, 1)
    em = propVale(1)

    call rcvalb(BEHInteg%materPara%schemePara%fami, &
                BEHInteg%materPara%schemePara%kpg, &
                BEHInteg%materPara%schemePara%ksp, &
                '+', &
                BEHInteg%materPara%jvMaterCode, &
                materPoin, 'ELAS', &
                0, ' ', [0.d0], &
                1, 'E', propVale, propCode, 1)
    ep = propVale(1)
!
    call verift(BEHInteg%materPara%schemePara%fami, &
                BEHInteg%materPara%schemePara%kpg, &
                BEHInteg%materPara%schemePara%ksp, &
                'T', &
                BEHInteg%materPara%jvMaterCode, &
                epsth_=deps_ther)
    deps_meca = deps(1)-deps_ther

    call nm1dis(BEHInteg%materPara, &
                option, compor(RELA_NAME), materPoin, &
                em, ep, sigm(1), deps_meca, vim, &
                sig, vi, dsde)
    codret = 0

    if (codret .eq. 0) then
        if (lSigm) sigp(1:ndimsi) = sig
        if (lVari) vip(1:nvi) = vi
        if (lMatr) dsidep(1:ndimsi, 1:ndimsi) = dsde
    end if
end subroutine
