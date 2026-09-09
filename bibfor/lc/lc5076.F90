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

subroutine lc5076(BEHInteg, &
                  fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, typmod, ndsde, &
                  dsidep, codret)

    use Behaviour_type
    use vmis_isot_nl_module, only: CONSTITUTIVE_LAW, Init, InitViscoPlasticity, Integrate
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"
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
!   RELATIONS VMIS_ISOT_NL ET VISC_ISOT_NL
! --------------------------------------------------------------------------------------------------
    aster_logical :: lMatr, lSigm, lVari
    integer(kind=8) :: ndimsi
    real(kind=8) :: vi(nvi)
    real(kind=8) :: sig(BEHInteg%behavPara%ndimsi), eps(BEHInteg%behavPara%ndimsi)
    real(kind=8) :: dsde(BEHInteg%behavPara%ndimsi, BEHInteg%behavPara%ndimsi)
    type(CONSTITUTIVE_LAW):: cl
! --------------------------------------------------------------------------------------------------
    ndimsi = BEHInteg%behavPara%ndimsi
    sig = 0
    vi = 0
    dsde = 0
    eps = epsm(1:ndimsi)+deps(1:ndimsi)

    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)

    if (lVari) vip = 0

    cl = Init(ndimsi, option, fami, kpg, ksp, imate, &
              nint(carcri(ITER_INTE_MAXI)), carcri(RESI_INTE))

    if (compor(RELA_NAME) (1:4) .eq. 'VISC') &
        call InitViscoPlasticity(cl, fami, kpg, ksp, imate, instap-instam)

    call Integrate(cl, eps, vim(1:nvi), sig, vi, dsde)

    codret = cl%exception

    if (codret .eq. 0) then
        if (lSigm) sigp(1:ndimsi) = sig
        if (lVari) vip(1:nvi) = vi
        if (lMatr) dsidep(1:ndimsi, 1:ndimsi) = dsde
    end if
end subroutine
