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
! aslint: disable=W1504,C1505,W1306,W0104
!
subroutine lc6076(BEHInteg, &
                  fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, typmod, ndsde, &
                  dsidep, codret)
!
    use Behaviour_type
    use vmis_isot_nl_module, only: CONSTITUTIVE_LAW, Init, InitViscoPlasticity, Integrate, &
                                   InitGradVari
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcgrad.h"
! --------------------------------------------------------------------------------------------------
    type(Behaviour_Integ)        :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam, instap
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8)                 :: sigp(nsig)
    real(kind=8)                 :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(out) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                        merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out) :: codret
! --------------------------------------------------------------------------------------------------
!   RELATION VMIS_ISOT_NL ET VISC_ISOT_NL + GRAD_VARI
! --------------------------------------------------------------------------------------------------
    aster_logical :: lMatr, lSigm, lVari
    integer(kind=8) :: ndimsi
    real(kind=8) :: sig(2*ndim), vi(nvi), vinl
    real(kind=8) :: deps_sig(2*ndim, 2*ndim), deps_vi(2*ndim), dphi_sig(2*ndim), dphi_vi
    real(kind=8) :: apg, lag, grad(ndim), eps_gene(neps), eps_meca(2*ndim)
    type(CONSTITUTIVE_LAW):: cl
! --------------------------------------------------------------------------------------------------
!
    ASSERT(neps*nsig .eq. ndsde)
    ASSERT(neps .eq. nsig)
    ASSERT(neps .ge. 2*ndim)
    ASSERT(neps .ge. 3*ndim+2)
!
    ndimsi = 2*ndim
    sig = 0
    vi = 0
    deps_sig = 0
    deps_vi = 0
    dphi_sig = 0
    dphi_vi = 0
    eps_gene = epsm(1:neps)+deps(1:neps)
    eps_meca = eps_gene(1:ndimsi)
    apg = eps_gene(ndimsi+1)
    lag = eps_gene(ndimsi+2)
    grad = eps_gene(ndimsi+3:ndimsi+2+ndim)

    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)

    if (lVari) vip = 0

    cl = Init(ndimsi, option, fami, kpg, ksp, imate, &
              nint(carcri(ITER_INTE_MAXI)), carcri(RESI_INTE))

    call InitGradVari(cl, fami, kpg, ksp, imate, lag, apg)

    if (compor(RELA_NAME) (1:4) .eq. 'VISC') &
        call InitViscoPlasticity(cl, fami, kpg, ksp, imate, instap-instam)

    call Integrate(cl, eps_meca, vim(1:nvi), sig, &
                   vi, deps_sig, vinl, dphi_sig, deps_vi, dphi_vi)

    codret = cl%exception
    if (codret .eq. 0) then
        if (lVari) vip(1:nvi) = vi
        call lcgrad(lSigm, lMatr, sig, apg, lag, grad, vinl, &
                    cl%mat%r, cl%mat%c, deps_sig, dphi_sig, deps_vi, dphi_vi, sigp, dsidep)
    end if
!
end subroutine
