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
! aslint: disable=W1504,W0104,W1306,C1505
!
subroutine lc6057(BEHInteg, &
                  fami, kpg, ksp, ndim, imate, &
                  carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, typmod, ndsde, &
                  dsidep, codret)
!
    use Behaviour_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcesgv.h"
#include "asterfort/lcmfga.h"
#include "asterfort/lcmfma.h"
! --------------------------------------------------------------------------------------------------
    type(Behaviour_Integ) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
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
    real(kind=8) :: sigp(nsig)
    real(kind=8) :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out):: codret
! --------------------------------------------------------------------------------------------------
!   RELATION ENDO_FISS_EXP
! --------------------------------------------------------------------------------------------------
    aster_logical         :: lMatr, lSigm, lVari
    real(kind=8)          :: sig(nsig), dsde(nsig, neps), vi(nvi)
! --------------------------------------------------------------------------------------------------

    ASSERT(nsig .eq. neps)
    ASSERT(neps*nsig .eq. ndsde)

    sig = 0
    vi = 0
    dsde = 0

    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)

    if (lVari) vip = 0

    call lcesgv(fami, kpg, ksp, ndim, neps, typmod, option, imate, lcmfma, lcmfga, &
                epsm, deps, vim, nint(carcri(1)), carcri(3), sig, &
                vi, dsde, codret)

    if (codret .eq. 0) then
        if (lSigm) sigp = sig
        if (lVari) vip = vi
        if (lMatr) dsidep = dsde
    end if

end subroutine
