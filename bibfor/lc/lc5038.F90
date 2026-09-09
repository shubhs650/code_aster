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
! aslint: disable=C1505

subroutine lc5038(BEHInteg, neps, nsig, nvi, option, sigp, vip, ndsde, dsidep, codret)
    use Behaviour_type
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
! --------------------------------------------------------------------------------------------------
    type(Behaviour_Integ) :: BEHinteg
    character(len=16), intent(in) :: option
    integer(kind=8), intent(in) :: neps, nsig
    integer(kind=8), intent(in) :: nvi
    real(kind=8)                 :: sigp(nsig)
    real(kind=8)                 :: vip(nvi)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out):: codret
! --------------------------------------------------------------------------------------------------
!   RELATION SANS
! --------------------------------------------------------------------------------------------------
    aster_logical :: lMatr, lSigm, lVari
    integer(kind=8) :: ndimsi
! --------------------------------------------------------------------------------------------------
    ndimsi = BEHInteg%behavPara%ndimsi

    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)

    codret = 0
    if (lSigm) sigp(1:ndimsi) = 0
    if (lVari) vip(1:nvi) = 0
    if (lMatr) dsidep(1:ndimsi, 1:ndimsi) = 0

end subroutine
