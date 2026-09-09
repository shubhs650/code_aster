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
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine nmcomp(BEHinteg, &
                      ndim, option, typmod, &
                      instam, instap, &
                      compor, carcri, multComp, &
                      neps, epsm_inp, deps_inp, &
                      nsig, sigm, &
                      vim, &
                      sigp, vip, &
                      ndsde, dsidep, codret, &
                      l_epsi_varc_)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        integer(kind=8), intent(in) :: ndim
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: typmod(2)
        real(kind=8), intent(in) :: instam, instap
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        character(len=16), intent(in) :: multComp
        integer(kind=8), intent(in) :: neps
        real(kind=8), intent(in) :: epsm_inp(neps), deps_inp(neps)
        integer(kind=8), intent(in) :: nsig
        real(kind=8), intent(in) :: sigm(nsig)
        real(kind=8), intent(in) :: vim(*)

        real(kind=8), intent(inout) :: sigp(nsig), vip(*)
        integer(kind=8), intent(in) :: ndsde
        real(kind=8), intent(inout) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                              merge(neps, 6, nsig*neps .eq. ndsde))
        integer(kind=8), intent(inout) :: codret
        aster_logical, optional, intent(in) :: l_epsi_varc_
    end subroutine nmcomp
end interface
