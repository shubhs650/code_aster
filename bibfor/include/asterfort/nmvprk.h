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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine nmvprk(BEHinteg, &
                      option, typmod, ndim, &
                      compor, carcri, &
                      instam, instap, &
                      neps, epsdt, depst, sigd, &
                      nvi, vind, sigf, &
                      vinf, dsde, iret, multComp_)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: typmod(2)
        real(kind=8), intent(in) :: instam, instap
        integer(kind=8), intent(in) :: nvi
        character(len=16), optional, intent(in) :: multComp_
        integer(kind=8) :: neps
        integer(kind=8) :: ndim
        real(kind=8) :: epsdt(neps)
        real(kind=8) :: depst(neps)
        real(kind=8) :: sigd(6)
        real(kind=8) :: vind(*)
        real(kind=8) :: sigf(6)
        real(kind=8) :: vinf(*)
        real(kind=8) :: dsde(6, *)
        integer(kind=8) :: iret
    end subroutine nmvprk
end interface
