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
!
interface
    subroutine lcejtu(BEHinteg, &
                      ndim, option, &
                      epsm, deps, sigm, sigp, &
                      dsidep, vim, vip, typmod, &
                      instam, instap)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: epsm(ndim), deps(ndim), sigm(6), vim(*)
        real(kind=8), intent(in) :: instam, instap
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option
        real(kind=8), intent(out) :: vip(*), sigp(6), dsidep(6, 6)
    end subroutine lcejtu
end interface
