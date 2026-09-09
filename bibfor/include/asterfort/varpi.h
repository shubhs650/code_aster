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
    subroutine varpi(ds_thm, &
                     p1, p1m, dp1, dp2, &
                     ep, surf, shut, &
                     phi0, dpi, sbjhm, &
                     wbjhm, epm, sbjh, wbjh)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        real(kind=8), intent(in) :: p1, p1m, dp2, dp1
        real(kind=8), intent(in) :: phi0
        real(kind=8), intent(in) :: ep, surf, shut, sbjh, wbjh
        real(kind=8), intent(out) :: dpi
        real(kind=8), intent(out) :: sbjhm, wbjhm, epm
    end subroutine varpi
end interface
