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
    subroutine nmvcpr_elem(modelZ, materFieldZ, materCodeZ, caraElemZ, &
                           poum, hval_incr, &
                           varcRefeZ, comporZ, &
                           vectElemZ)
        character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, materCodeZ
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: varcRefeZ, comporZ
        character(len=19), intent(in) :: hval_incr(*)
        character(len=*), intent(in) :: vectElemZ
    end subroutine nmvcpr_elem
end interface
