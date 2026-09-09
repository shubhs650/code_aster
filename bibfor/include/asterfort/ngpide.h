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

interface
    subroutine ngpide(compor, npg, neps, nddl, b, &
                ddlm, ddld, ddl0, ddl1, dtau, copilo, neps_meca)

#include "asterfort/Behaviour_type.h"

        character(len=16) ::compor(COMPOR_SIZE)
        integer(kind=8), intent(in):: npg
        integer(kind=8), intent(in):: neps
        integer(kind=8), intent(in):: nddl
        real(kind=8), intent(in) :: b(neps, npg, nddl)
        real(kind=8), intent(in) :: ddlm(nddl)
        real(kind=8), intent(in) :: ddld(nddl)
        real(kind=8), intent(in) :: ddl0(nddl)
        real(kind=8), intent(in) :: ddl1(nddl)
        real(kind=8), intent(in) :: dtau
        real(kind=8), intent(out) :: copilo(5, npg)
        integer(kind=8), intent(in), optional:: neps_meca
    end subroutine
end interface
