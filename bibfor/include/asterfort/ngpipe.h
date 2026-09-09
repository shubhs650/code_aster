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
    subroutine ngpipe(BEHInteg, typmod, compor, ndim, npg, neps, nddl, b, ni2ldc, &
                ddlm, ddld, ddl0, ddl1, lgpg, sigm_user, vim, &
                dtau, etamin, etamax, copilo)

        use Behaviour_type
#include "asterfort/Behaviour_type.h"

        type(Behaviour_Integ) :: BEHInteg
        character(len=8) :: typmod(2)
        character(len=16) ::compor(COMPOR_SIZE)
        integer(kind=8), intent(in):: ndim
        integer(kind=8), intent(in):: npg
        integer(kind=8), intent(in):: neps
        integer(kind=8), intent(in):: nddl
        real(kind=8), intent(in) :: b(neps, npg, nddl)
        real(kind=8), intent(in) :: ni2ldc(neps, npg)
        real(kind=8), intent(in) :: ddlm(nddl)
        real(kind=8), intent(in) :: ddld(nddl)
        real(kind=8), intent(in) :: ddl0(nddl)
        real(kind=8), intent(in) :: ddl1(nddl)
        integer(kind=8), intent(in) :: lgpg
        real(kind=8), intent(in) :: sigm_user(neps, npg)
        real(kind=8), intent(in) :: vim(lgpg, npg)
        real(kind=8), intent(in) :: dtau
        real(kind=8), intent(in) :: etamin
        real(kind=8), intent(in) :: etamax
        real(kind=8), intent(out) :: copilo(5, npg)
    end subroutine
end interface
