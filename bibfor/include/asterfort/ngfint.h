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
    subroutine ngfint(BEHinteg, &
                      option, typmod, ndim, nddl, neps, &
                      npg, w, b, compor, &
                      lgpg, carcri, instam, instap, &
                      ddlm, ddld, ni2ldc, sigmam, &
                      vim, sigmap, vip, fint, matsym, matuu, matns, &
                      lMatr, lVect, lSigm, codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        character(len=16) :: option
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        integer(kind=8) :: ndim, nddl, neps, npg, lgpg
        real(kind=8) :: w(neps, npg), ni2ldc(neps, npg), b(neps, npg, nddl)
        real(kind=8) :: instam, instap
        real(kind=8) :: ddlm(nddl), ddld(nddl)
        real(kind=8) :: sigmam(neps, npg), sigmap(neps, npg)
        real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg), fint(nddl)
        real(kind=8), intent(out):: matuu((nddl*(nddl+1))/2)
        real(kind=8), intent(out), target :: matns(nddl, nddl)
        aster_logical, intent(in) :: matsym, lMatr, lVect, lSigm
        integer(kind=8), intent(out) :: codret
    end subroutine ngfint
end interface
