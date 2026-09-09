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
!
#include "asterf_types.h"
!
interface
    subroutine dinon3(neq, ul, dul, utl, nno,&
                      nbcomp, varimo, raide, nbpar, param,&
                      okdire, varipl, dforcz)
        integer(kind=8), intent(in) :: nbpar
        integer(kind=8), intent(in) :: nbcomp
        integer(kind=8), intent(in) :: neq
        real(kind=8), intent(in) :: ul(neq)
        real(kind=8), intent(in) :: dul(neq)
        real(kind=8), intent(in) :: utl(neq)
        integer(kind=8), intent(in) :: nno
        real(kind=8), intent(in) :: varimo(nbcomp*3)
        real(kind=8), intent(inout) :: raide(nbcomp)
        real(kind=8), intent(in) :: param(6, nbpar)
        aster_logical, intent(in) :: okdire(6)
        real(kind=8), intent(out) :: varipl(nbcomp*3)
        real(kind=8), intent(out), optional :: dforcz(6)
    end subroutine dinon3
end interface
