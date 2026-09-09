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

subroutine pinorm(a, b, c, dtau, copilo)

    implicit none

#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/zerop2.h"

    real(kind=8), intent(in):: a(:)
    real(kind=8), intent(in):: b(:)
    real(kind=8), intent(in):: c
    real(kind=8), intent(in):: dtau
    real(kind=8), intent(out):: copilo(:)
! --------------------------------------------------------------------------------------------------
!  Solve path-following equation || eta*a + b || = c
! --------------------------------------------------------------------------------------------------
!  dtau: path-following increment (to adjust to former copilo interface)
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: nsol
    real(kind=8):: sol(2), a2, p0, p1
! --------------------------------------------------------------------------------------------------
    ASSERT(size(copilo) .eq. 5)

    ! Too small (or negative) path-following increment
    if (c .le. 0) then
        nsol = 0

        ! quasi-constant path-following function (vector a almost equal to zero)
    else if (norm2(a) .le. c/r8gaem()) then
        nsol = merge(-1, 0, norm2(b) .le. c)

        ! Usual quadratic function
    else

        ! Solution of P2 polynom in decreasing order
        a2 = dot_product(a, a)
        p0 = (dot_product(b, b)-c**2)/a2
        p1 = 2*dot_product(a, b)/a2
        call zerop2(p1, p0, sol, nsol)

        ! Interval reduced to a single point -> no regular solution
        if (nsol .eq. 1) nsol = 0
    end if

    ! Path-following coefficients (in old format copilo)
    copilo = r8vide()
    select case (nsol)
    case (0)
        copilo(5) = 0.d0

    case (2)
        copilo(1) = dtau-sol(1)
        copilo(2) = 1.d0
        copilo(3) = dtau+sol(2)
        copilo(4) = -1.d0

    case default
        ASSERT(ASTER_FALSE)
    end select

end subroutine
