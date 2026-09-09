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

subroutine longeleori(jdno, jdco, ino1, ino2, longseg, x3)
!
!
! --------------------------------------------------------------------------------------------------
!
!                           CALCULE LA LONGEUR DES SEGMENTS
!
!   OUT
!       longseg   : longueur du segment
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: jdno, jdco, ino1, ino2
    real(kind=8) :: longseg, x3(3)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: no1, no2, ii
    real(kind=8) :: x1(3), x2(3)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    no1 = zi(jdno-1+ino1)
    no2 = zi(jdno-1+ino2)
    do ii = 1, 3
        x1(ii) = zr(jdco+(no1-1)*3+ii-1)
        x2(ii) = zr(jdco+(no2-1)*3+ii-1)
        x3(ii) = x2(ii)-x1(ii)
    end do
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    longseg = sqrt(ddot(b_n, x3, b_incx, x3, b_incy))
!
end subroutine
