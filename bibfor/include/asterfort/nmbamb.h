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
    subroutine nmbamb(ndim, nno, npg, geom, section, &
                    dxi_ff, w_ref, nddl, neps, b, &
                    w, ni2ldc)
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: nno
        integer(kind=8), intent(in) :: npg
        real(kind=8), intent(in)    :: geom(ndim, nno)
        real(kind=8), intent(in)    :: dxi_ff(nno,npg)
        real(kind=8), intent(in)    :: w_ref(npg)
        real(kind=8), intent(in)    :: section
        integer(kind=8), intent(out):: nddl
        integer(kind=8), intent(out):: neps
        real(kind=8), intent(out), allocatable :: b(:, :, :)
        real(kind=8), intent(out), allocatable :: w(:, :)
        real(kind=8), intent(out), allocatable :: ni2ldc(:, :)
    end subroutine nmbamb
end interface
