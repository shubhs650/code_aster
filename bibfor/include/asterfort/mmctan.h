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
interface
#include "MeshTypes_type.h"
    subroutine mmctan(numema, alias, nno, ndim, coorma,&
                      coorno, epsmax, tau1, tau2, err_appa)
        integer(kind=8) :: numema
        character(len=8) :: alias
        integer(kind=8) :: nno
        integer(kind=8) :: ndim
        real(kind=8) :: coorma(3,MT_NNOMAX2D)
        real(kind=8) :: coorno(3)
        real(kind=8) :: epsmax
        real(kind=8) :: tau1(3)
        real(kind=8) :: tau2(3)
        integer(kind=8), intent(inout) :: err_appa
    end subroutine mmctan
end interface
