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

subroutine tophys_ms(dplmod, psidel, coef, xgene, xphys, nbddl_)
    implicit none
!    Convert into the physical basis some data in generalized coordinates
!
!    Special treatment for multi supported systems with several
!    distinct accelerograms (number = nbexci)
!
!    See tophys.F90 for the regular treatment
!
!-----------------------------------------------------------------------
!
#include "jeveux.h"
!-----------------------------------------------------------------------
    real(kind=8), pointer  :: dplmod(:)
    real(kind=8), pointer  :: psidel(:)
    real(kind=8), pointer  :: coef(:)
    real(kind=8), pointer  :: xgene(:)
    real(kind=8), intent(out) :: xphys(:)
    integer(kind=8), intent(in), optional :: nbddl_
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, nbmode, nbexci, nbddl
!-----------------------------------------------------------------------
    if (present(nbddl_)) then
        nbddl = nbddl_
    else
        nbddl = 3
    end if
!
    nbmode = size(dplmod)/nbddl
    nbexci = size(psidel)/nbddl

    do j = 1, nbddl
        xphys(j) = 0.d0
    end do
!
    do j = 1, nbddl
        do i = 1, nbmode
            xphys(j) = xphys(j)+dplmod((i-1)*nbddl+j)*xgene(i)
        end do
    end do

    do j = 1, nbddl
        do i = 1, nbexci
            xphys(j) = xphys(j)+psidel((i-1)*nbddl+j)*coef(i)
        end do
    end do

end subroutine
