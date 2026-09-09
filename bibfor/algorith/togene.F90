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

subroutine togene(dplmod, fphys, fgene, coef, nbddl_)
    implicit none
!    Convert into the modal basis a force defined in physical coordinates
!
!    Attention : cumulated fgene, fgene needs thus to be initialized
!                before calling this routine
!-----------------------------------------------------------------------
!
#include "jeveux.h"
!-----------------------------------------------------------------------
    real(kind=8), pointer  :: dplmod(:)
    real(kind=8), intent(in)  :: fphys(:)
    real(kind=8), intent(out) :: fgene(:)
    real(kind=8), optional, intent(in)  :: coef
    integer(kind=8), intent(in), optional :: nbddl_
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, nbmode, nbddl
    real(kind=8) :: coef_m
!-----------------------------------------------------------------------
    coef_m = 1.d0
    if (present(coef)) coef_m = coef
    if (present(nbddl_)) then
        nbddl = nbddl_
    else
        nbddl = 3
    end if

    nbmode = size(dplmod)/nbddl

    do j = 1, nbddl
        do i = 1, nbmode
            fgene(i) = fgene(i)+coef_m*dplmod((i-1)*nbddl+j)*fphys(j)
        end do
    end do

end subroutine
