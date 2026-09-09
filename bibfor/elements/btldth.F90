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
subroutine btldth(nb1, btild, wgt, &
                  hasTemp, tempKpg, &
                  young, nu, alpha, &
                  forthi)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nb1
    real(kind=8), intent(in) :: btild(5, 42), wgt
    aster_logical, intent(in) :: hasTemp
    real(kind=8), intent(in) :: tempKpg
    real(kind=8), intent(in) :: young, nu, alpha
    real(kind=8), intent(out) :: forthi(42)
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Get temperature and dilatation coefficient at current integration point
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, k
    real(kind=8) :: vecthr(2)
!
! --------------------------------------------------------------------------------------------------
!
    forthi = 0.d0
    if (hasTemp) then
        vecthr(1) = young*alpha*tempKpg/(1.d0-nu)
        vecthr(2) = vecthr(1)
        do i = 1, 5*nb1+2
            forthi(i) = 0.d0
            do k = 1, 2
                forthi(i) = forthi(i)+btild(k, i)*vecthr(k)*wgt
            end do
        end do
    end if
!
end subroutine
