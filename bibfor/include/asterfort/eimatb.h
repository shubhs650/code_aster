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

#include "asterf_types.h"
!
interface
    subroutine eimatb(nomte, ndim, axi, nno1, nno2, npg, &
                      wref, vff1, vff2, dffr2, geom, &
                      anglNautPg, b, wg, ni2ldc)
        character(len=16) :: nomte
        aster_logical, intent(in) :: axi
        integer(kind=8), intent(in) :: ndim, nno1, nno2, npg
        real(kind=8), intent(in) :: vff1(nno1, npg), vff2(nno2, npg), geom(ndim, nno2)
        real(kind=8), intent(in) :: wref(npg)
        real(kind=8), intent(in) :: dffr2(ndim-1, nno2, npg)
        real(kind=8), intent(in) :: anglNautPg(merge(1, 3, ndim .eq. 2), nno2)
        real(kind=8), intent(out) :: b(2*ndim, npg, 2*ndim*nno1+ndim*nno2)
        real(kind=8), intent(out) :: wg(2*ndim, npg)
        real(kind=8), intent(out) :: ni2ldc(2*ndim, npg)
    end subroutine
end interface
