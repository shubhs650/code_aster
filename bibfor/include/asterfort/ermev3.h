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
    subroutine ermev3(nno, ipg, ivf, isig, nbcmp,&
                      dfdx, dfdy, dfdz, dsx, dsy,&
                      dsz, norme)
        integer(kind=8) :: nno
        integer(kind=8) :: ipg
        integer(kind=8) :: ivf
        integer(kind=8) :: isig
        integer(kind=8) :: nbcmp
        real(kind=8) :: dfdx(MT_NNOMAX3D)
        real(kind=8) :: dfdy(MT_NNOMAX3D)
        real(kind=8) :: dfdz(MT_NNOMAX3D)
        real(kind=8) :: dsx
        real(kind=8) :: dsy
        real(kind=8) :: dsz
        real(kind=8) :: norme
    end subroutine ermev3
end interface
