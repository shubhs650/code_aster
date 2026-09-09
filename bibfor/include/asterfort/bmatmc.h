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
#include "FE_module.h"
!
!
interface
    subroutine bmatmc(igau, nbsig, xyz, ipoids, ivf,&
                      idfde, nno, nharm, jacob, b)
       integer(kind=8), intent(in) :: nbsig, igau, ipoids, ivf, idfde, nno
       real(kind=8), intent(in) :: nharm, xyz(1)
       real(kind=8), intent(out) :: jacob, b(nbsig, MAX_BV_CG)
    end subroutine bmatmc
end interface
