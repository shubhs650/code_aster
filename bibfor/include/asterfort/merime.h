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
    subroutine merime(modelz, nbLoad, listLoadK24, &
                      materFieldZ, materCodeZ, caraElemz, &
                      time, comporMultz, matrElemz, modeFourier, &
                      jvBaseZ, listElemCalcz, &
                      hasExteStatVari_, onlyDirichlet_)
        character(len=*), intent(in) :: modelz
        integer(kind=8), intent(in) :: nbLoad
        character(len=24), pointer :: listLoadK24(:)
        character(len=*), intent(in) :: materFieldZ, materCodeZ, caraElemz
        real(kind=8), intent(in) :: time
        character(len=*), intent(in) :: comporMultz, matrElemz
        integer(kind=8), intent(in) :: modeFourier
        character(len=*), intent(in) :: jvBaseZ, listElemCalcz
        aster_logical, intent(in), optional :: hasExteStatVari_, onlyDirichlet_
    end subroutine merime
end interface
