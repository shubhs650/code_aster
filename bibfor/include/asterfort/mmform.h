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
#include "MeshTypes_type.h"
interface
    subroutine mmform(ndim  ,&
                      typmae, typmam,&
                      nne   , nnm   ,&
                      xpc   , ypc   , xpr  , ypr,&
                      ffe   , dffe  , ddffe,&
                      ffm   , dffm  , ddffm,&
                      ffl   , dffl  , ddffl)
        integer(kind=8), intent(in) :: ndim
        character(len=8), intent(in) :: typmae, typmam
        integer(kind=8), intent(in) :: nne, nnm
        real(kind=8), intent(in) :: xpc, ypc, xpr, ypr
        real(kind=8), intent(out) :: ffe(9), dffe(2, MT_NNOMAX2D), ddffe(3, MT_NNOMAX2D)
        real(kind=8), intent(out) :: ffm(9), dffm(2, MT_NNOMAX2D), ddffm(3, MT_NNOMAX2D)
        real(kind=8), intent(out) :: ffl(9), dffl(2, MT_NNOMAX2D), ddffl(3, MT_NNOMAX2D)
    end subroutine mmform
end interface
