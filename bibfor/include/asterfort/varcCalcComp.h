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
    subroutine varcCalcComp(modelZ, chsithz, &
                            l_temp, l_hydr, l_ptot, &
                            l_sech, l_epsa, &
                            nbFieldIn, nbFieldOut, &
                            lpain, lchin, &
                            lpaout, lchout, &
                            vectElemZ)
        character(len=*), intent(in) :: modelZ, chsithz
        aster_logical, intent(in)  :: l_temp, l_hydr, l_ptot, l_sech, l_epsa
        integer(kind=8), intent(in) :: nbFieldIn, nbFieldOut
        character(len=8), intent(in) :: lpain(*), lpaout(*)
        character(len=19), intent(in) :: lchin(*)
        character(len=19), intent(inout) :: lchout(*)
        character(len=*), intent(in) :: vectElemZ
    end subroutine varcCalcComp
end interface
