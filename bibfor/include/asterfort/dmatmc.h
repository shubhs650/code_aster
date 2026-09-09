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
    subroutine dmatmc(materPara, poum, time, &
                      tensSize, dr_, l_modi_cp, di_)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materPara
        character(len=*), intent(in) :: poum
        real(kind=8), intent(in) :: time
        integer(kind=8), intent(in) :: tensSize
        real(kind=8), optional, intent(out) :: dr_(tensSize, tensSize), di_(tensSize, tensSize)
        aster_logical, optional, intent(in) :: l_modi_cp
    end subroutine dmatmc
end interface
