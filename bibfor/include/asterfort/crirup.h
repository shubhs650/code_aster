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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine crirup(materPara, &
                      ndim, npg, lgpg, &
                      option, sigp, vip, vim, &
                      instam, instap)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materPara
        character(len=16), intent(in) :: option
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: ndim
        real(kind=8) :: sigp(2*ndim, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: instam
        real(kind=8) :: instap
    end subroutine crirup
end interface
