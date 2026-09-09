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
    subroutine nm1vil(materPara, &
                      relaComp, carcri, &
                      materPoin, &
                      instam, instap, tm, tp, &
                      deps, sigm, vim, &
                      defam, defap, sigp, vip, &
                      dsidep, iret, nbvalc)
        use MaterialPara_type
        type(Material_Para), intent(inout) :: materPara
        character(len=16), intent(in) :: relaComp
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        character(len=*), intent(in) :: materPoin
        integer(kind=8) :: iret, nbvalc
        real(kind=8) :: instam, instap
        real(kind=8) :: tm, tp
        real(kind=8) :: irram, irrap
        real(kind=8) :: deps
        real(kind=8) :: sigm, vim(nbvalc)
        real(kind=8) :: defam, defap
        real(kind=8) :: sigp, vip(nbvalc), dsidep, alpha
    end subroutine nm1vil
end interface
