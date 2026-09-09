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
    subroutine pmsta1(sigmPrev, sigmCurr, epsiIncr, &
                      nbVari, nbVariTabl, &
                      vim, vip, &
                      tablNbPara, tablType, &
                      tablParaName, tablParaType, tablVale, &
                      lLoadGrad, variName, sddisc, &
                      liccvg, lIterNewtMaxi, conver, newtLoopAction)
        real(kind=8), intent(in) :: sigmPrev(6), sigmCurr(6), epsiIncr(9)
        integer(kind=8), intent(in)  :: nbVari, nbVariTabl
        real(kind=8), intent(in) :: vim(nbVari), vip(nbVari)
        integer(kind=8), intent(in) :: tablNbPara, tablType
        character(len=16), dimension(:), intent(in) :: tablParaName
        character(len=8), dimension(:), intent(in) :: tablParaType
        real(kind=8), dimension(:), intent(inout) :: tablVale
        aster_logical, intent(in) :: lLoadGrad
        character(len=8), intent(in) :: variName(nbVari)
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: liccvg(5)
        aster_logical, intent(in) :: lIterNewtMaxi, conver
        integer(kind=8), intent(out) :: newtLoopAction
    end subroutine pmsta1
end interface
