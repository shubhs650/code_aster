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
    subroutine pmstab(sigmPrev, sigmCurr, epsiPrev, epsiIncr, &
                      nbVari, vim, vip, &
                      timePrev, timeCurr, iterNewt, &
                      tablName, tablType, tablNbPara, &
                      tablParaName, tablVale, &
                      lLoadGrad, valeImpo, lPrintMatr, dsidep, variName, &
                      nbVariTabl)
        real(kind=8), intent(in) :: sigmCurr(6), epsiIncr(9)
        real(kind=8), intent(out) :: sigmPrev(6), epsiPrev(9)
        integer(kind=8), intent(in) :: nbVari
        real(kind=8), intent(out) :: vim(nbVari)
        real(kind=8), intent(in) :: vip(nbVari)
        real(kind=8), intent(out) :: timePrev
        real(kind=8), intent(in) :: timeCurr
        integer(kind=8), intent(in) :: iterNewt
        character(len=8), intent(in) :: tablName
        integer(kind=8), intent(in) :: tablType
        integer(kind=8), intent(in) :: tablNbPara
        character(len=16), dimension(:), intent(in) :: tablParaName
        real(kind=8), dimension(:), intent(inout) :: tablVale
        aster_logical, intent(in) :: lLoadGrad
        real(kind=8), intent(in) :: valeImpo(9)
        aster_logical, intent(in) :: lPrintMatr
        real(kind=8), intent(in) :: dsidep(*)
        character(len=8), intent(in) :: variName(nbVari)
    end subroutine pmstab
end interface
