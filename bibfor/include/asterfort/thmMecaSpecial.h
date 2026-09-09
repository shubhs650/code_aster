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
    subroutine thmMecaSpecial(ds_thm, option, lMatr, meca, &
                              nl, &
                              ndim, typmod, carcri, &
                              addeme, adcome, addep1, addep2, &
                              dimdef, dimcon, &
                              defgem, deps, &
                              congem, vintm, &
                              congep, vintp, &
                              time_prev, time_curr, &
                              dsde, ther_meca, retcom)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        character(len=16), intent(in) :: option, meca
        aster_logical, intent(in) :: lMatr
        real(kind=8), intent(in) :: nl
        character(len=8), intent(in) :: typmod(2)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        integer(kind=8), intent(in) :: ndim, dimdef, dimcon, addeme, adcome, addep1, addep2
        real(kind=8), intent(in) :: vintm(*)
        real(kind=8), intent(in) :: defgem(dimdef), deps(6), congem(dimcon)
        real(kind=8), intent(inout) :: congep(dimcon)
        real(kind=8), intent(inout) :: vintp(*)
        real(kind=8), intent(in) :: time_prev, time_curr
        real(kind=8), intent(inout) :: dsde(dimcon, dimdef)
        real(kind=8), intent(out) :: ther_meca(6)
        integer(kind=8), intent(out) :: retcom
    end subroutine thmMecaSpecial
end interface
