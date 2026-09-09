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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine comthm(ds_thm, &
                      lMatr, lSigm, lVari, lMatrPred, &
                      option, typmod, &
                      ndim, nbvari, &
                      dimdef, dimcon, &
                      adcome, adcote, adcp11, adcp12, adcp21, adcp22, adco2nd, &
                      addeme, addete, addep1, addep2, adde2nd, &
                      kpi, npg, &
                      carcri, &
                      defgem, defgep, &
                      congem, congep, &
                      vintm, vintp, &
                      time_prev, time_curr, &
                      dsde, gravity, retcom)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        aster_logical, intent(in) :: lMatr, lSigm, lVari, lMatrPred
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: typmod(2)
        integer(kind=8), intent(in) :: ndim, nbvari
        integer(kind=8), intent(in) :: dimdef, dimcon
        integer(kind=8), intent(in) :: adcome, adcote, adcp11, adcp12, adcp21, adcp22, adco2nd
        integer(kind=8), intent(in) :: addeme, addete, addep1, addep2, adde2nd
        integer(kind=8), intent(in) :: kpi, npg
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: defgem(1:dimdef), defgep(1:dimdef)
        real(kind=8), intent(in) :: congem(1:dimcon)
        real(kind=8), intent(inout) :: congep(1:dimcon)
        real(kind=8), intent(in) :: vintm(1:nbvari)
        real(kind=8), intent(inout) :: vintp(nbvari)
        real(kind=8), intent(in) :: time_prev, time_curr
        real(kind=8), intent(inout) :: dsde(1:dimcon, 1:dimdef)
        real(kind=8), intent(out) :: gravity(3)
        integer(kind=8), intent(out) :: retcom
    end subroutine comthm
end interface
