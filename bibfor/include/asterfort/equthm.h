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
    subroutine equthm(ds_thm, option, &
                      lMatr, lSigm, &
                      lVari, lMatrPred, &
                      typmod, parm_theta, &
                      ndim, nbvari, &
                      kpi, npg, &
                      dimdef, dimcon, &
                      mecani, press1, press2, tempe, second, &
                      carcri, &
                      defgem, defgep, &
                      congem, congep, &
                      vintm, vintp, &
                      time_prev, time_curr, time_incr, &
                      r, drds, dsde, retcom)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        character(len=16), intent(in) :: option
        aster_logical, intent(in) :: lMatr, lSigm, lVari, lMatrPred
        character(len=8), intent(in) :: typmod(2)
        real(kind=8), intent(in)  :: parm_theta
        integer(kind=8), intent(in) :: ndim, nbvari
        integer(kind=8), intent(in) :: npg, kpi
        integer(kind=8), intent(in) :: dimdef, dimcon
        integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: defgem(dimdef), defgep(dimdef)
        real(kind=8), intent(inout) :: congem(dimcon), congep(dimcon)
        real(kind=8), intent(in) :: vintm(nbvari)
        real(kind=8), intent(inout) :: vintp(nbvari)
        real(kind=8), intent(in) :: time_prev, time_curr, time_incr
        real(kind=8), intent(out) :: r(dimdef+1)
        real(kind=8), intent(out) :: drds(dimdef+1, dimcon), dsde(dimcon, dimdef)
        integer(kind=8), intent(out) :: retcom
    end subroutine equthm
end interface
