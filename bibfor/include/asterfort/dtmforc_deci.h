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
!
interface
    subroutine dtmforc_deci(nl_ind, sd_dtm_, sd_nl_, buffdtm, buffnl,&
                            time  , step   ,depl   , fext)
        integer(kind=8)               , intent(in)  :: nl_ind
        character(len=*)      , intent(in)  :: sd_dtm_
        character(len=*)      , intent(in)  :: sd_nl_
        integer(kind=8)     , pointer  :: buffdtm  (:)
        integer(kind=8)     , pointer  :: buffnl   (:)
        real(kind=8)          , intent(in)  :: time
        real(kind=8)          , intent(in)  :: step
        real(kind=8), pointer  :: depl     (:)
        real(kind=8), pointer :: fext     (:)
    end subroutine dtmforc_deci
end interface
