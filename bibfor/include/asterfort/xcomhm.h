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
interface
    subroutine xcomhm(ds_thm, &
                      option, time_curr, &
                      ndim, dimdef, dimcon, nbvari, &
                      addeme, adcome, addep1, adcp11, &
                      addep2, addete, defgem, &
                      defgep, congem, congep, vintm, &
                      vintp, dsde, gravity, retcom, &
                      kpg, npg, dimenr, &
                      yaenrh, adenhy, nfh)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        integer(kind=8) :: dimenr
        integer(kind=8) :: nbvari
        integer(kind=8) :: dimcon
        integer(kind=8) :: dimdef
        integer(kind=8) :: ndim
        character(len=16) :: option
        real(kind=8) :: time_curr
        integer(kind=8) :: addeme
        integer(kind=8) :: adcome
        integer(kind=8) :: addep1
        integer(kind=8) :: adcp11
        integer(kind=8) :: addep2
        integer(kind=8) :: addete
        real(kind=8) :: defgem(1:dimdef)
        real(kind=8) :: defgep(1:dimdef)
        real(kind=8) :: congem(1:dimcon)
        real(kind=8) :: congep(1:dimcon)
        real(kind=8) :: vintm(1:nbvari)
        real(kind=8) :: vintp(1:nbvari)
        real(kind=8) :: dsde(1:dimcon, 1:dimenr)
        real(kind=8) :: gravity(3)
        integer(kind=8) :: retcom
        integer(kind=8) :: kpg, npg
        integer(kind=8) :: idecpg
        integer(kind=8) :: yaenrh
        integer(kind=8) :: adenhy
        integer(kind=8) :: nfh
    end subroutine xcomhm
end interface
