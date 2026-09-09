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
    subroutine radipg(sig1, sig2, npg, nbsig, radia, &
                      cosang, ind, relaComp, imate, nvi, &
                      vari1, vari2)
        real(kind=8) :: sig1(*)
        real(kind=8) :: sig2(*)
        integer(kind=8) :: npg
        integer(kind=8) :: nbsig
        real(kind=8) :: radia(*)
        real(kind=8) :: cosang(*)
        integer(kind=8) :: ind
        character(len=16) :: relaComp
        integer(kind=8) :: imate
        integer(kind=8) :: nvi
        real(kind=8) :: vari1(*)
        real(kind=8) :: vari2(*)
    end subroutine radipg
end interface
