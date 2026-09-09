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
    subroutine nmveei(materPara, &
                      carcri, compor, ndim, typmod, &
                      instam, instap, &
                      epsm, deps, sigm, nvi, vim, option, &
                      sigp, vip, dsidep, iret)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materPara
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        character(len=8), intent(in) :: typmod(2)
        integer(kind=8), intent(in) :: nvi, ndim
        character(len=16), intent(in) ::  option
        real(kind=8) :: instam, instap
        real(kind=8) :: epsm(6)
        real(kind=8) :: deps(6)
        real(kind=8) :: sigm(6)
        real(kind=8) :: vim(*)
        real(kind=8) :: sigp(6)
        real(kind=8) :: vip(*)
        real(kind=8) :: dsidep(6, 6)
        integer(kind=8) :: iret
    end subroutine nmveei
end interface
