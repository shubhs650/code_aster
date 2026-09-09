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
    subroutine nmorth(materPara, ndim, poum, &
                      dEpsiIn, sigmPrev, option, &
                      sigmCurr, dsidep)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materPara
        integer(kind=8), intent(in) :: ndim
        character(len=*), intent(in) :: poum
        real(kind=8), intent(in) :: dEpsiIn(2*ndim)
        real(kind=8), intent(in) :: sigmPrev(2*ndim)
        character(len=16), intent(in):: option
        real(kind=8), intent(out) :: sigmCurr(2*ndim)
        real(kind=8), intent(out) :: dsidep(2*ndim, 2*ndim)
    end subroutine nmorth
end interface
