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
    subroutine pipeds(materPara, ndim, typmod, &
                      tau, &
                      vim, epsm, epspc, epsdc, etamin, etamax, &
                      a0, a1, a2, a3, etas)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materPara
        character(len=8), intent(in) :: typmod(2)
        integer(kind=8), intent(in) :: ndim
        real(kind=8) :: tau
        real(kind=8) :: vim(2)
        real(kind=8) :: epsm(6)
        real(kind=8) :: epspc(6)
        real(kind=8) :: epsdc(6)
        real(kind=8) :: etamin
        real(kind=8) :: etamax
        real(kind=8) :: a0
        real(kind=8) :: a1
        real(kind=8) :: a2
        real(kind=8) :: a3
        real(kind=8) :: etas
    end subroutine pipeds
end interface
