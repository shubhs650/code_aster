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
    subroutine sigtmc(materPara, time, &
                      nbsig, npg, ndim, &
                      indxVarcStrain, sigmVarc)
        use MaterialPara_type
        real(kind=8), intent(in) :: time
        type(Material_Para), intent(inout) :: materPara
        integer(kind=8), intent(in) :: nbsig, npg, ndim
        integer(kind=8), intent(in) :: indxVarcStrain
        real(kind=8), intent(out) :: sigmVarc(nbsig*npg)
    end subroutine sigtmc
end interface
