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
    subroutine te0550_implement(option, fami, nno, npg, ndim_sp, &
        wref, vff, dxi_ff, aire, geom, pesa, &
        mate, fext)

    character(len=*) :: option
    character(len=8) :: fami
    integer(kind=8), intent(in):: nno, npg, ndim_sp, mate
    real(kind=8), intent(in):: geom(ndim_sp, nno), wref(npg), vff(nno,npg), dxi_ff(nno,npg)
    real(kind=8), intent(in):: aire, pesa(0:ndim_sp)
    real(kind=8), intent(out):: fext(ndim_sp, nno)

    end subroutine te0550_implement
end interface
