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

function lonelesp(geom, nno, no1, no2)
!
!
! --------------------------------------------------------------------------------------------------
!
!                           CALCULE LA LONGEUR DES ELEMENTS DE POUTRE
!                                    DES ELEMENTS 3D_SOLPIEU
!
!   OUT
!       lonelesp   : longueur de la poutre
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
!
    integer(kind=8)  :: nno, no1, no2
    real(kind=8)     :: geom(3, nno)
    real(kind=8)     :: lonelesp
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)     :: i
    real(kind=8)        :: xl
!
! --------------------------------------------------------------------------------------------------
!
    xl = 0.d0
    do i = 1, 3
        xl = xl+(geom(i, no1)-geom(i, no2))**2
    end do
    lonelesp = sqrt(xl)
!
end function
