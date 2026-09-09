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
    subroutine massup(jvMaterCode, &
                      option, ndim, dlns, nno, nnos, &
                      npg, ipoids, idfde, &
                      geom, vff1, imatuu, icodre, igeom, &
                      ivf)
        character(len=16) :: option
        integer(kind=8) :: ndim, nno, nnos, npg
        integer(kind=8) :: dlns
        integer(kind=8) :: ipoids, jvMaterCode
        integer(kind=8) :: idfde
        real(kind=8) :: geom(ndim, nno)
        real(kind=8) :: vff1(nno, npg)
        integer(kind=8) :: imatuu
        integer(kind=8) :: icodre(1)
        integer(kind=8) :: igeom
        integer(kind=8) :: ivf
    end subroutine massup
end interface
