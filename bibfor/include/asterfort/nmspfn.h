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
    subroutine nmspfn(ndim, nno, nddl, &
                  nno_p, nno_s, npg, &
                  wref, vff_s, vf_p, pgl, geom, nsigm, mate, matpou, &
                  sigma, fint)
        integer(kind=8) :: ndim, nno, nddl, nno_p, nno_s, npg
        integer(kind=8) :: nsigm, mate
        real(kind=8) :: geom(3*nno), wref(npg), vf_p(ndim, npg), vff_s(nno, npg)
        real(kind=8) :: pgl(3, 3)
        real(kind=8) :: fint(nddl)
        real(kind=8) :: sigma(nsigm, npg)
        character(len=8), intent(in)   :: matpou
    end subroutine nmspfn
end interface
