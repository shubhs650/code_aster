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
#include "asterf_types.h"
!
interface
    subroutine nmspse(ndim, nno, nddl, &
                  nno_p, nno_s, nddl_s, npg, &
                  vff_s, vf_p, pgl, geom, neps, mate, matpou, &
                  deplm, urpg)
        integer(kind=8) :: ndim
        integer(kind=8) :: nno
        integer(kind=8) :: nddl
        integer(kind=8) :: nno_p
        integer(kind=8) :: nno_s
        integer(kind=8) :: nddl_s
        integer(kind=8) :: npg
        real(kind=8) :: vf_p(ndim, npg)
        real(kind=8) :: vff_s(nno, npg)
        real(kind=8) :: pgl(3, 3)
        real(kind=8) :: geom(3*nno)
        integer(kind=8) :: neps
        integer(kind=8) :: mate
        character(len=8)  :: matpou
        real(kind=8) :: deplm(nddl)
        real(kind=8) :: urpg(neps, npg)
    end subroutine nmspse
end interface
