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
    subroutine matthm(ds_thm, &
                      ndim, axi, nno1, nno2, dimuel, &
                      dimdef, iu, ip, ipf, iq, &
                      addep1, &
                      addlh1, vff1, vff2, dffr2, wref, &
                      geom, wi, q)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        integer(kind=8) :: dimdef
        integer(kind=8) :: dimuel
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        integer(kind=8) :: ndim
        aster_logical :: axi
        integer(kind=8) :: iu(3, 18)
        integer(kind=8) :: ip(2, 9)
        integer(kind=8) :: ipf(2, 2, 9)
        integer(kind=8) :: iq(2, 2, 9)
        integer(kind=8) :: addep1
        integer(kind=8) :: addlh1
        real(kind=8) :: vff1(nno1)
        real(kind=8) :: vff2(nno2)
        real(kind=8) :: dffr2(ndim-1, nno2)
        real(kind=8) :: wref
        real(kind=8) :: geom(ndim, nno2)
        real(kind=8) :: wi
        real(kind=8) :: q(dimdef, dimuel)
    end subroutine matthm
end interface
