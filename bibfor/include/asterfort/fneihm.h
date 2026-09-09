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
    subroutine fneihm(ds_thm, &
                      nno1, nno2, &
                      npi, npg, wref, iu, ip, &
                      ipf, iq, vff1, vff2, dffr2, &
                      geom, congem, r, vectu, &
                      mecani, press1, press2, dimdef, &
                      dimcon, dimuel, ndim, axi)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        integer(kind=8) :: ndim
        integer(kind=8) :: dimuel
        integer(kind=8) :: dimcon
        integer(kind=8) :: dimdef
        integer(kind=8) :: npg
        integer(kind=8) :: npi
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        real(kind=8) :: wref(npg)
        integer(kind=8) :: iu(3, 18)
        integer(kind=8) :: ip(2, 9)
        integer(kind=8) :: ipf(2, 2, 9)
        integer(kind=8) :: iq(2, 2, 9)
        real(kind=8) :: vff1(nno1, npi)
        real(kind=8) :: vff2(nno2, npi)
        real(kind=8) :: dffr2(ndim-1, nno2, npi)
        real(kind=8) :: geom(ndim, nno2)
        real(kind=8) :: congem(dimcon, npi)
        real(kind=8) :: r(dimdef)
        real(kind=8) :: vectu(dimuel)
        integer(kind=8) :: mecani(8)
        integer(kind=8) :: press1(9)
        integer(kind=8) :: press2(9)
        aster_logical :: axi
    end subroutine fneihm
end interface
