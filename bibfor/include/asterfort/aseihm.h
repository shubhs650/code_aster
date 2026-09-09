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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine aseihm(ds_thm, option, &
                      lSigm, lVari, lMatr, lVect, &
                      l_axi, ndim, nno1, nno2, &
                      npi, npg, dimuel, dimdef, dimcon, &
                      nbvari, iu, ip, ipf, &
                      iq, mecani, press1, press2, tempe, &
                      vff1, vff2, dffr2, time_prev, time_curr, &
                      deplm, deplp, sigm, sigp, varim, &
                      varip, wref, geom, &
                      compor, vectu, matuu, &
                      retcom)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        aster_logical, intent(in) :: lSigm, lVari, lMatr, lVect
        integer(kind=8) :: nbvari
        integer(kind=8) :: dimcon
        integer(kind=8) :: dimdef
        integer(kind=8) :: dimuel
        integer(kind=8) :: npi
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        integer(kind=8) :: ndim
        character(len=16) :: option
        aster_logical :: l_axi
        integer(kind=8) :: npg
        integer(kind=8) :: iu(3, 18)
        integer(kind=8) :: ip(2, 9)
        integer(kind=8) :: ipf(2, 2, 9)
        integer(kind=8) :: iq(2, 2, 9)
        integer(kind=8) :: mecani(8)
        integer(kind=8) :: press1(9)
        integer(kind=8) :: press2(9)
        integer(kind=8) :: tempe(5)
        real(kind=8) :: vff1(nno1, npi)
        real(kind=8) :: vff2(nno2, npi)
        real(kind=8) :: dffr2(ndim-1, nno2, npi)
        real(kind=8) :: time_prev
        real(kind=8) :: time_curr
        real(kind=8) :: deplm(dimuel)
        real(kind=8) :: deplp(dimuel)
        real(kind=8) :: sigm(dimcon, npi)
        real(kind=8) :: sigp(dimcon, npi)
        real(kind=8) :: varim(nbvari, npi)
        real(kind=8) :: varip(nbvari, npi)
        real(kind=8) :: wref(npi)
        real(kind=8) :: geom(ndim, nno2)
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8) :: vectu(dimuel)
        real(kind=8) :: matuu(dimuel*dimuel)
        integer(kind=8) :: retcom
    end subroutine aseihm
end interface
