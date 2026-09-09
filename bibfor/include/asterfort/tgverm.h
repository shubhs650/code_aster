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
#include "asterfort/Behaviour_type.h"
#include "asterf_types.h"
!
interface
    subroutine tgverm(option, carcri, compor, nno1, nno2, &
                      nno3, geom, ndim, nddl, deplp, &
                      sdepl, vu, vg, vp, vectu, &
                      svect, ncont, contp, scont, nvari, &
                      varip, svari, matuu, smatr, matsym, &
                      epsilo, epsilp, epsilg, varia, iret)
        integer(kind=8) :: nddl
        character(len=16) :: option
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        integer(kind=8) :: nno1
        integer(kind=8) :: nno2
        integer(kind=8) :: nno3
        real(kind=8) :: geom(*)
        integer(kind=8) :: ndim
        real(kind=8) :: deplp(*)
        real(kind=8) :: sdepl(*)
        integer(kind=8) :: vu(3, 27)
        integer(kind=8) :: vg(27)
        integer(kind=8) :: vp(27)
        real(kind=8) :: vectu(*)
        real(kind=8) :: svect(*)
        integer(kind=8) :: ncont
        real(kind=8) :: contp(*)
        real(kind=8) :: scont(*)
        integer(kind=8) :: nvari
        real(kind=8) :: varip(*)
        real(kind=8) :: svari(*)
        real(kind=8) :: matuu(*)
        real(kind=8) :: smatr(*)
        real(kind=8) :: varia(*)
        aster_logical :: matsym
        real(kind=8) :: epsilo
        real(kind=8) :: epsilp
        real(kind=8) :: epsilg
        integer(kind=8) :: iret
    end subroutine tgverm
end interface
