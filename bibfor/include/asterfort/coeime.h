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
    subroutine coeime(ds_thm, &
                      lSigm, lVari, lMatr, &
                      option, &
                      ndim, dimdef, dimcon, &
                      addeme, addep1, &
                      nbvari, npg, npi, &
                      defgep, defgem, sigm, sigp, varim, &
                      varip, ouvh, tlint, drde, kpi, &
                      retcom)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        character(len=16), intent(in) :: option
        aster_logical, intent(in) :: lSigm, lVari, lMatr
        integer(kind=8), intent(in) :: ndim, dimcon, dimdef
        integer(kind=8), intent(in) :: addeme, addep1, npg, kpi, npi, nbvari
        real(kind=8), intent(in) :: defgem(dimdef), defgep(dimdef)
        real(kind=8), intent(in) :: sigm(dimcon)
        real(kind=8), intent(inout) :: sigp(dimcon)
        real(kind=8), intent(in) :: varim(nbvari)
        real(kind=8), intent(inout) :: varip(nbvari)
        real(kind=8), intent(out) :: ouvh, tlint
        real(kind=8), intent(inout) :: drde(dimdef, dimdef)
        integer(kind=8), intent(out) :: retcom
    end subroutine coeime
end interface
