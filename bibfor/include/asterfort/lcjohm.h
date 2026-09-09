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
    subroutine lcjohm(materParaFPG1, &
                      lSigm, lMatr, lVari, &
                      kpg, npg, &
                      addeme, advico, ndim, dimdef, &
                      dimcon, nbvari, defgem, defgep, varim, &
                      varip, sigm, sigp, drde, ouvh, &
                      retcom)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materParaFPG1
        integer(kind=8) :: nbvari
        integer(kind=8) :: dimcon
        integer(kind=8) :: dimdef
        aster_logical, intent(in) :: lSigm, lMatr, lVari
        integer(kind=8) :: kpg
        integer(kind=8) :: jvMaterCode
        integer(kind=8) :: npg
        integer(kind=8) :: addeme
        integer(kind=8) :: advico
        integer(kind=8) :: ndim
        real(kind=8) :: defgem(dimdef)
        real(kind=8) :: defgep(dimdef)
        real(kind=8) :: varim(nbvari)
        real(kind=8) :: varip(nbvari)
        real(kind=8) :: sigm(dimcon)
        real(kind=8) :: sigp(dimcon)
        real(kind=8) :: drde(dimdef, dimdef)
        real(kind=8) :: ouvh
        integer(kind=8) :: retcom
    end subroutine lcjohm
end interface
