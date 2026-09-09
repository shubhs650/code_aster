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
    subroutine resi_ther(l_stat, &
                         modelZ, caraElemZ, materCodeZ, &
                         timePara, timeMapZ, varcCurrZ, &
                         comporTherZ, tempIterZ, &
                         tempPrevZ, hydrPrevZ, hydrCurrZ, &
                         resuElemZ, vectElemZ, jvBase)
        aster_logical, intent(in) :: l_stat
        character(len=*), intent(in) :: modelZ, caraElemZ, materCodeZ
        real(kind=8), intent(in) :: timePara(2)
        character(len=*), intent(in) :: tempIterZ, comporTherZ, varcCurrZ
        character(len=*), intent(in) :: tempPrevZ, hydrPrevZ, hydrCurrZ, timeMapZ
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
        character(len=1), intent(in) :: jvBase
    end subroutine resi_ther
end interface
