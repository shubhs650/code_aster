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
    subroutine mbgchg(plateOrie, &
                      option, fami, &
                      nddl, nno, ncomp, kpg, &
                      jvMaterc, jvSief, &
                      ipoids, jvPesa, jvGeom, jvVect, &
                      vff, dff, &
                      h, preten)
        use plate_type
        type(plateOrie_Para), intent(in) :: plateOrie
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: fami
        integer(kind=8), intent(in) :: nddl, nno, ncomp
        integer(kind=8), intent(in) :: kpg
        integer(kind=8), intent(in) :: ipoids, jvGeom, jvMaterc, jvPesa
        integer(kind=8), intent(in) :: jvVect, jvSief
        real(kind=8), intent(in) :: dff(2, nno), vff(nno)
        real(kind=8), intent(in) :: h, preten
    end subroutine mbgchg
end interface
