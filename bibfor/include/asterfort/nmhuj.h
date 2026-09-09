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
!
interface
    subroutine nmhuj(BEHInteg, &
                     fami, kpg, ksp, typmod, jvMaterCode, &
                     carcri, epsd, &
                     deps, sigd, vind, option, sigf, &
                     vinf, dsde, iret)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=8) :: typmod(2)
        integer(kind=8) :: jvMaterCode
        real(kind=8) :: carcri(CARCRI_SIZE)
        real(kind=8) :: epsd(6)
        real(kind=8) :: deps(6)
        real(kind=8) :: sigd(6)
        real(kind=8) :: vind(50)
        character(len=16) :: option
        real(kind=8) :: sigf(6)
        real(kind=8) :: vinf(50)
        real(kind=8) :: dsde(6, 6)
        integer(kind=8) :: iret
    end subroutine nmhuj
end interface
