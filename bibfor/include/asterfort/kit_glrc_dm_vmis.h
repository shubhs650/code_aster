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
    subroutine kit_glrc_dm_vmis(plateCara, plateOrie, &
                                jvMaterCode, compor, epsm, deps, vim, &
                                option, sigm, sig, vip, dsidep, &
                                carcri, iret)
        use plate_type
        type(plateCara_Para), intent(in) :: plateCara
        type(plateOrie_Para), intent(in) :: plateOrie
        integer(kind=8) :: imate
        character(len=16) :: compor
        real(kind=8) :: epsm(6)
        real(kind=8) :: deps(6)
        real(kind=8) :: vim(*)
        character(len=16) :: option
        real(kind=8) :: sigm(*)
        real(kind=8) :: sig(*)
        real(kind=8) :: vip(*)
        real(kind=8) :: dsidep(6, *)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        integer(kind=8) :: iret
    end subroutine kit_glrc_dm_vmis
end interface
