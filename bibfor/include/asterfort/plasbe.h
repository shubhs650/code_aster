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
    subroutine plasbe(BEHinteg, &
                      fami, kpg, ksp, typmod, imat, l_epsi_varc, &
                      carcri, epsdt, depst, sigd, nvi, vind, &
                      opt, sigf, vinf, dsde, &
                      irteti)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        aster_logical, intent(in) :: l_epsi_varc
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=8) :: typmod(*)
        integer(kind=8) :: imat, nvi
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8) :: epsdt(6)
        real(kind=8) :: depst(6)
        real(kind=8) :: sigd(6)
        real(kind=8) :: vind(nvi)
        character(len=16) :: opt
        real(kind=8) :: sigf(6)
        real(kind=8) :: vinf(nvi)
        real(kind=8) :: dsde(6, 6)
        integer(kind=8) :: irteti
    end subroutine plasbe
end interface
