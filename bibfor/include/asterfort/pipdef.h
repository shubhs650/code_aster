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
    subroutine pipdef(typmod, &
                    ndim, nno, kpg, jv_poids, jv_vff, &
                    jv_dfde, geom, deplm, &
                    ddepl, depl0, depl1, &
                    epsm, deps_cst, deps_pil)
        character(len=8), intent(in) :: typmod(:)
        integer(kind=8), intent(in) :: ndim, nno, kpg
        integer(kind=8),intent(in) :: jv_poids, jv_vff, jv_dfde
        real(kind=8),intent(in) :: geom(:,:), deplm(:,:), ddepl(:,:), depl0(:,:), depl1(:,:)
        real(kind=8), intent(out) :: epsm(:), deps_cst(:), deps_pil(:)
    end subroutine pipdef
end interface
