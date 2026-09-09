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

interface
    subroutine te0544_implement(typmod, compor, ndim, nno, npg, &
                    jv_poids, jv_vff, jv_dfde, geom_i, &
                    deplm, ddepl, depl0, depl1, dtau, copilo)
#include "asterfort/Behaviour_type.h"
        character(len=8), intent(in):: typmod(:)
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        integer(kind=8) :: ndim, nno, npg
        integer(kind=8) :: jv_poids, jv_vff, jv_dfde
        real(kind=8) :: geom_i(ndim, nno), deplm(ndim, nno), ddepl(ndim, nno)
        real(kind=8) :: depl0(ndim,nno), depl1(ndim,nno)
        real(kind=8) :: dtau
        real(kind=8) :: copilo(5, npg)
    end subroutine
end interface
