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

subroutine pidefo(epsm, epsd_cste, epsd_pilo, dtau, copilo)

    implicit none
#include "asterc/r8prem.h"
#include "asterfort/pinorm.h"

    real(kind=8), intent(in) :: epsm(:)
    real(kind=8), intent(in) :: epsd_cste(:)
    real(kind=8), intent(in) :: epsd_pilo(:)
    real(kind=8), intent(in) :: dtau
    real(kind=8), intent(out) :: copilo(:)
!
! --------------------------------------------------------------------------------------------------
!
! routine meca_non_line (pilotage)
!
! calcul des coefficients de pilotage pour pred_elas/deformation
!
! --------------------------------------------------------------------------------------------------
! in  epsm  : déformation en t-
! in  epsd_cste  : incrément de déformation pour forces fixes
! in  epsd_pilo  : incrément de déformation pour forces pilotees
! out copilo : coefficients de pilotage
! --------------------------------------------------------------------------------------------------
    real(kind=8):: epsmno
! --------------------------------------------------------------------------------------------------

    epsmno = norm2(epsm)

    if (epsmno .gt. r8prem()) then
        ! Usual path-following function dot(epsm, epsd) = dtau
        copilo(1) = dot_product(epsm, epsd_cste)/epsmno
        copilo(2) = dot_product(epsm, epsd_pilo)/epsmno
    else
        ! Mode failsafe for zero initial strain -> Solution or norm(deps)=dtau
        call pinorm(epsd_pilo, epsd_cste, dtau, dtau, copilo)
    end if

end subroutine
