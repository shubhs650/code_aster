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
    subroutine pi0000(BEHInteg, compor, typmod, ndim, &
                  epsm, epsd_cste, epsd_pilo, &
                  sigm, vim, dtau, etamin, etamax, copilo)

        use Behaviour_type
#include "asterfort/Behaviour_type.h"

    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    character(len=8), intent(in) :: typmod(2)
    integer(kind=8), intent(in):: ndim
    real(kind=8), intent(in) :: epsm(:)
    real(kind=8), intent(in) :: epsd_cste(:)
    real(kind=8), intent(in) :: epsd_pilo(:)
    real(kind=8), intent(in) :: sigm(:)
    real(kind=8), intent(in) :: vim(:)
    real(kind=8), intent(in) :: dtau
    real(kind=8), intent(in) :: etamin
    real(kind=8), intent(in) :: etamax
    real(kind=8), intent(out) :: copilo(:)
    end subroutine
end interface
