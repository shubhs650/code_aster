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
subroutine pi9041(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                  sigm, vim, dtau, etamin, etamax, copilo)
!
    use Behaviour_type
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/pipetc.h"
#include "asterfort/utmess.h"

    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=8), intent(in) :: typmod(2)
    integer(kind=8) :: ndim
    real(kind=8) :: epsm(:)
    real(kind=8) :: epsd_cste(:)
    real(kind=8) :: epsd_pilo(:)
    real(kind=8) :: sigm(:)
    real(kind=8) :: vim(:)
    real(kind=8) :: dtau
    real(kind=8) :: etamin
    real(kind=8) :: etamax
    real(kind=8), intent(out) :: copilo(:)
! --------------------------------------------------------------------------------------------------
!  Path-following pred_elas for endo_orth_beton
! --------------------------------------------------------------------------------------------------
    real(kind=8):: mu_cste(3), su_cste(3), mu_pilo(3), su_pilo(3)
! --------------------------------------------------------------------------------------------------

    ! Controls
    ASSERT(size(epsm) .ge. 2*ndim)

    ! Larger arrays
    su_cste = 0
    su_pilo = 0
    mu_cste = 0
    mu_pilo = 0
    su_cste(1:ndim) = epsm(1:ndim)+epsd_cste(1:ndim)
    su_pilo(1:ndim) = epsd_pilo(1:ndim)
    mu_cste(1:ndim) = epsm(ndim+1:2*ndim)+epsd_cste(ndim+1:2*ndim)
    mu_pilo(1:ndim) = epsd_pilo(ndim+1:2*ndim)

    ! Compute path-following coefficients
    call pipetc(BEHInteg%materPara%jvMaterCode, &
                su_cste, su_pilo, mu_cste, mu_pilo, &
                vim, dtau, copilo)

end subroutine
