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
! aslint: disable=W0413
! comparaison aver r8gaem et -r8gaem uniquement
!
subroutine pi6036(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                  sigm, vim, dtau, etamin, etamax, copilo)
!
    use Behaviour_type
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/pieigv.h"
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
!  Path-following pred_elas for grad_vari / endo_isot_beton
! --------------------------------------------------------------------------------------------------
    aster_logical:: lBounds
    integer(kind=8):: ndimsi
! --------------------------------------------------------------------------------------------------

    ndimsi = BEHInteg%behavPara%ndimsi
    ASSERT(ndimsi .eq. 4 .or. ndimsi .eq. 6)

    ! Bounds control
    lBounds = etamin .ne. -r8gaem() .and. etamax .ne. r8gaem()
    if (.not. lBounds) call utmess('F', 'MECANONLINE_60', sk='ENDO_ISOT_BETON')

    ! Compute path-following coefficients
    call pieigv(size(epsm), dtau, BEHInteg%materPara%jvMaterCode, &
                vim, epsm, epsd_cste, epsd_pilo, typmod, &
                etamin, etamax, copilo)

end subroutine
