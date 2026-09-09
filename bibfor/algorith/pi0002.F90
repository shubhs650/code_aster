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
subroutine pi0002(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                  sigm, vim, dtau, etamin, etamax, copilo, relaComp)
!
    use Behaviour_type
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/pipepl.h"
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
    character(len=16), intent(in):: relaComp
! --------------------------------------------------------------------------------------------------
!  Path-following pred_elas for vmis_isot_line and vmis_isot_trac
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: ndimsi
    real(kind=8):: epsm_6(6), epsp_6(6), epsd_6(6), sigm_6(6)
! --------------------------------------------------------------------------------------------------

    ndimsi = BEHInteg%behavPara%ndimsi
    ASSERT(ndimsi .eq. 4 .or. ndimsi .eq. 6)
    ASSERT(relaComp .eq. 'VMIS_ISOT_LINE' .or. relaComp .eq. 'VMIS_ISOT_TRAC')

    ! Larger arrays
    epsm_6 = 0
    epsp_6 = 0
    epsd_6 = 0
    sigm_6 = 0
    epsm_6(1:ndimsi) = epsm
    epsp_6(1:ndimsi) = epsm+epsd_cste
    epsd_6(1:ndimsi) = epsd_pilo
    sigm_6(1:ndimsi) = sigm

    ! Compute path-following coefficients
    call pipepl(BEHInteg%materPara, ndim, relaComp, typmod, &
                dtau, sigm_6, vim, epsp_6-epsm_6, epsd_6, &
                copilo(1), copilo(2), copilo(3), copilo(4), copilo(5))

end subroutine
