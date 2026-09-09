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

subroutine pi5076(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                  sigm, vim, dtau, etamin, etamax, copilo)

    use Behaviour_type
    use vmis_isot_nl_module, only: CONSTITUTIVE_LAW, Init, PathFollowing
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"
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
!  Path-following pred_elas for vmis_isot_nl
! --------------------------------------------------------------------------------------------------
! in  BEHInteg : constitutive law object
! in  typmod : type de modelisation
! in  ndim   : dimension de l'espace
! in  epsm   : deformations au temps moins
! in  epsd_cste   : correction de deformations dues aux charges fixes
! in  epsd_pilo   : correction de deformations dues aux charges pilotees
! in  sigm   : contraintes avec sqrt(2)
! in  vim    : variables internes en t-
! in  dtau   : 2nd membre de l'equation f(eta)=tau
! in  etamin : borne inf. pilotage
! in  etamax : borne sup. pilotage
! out copilo : coefficients de piotage
! --------------------------------------------------------------------------------------------------
    character(len=16), parameter:: option = 'PILO_PRED_ELAS'
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: ndimsi
    type(CONSTITUTIVE_LAW) :: ldc
! --------------------------------------------------------------------------------------------------

    ndimsi = BEHInteg%behavPara%ndimsi
    ASSERT(ndimsi .eq. 1)

    ! Constitutive law initialisation
    ldc = Init(ndimsi, option, 'NONE', &
               BEHInteg%materPara%schemePara%kpg, &
               BEHInteg%materPara%schemePara%ksp, &
               BEHInteg%materPara%jvMaterCode, &
               100, 1.d-6)

    ! Compute path-following coefficients
    call PathFollowing(ldc, dtau, vim, epsm+epsd_cste, epsd_pilo, copilo)
    if (ldc%exception .ne. 0) call utmess('F', 'PILOTAGE_83')

end subroutine
