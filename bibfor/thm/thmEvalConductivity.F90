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
subroutine thmEvalConductivity(ds_thm, &
                               ndim, &
                               satur, phi, &
                               lambs, dlambs, lambp, dlambp, &
                               tlambt, tlamct, tdlamt)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterfort/rcvala.h"
#include "asterfort/tdlamb.h"
#include "asterfort/telamb.h"
#include "asterfort/THM_type.h"
#include "asterfort/tlambc.h"
!
    type(THM_DS), intent(in) :: ds_thm
    integer(kind=8), intent(in) :: ndim
    real(kind=8), intent(in) :: satur, phi
    real(kind=8), intent(out) :: lambs, dlambs
    real(kind=8), intent(out) :: lambp, dlambp
    real(kind=8), intent(out) :: tlambt(ndim, ndim)
    real(kind=8), intent(out) :: tlamct(ndim, ndim)
    real(kind=8), intent(out) :: tdlamt(ndim, ndim)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Evaluate thermal conductivity
!n
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  ndim             : dimension of space
! In  satur            : saturation
! In  phi              : porosity
! Out lambs            : thermal conductivity depending on saturation
! Out dlambs           : derivative of thermal conductivity depending on saturation
! Out lambp            : thermal conductivity depending on porosity
! Out dlambp           : derivative of thermal conductivity depending on porosity
! Out tlambt           : tensor of thermal conductivity
! Out tlamct           : tensor of thermal conductivity (constant part)
! Out tdlamt           : tensor of dnerivatives for thermal conductivity
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPara = 2
    real(kind=8) :: paraVale(nbPara)
    character(len=4), parameter :: paraName(nbPara) = (/'SAT ', 'PORO'/)
    integer(kind=8), parameter :: nbProp = 4
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'LAMB_S  ', 'D_LB_S  ', &
                                                         'LAMB_PHI', 'D_LB_PHI'/)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: anglNaut(3)
!
! --------------------------------------------------------------------------------------------------
!
    anglNaut = ds_thm%ds_behaviour%BEHInteg%materPara%lcsPara%lcsAngle
    lambs = 1.d0
    dlambs = 0.d0
    lambp = 1.d0
    dlambp = 0.d0
    propVale = 0.d0
    propVale(1) = 1.d0
    propVale(3) = 1.d0

! - Get parameters depending on porosity and saturation
    paraVale(1) = satur
    paraVale(2) = phi
    call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                ' ', 'THM_DIFFU', &
                nbPara, paraName, paraVale, &
                nbProp, propName, propVale, &
                propCode, 0, nan='NON')
    lambs = propVale(1)
    dlambs = propVale(2)
    lambp = propVale(3)
    dlambp = propVale(4)

! - Compute tensor of thermal conductivity
    call telamb(ds_thm, anglNaut, ndim, tlambt)

! - Compute tensor of thermal conductivity (constant part)
    call tlambc(ds_thm, anglNaut, ndim, tlamct)

! - Compute tensor of derivatives (by temperature) for thermal conductivity
    call tdlamb(ds_thm, anglNaut, ndim, tdlamt)
!
end subroutine
