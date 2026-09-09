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
subroutine te0562(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrefv.h"
#include "asterfort/jevech.h"
#include "asterfort/massup.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_GRAD_VARI, 3D_GRAD_INCO
!           D_PLAN_GRAD_VARI, D_PLAN_GRAD_INCO
!           AXIS_GRAD_VARI, AXIS_GRAD_INCO
!
! Options: MASS_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "MASS"
    integer(kind=8) :: nb_DOF
    integer(kind=8) :: nnoQ, npg, imatuu, ndim, nnos, jv_ganoQ, icodr1(1)
    integer(kind=8) :: jv_poids, jv_vfQ, jv_dfdeQ, jvGeom, jvMaterc
    integer(kind=8) :: nnoL, jv_vfL, jv_dfdeL, jv_ganoL
!
! --------------------------------------------------------------------------------------------------
!

! - Get parameters of element
    call elrefv(fami, ndim, &
                nnoL, nnoQ, nnos, &
                npg, jv_poids, &
                jv_vfL, jv_vfQ, &
                jv_dfdeL, jv_dfdeQ, &
                jv_ganoL, jv_ganoQ)
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)

! - Input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PMATUUR', 'E', imatuu)

! - nb_DOF: displacements (2 or 3) + LAMBDA + VAR_REG
    nb_DOF = ndim+2

! - Compute mass matrix
    call massup(zi(jvMaterc), &
                option, ndim, nb_DOF, nnoQ, nnoL, &
                npg, jv_poids, jv_dfdeQ, &
                zr(jvGeom), zr(jv_vfQ), imatuu, icodr1, jvGeom, &
                jv_vfQ)
!
end subroutine
