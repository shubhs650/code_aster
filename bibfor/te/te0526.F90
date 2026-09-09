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
subroutine te0526(option, nomte)

    implicit none

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/elref2.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ngpide.h"
#include "asterfort/nmgvmb.h"
#include "asterfort/teattr.h"
!
    character(len=16), intent(in) :: option, nomte
! --------------------------------------------------------------------------------------------------
! Elementary computation
!
! Elements: 3D_GRAD_INCO, 3D_GRAD_VARI
!           D_PLAN_GRAD_INCO, D_PLAN_GRAD_VARI
!
! Options: PILO_PRED_DEFO
!
! --------------------------------------------------------------------------------------------------
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter:: ntrou_max = 10
! --------------------------------------------------------------------------------------------------
    aster_logical :: axi
    character(len=8) :: typmod(2), lielrf(ntrou_max)
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: ntrou, ndim, nnoL, nnoQ, npg, nddl, neps
    integer(kind=8) :: jv_poids, jv_vffL, jv_vffQ, jv_dfdeL, jv_dfdeQ, jv_geom
    integer(kind=8) :: jv_copil, jv_deplm, jv_ddepl, jv_depl0, jv_depl1, jv_dtau
    real(kind=8), allocatable:: b(:, :, :), w(:, :), ni2ldc(:, :)
! --------------------------------------------------------------------------------------------------
!
! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('C', 'TYPMOD2', typmod(2), vattr_missing=' ')
    axi = lteatt('AXIS', 'OUI')

! - Get parameters of element
    call elref2(nomte, ntrou_max, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(2), fami=fami, ndim=ndim, nno=nnoL, npg=npg, &
                     jpoids=jv_poids, jvf=jv_vffL, jdfde=jv_dfdeL)
    call elrefe_info(elrefe=lielrf(1), fami=fami, ndim=ndim, nno=nnoQ, npg=npg, &
                     jpoids=jv_poids, jvf=jv_vffQ, jdfde=jv_dfdeQ)

! - Option parameters
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PDEPLMR', 'L', jv_deplm)
    call jevech('PDDEPLR', 'L', jv_ddepl)
    call jevech('PDEPL0R', 'L', jv_depl0)
    call jevech('PDEPL1R', 'L', jv_depl1)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCDTAU', 'L', jv_dtau)
    call jevech('PCOPILO', 'E', jv_copil)

    ! Kinematics
    call nmgvmb(ndim, nnoQ, nnoL, npg, axi, &
                zr(jv_geom), zr(jv_vffQ), zr(jv_vffL), jv_dfdeQ, jv_dfdeL, &
                jv_poids, nddl, neps, b, w, ni2ldc)

    ! Computation of path-following coefficients
    call ngpide(compor, npg, neps, nddl, b, &
                zr(jv_deplm), zr(jv_ddepl), zr(jv_depl0), zr(jv_depl1), &
                zr(jv_dtau), zr(jv_copil), neps_meca=2*ndim)

    deallocate (b, w, ni2ldc)
end subroutine
