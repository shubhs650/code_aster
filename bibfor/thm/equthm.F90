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
! aslint: disable=W1504
!
subroutine equthm(ds_thm, option, &
                  lMatr, lSigm, &
                  lVari, lMatrPred, &
                  typmod, parm_theta, &
                  ndim, nbvari, &
                  kpi, npg, &
                  dimdef, dimcon, &
                  mecani, press1, press2, tempe, second, &
                  carcri, &
                  defgem, defgep, &
                  congem, congep, &
                  vintm, vintp, &
                  time_prev, time_curr, time_incr, &
                  r, drds, dsde, retcom)
!
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comthm.h"
#include "asterfort/thmComputeMatrix.h"
#include "asterfort/thmComputeResidual.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    character(len=16), intent(in) :: option
    aster_logical, intent(in) :: lMatr, lSigm, lVari, lMatrPred
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(in) :: parm_theta
    integer(kind=8), intent(in) :: ndim, nbvari
    integer(kind=8), intent(in) :: npg, kpi
    integer(kind=8), intent(in) :: dimdef, dimcon
    integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: defgem(dimdef), defgep(dimdef)
    real(kind=8), intent(inout) :: congem(dimcon), congep(dimcon)
    real(kind=8), intent(in) :: vintm(nbvari)
    real(kind=8), intent(inout) :: vintp(nbvari)
    real(kind=8), intent(in) :: time_prev, time_curr, time_incr
    real(kind=8), intent(out) :: r(dimdef+1)
    real(kind=8), intent(out) :: drds(dimdef+1, dimcon), dsde(dimcon, dimdef)
    integer(kind=8), intent(out) :: retcom
!
! --------------------------------------------------------------------------------------------------
!
! THM - Compute
!
! Compute generalized stresses and derivatives at current Gauss point - Unsteady version
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  option           : name of option- to compute
! In  typmod           : type of modelization (TYPMOD2)
! In  parm_theta       : parameter PARM_THETA
! In  ndim             : dimension of space (2 or 3)
! In  nbvari           : total number of internal state variables
! In  kpi              : current Gauss point
! In  npg              : number of Gauss points
! In  dimdef           : dimension of generalized strains vector
! In  dimcon           : dimension of generalized stresses vector
! In  mecani           : parameters for mechanic
! In  press1           : parameters for hydraulic (capillary pressure)
! In  press2           : parameters for hydraulic (gaz pressure)
! In  tempe            : parameters for thermic
! In  second           : parameters for second gradient
! In  defgem           : generalized strains - At begin of current step
! In  defgep           : generalized strains - At end of current step
! IO  congem           : generalized stresses - At begin of current step
! IO  congep           : generalized stresses - At end of current step
! In  vintm            : internal state variables - At begin of current step
! IO  vintp            : internal state variables - At end of current step
! In  time_prev        : time at beginning of step
! In  time_curr        : time at end of step
! In  time_incr        : time increment
! Out r                : non-linear residual vector
! Out drds             : derivative matrix (residual/stress)
! Out dsde             : derivative matrix (stress/strain (behaviour only)
! Out retcom           : return code for error
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i
    integer(kind=8) :: addeme, addete, addep1, addep2, adde2nd
    integer(kind=8) :: adcp11, adcp12, adcp21, adcp22
    integer(kind=8) :: adcome, adcote, adco2nd
    real(kind=8) :: gravity(3)
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
!
! --------------------------------------------------------------------------------------------------
!
    drds(1:dimdef+1, 1:dimcon) = 0.d0
    dsde(1:dimcon, 1:dimdef) = 0.d0
    r(1:dimdef+1) = 0.d0
    gravity = 0.d0
    retcom = 0

! - Address in generalized strains vector
    addeme = mecani(2)
    addete = tempe(2)
    addep1 = press1(3)
    addep2 = press2(3)
    adde2nd = second(2)

! - Address in generalized stresses vector
    adcome = mecani(3)
    adcote = tempe(3)
    adcp11 = press1(4)
    adcp12 = press1(5)
    adcp21 = press2(4)
    adcp22 = press2(5)
    adco2nd = second(3)

! - Add sqrt(2) for stresses
    if (ds_thm%ds_elem%l_dof_meca) then
        do i = 4, 6
            congem(adcome+i-1) = congem(adcome+i-1)*rac2
            congem(adcome+6+i-1) = congem(adcome+6+i-1)*rac2
        end do
    end if

! - Initialization of stresses
    if (lSigm) then
        do i = 1, dimcon
            congep(i) = congem(i)
        end do
    end if

! - Compute generalized stresses and derivatives at current Gauss point
    call comthm(ds_thm, &
                lMatr, lSigm, &
                lVari, lMatrPred, &
                option, typmod, &
                ndim, nbvari, &
                dimdef, dimcon, &
                adcome, adcote, adcp11, adcp12, adcp21, adcp22, adco2nd, &
                addeme, addete, addep1, addep2, adde2nd, &
                kpi, npg, &
                carcri, &
                defgem, defgep, &
                congem, congep, &
                vintm, vintp, &
                time_prev, time_curr, &
                dsde, gravity, retcom)
    if (retcom .ne. 0) then
        goto 99
    end if

! - Compute non-linear residual
    if (lSigm) then
        call thmComputeResidual(ds_thm, parm_theta, gravity, &
                                ndim, &
                                dimdef, dimcon, &
                                mecani, press1, press2, tempe, second, &
                                congem, congep, &
                                time_incr, &
                                r)
    end if

! - Compute jacobian
    if (lMatr) then
        call thmComputeMatrix(ds_thm, parm_theta, gravity, &
                              ndim, &
                              dimdef, dimcon, &
                              mecani, press1, press2, tempe, second, &
                              congem, congep, &
                              time_incr, &
                              drds)
    end if

! - Add sqrt(2) for stresses
    if (ds_thm%ds_elem%l_dof_meca .and. lSigm) then
        do i = 4, 6
            congep(adcome+i-1) = congep(adcome+i-1)/rac2
            congep(adcome+6+i-1) = congep(adcome+6+i-1)/rac2
        end do
    end if
!
99  continue
!
end subroutine
