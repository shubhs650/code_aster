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
subroutine dtmforc_deci(nl_ind, sd_dtm_, sd_nl_, buffdtm, buffnl, &
                        time, step, depl, fext)
    implicit none
!
!
! dtmforc_decr : Calculates a "discrete model with isotropic behavior's"
!                force at the current step (t)
!
!       nl_ind           : nonlinearity index (for sd_nl access)
!       sd_dtm_, buffdtm : dtm data structure and its buffer
!       sd_nl_ , buffnl  : nl  data structure and its buffer
!       time, step       : current t and integration dt
!       depl             : structural modal displacement at "t"
!       fext             : projected total non-linear force
!
!
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dinon3.h"
#include "asterfort/dtmget.h"
#include "asterfort/fointe.h"
#include "asterfort/gloloc.h"
#include "asterfort/jeveuo.h"
#include "asterfort/locglo.h"
#include "asterfort/nlget.h"
#include "asterfort/tophys.h"
#include "asterfort/tophys_ms.h"
#include "asterfort/togene.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
!
!   -0.1- Input/output arguments
    integer(kind=8), intent(in) :: nl_ind
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_nl_
    integer(kind=8), pointer :: buffdtm(:)
    integer(kind=8), pointer :: buffnl(:)
    real(kind=8), intent(in) :: time
    real(kind=8), intent(in) :: step
    real(kind=8), pointer :: depl(:)
    real(kind=8), pointer :: fext(:)
!
!   -0.2- Local variables
    aster_logical :: multi_support, l_rota, okdire(6)
    integer(kind=8) :: iex, nbexci, ier, nbno, nbddl, iddl, nbmode
    integer(kind=8) :: ino, start, finish
    real(kind=8) :: sina, cosa, sinb, cosb, sing
    real(kind=8) :: cosg, depglo(6), deploc(12)
    real(kind=8) :: dflocal(6), dfglob(6), dul(12), deplocm(12)
    real(kind=8) :: coeflo(6, 4), raide(6), fglob(6)
    real(kind=8) :: floc_lin(6), fglob_lin(6), vdeploc(6), fglob_correc(6)
    real(kind=8) :: varmoi(18), varplu(18), origob(3)
    character(len=8) :: sd_dtm, sd_nl, monmot, obst_typ
    character(len=19) :: nomres
!
    integer(kind=8), pointer :: vindx(:) => null()
    real(kind=8), pointer :: coedep(:) => null()
    real(kind=8), pointer :: coevit(:) => null()
    real(kind=8), pointer :: psidel(:) => null()
    real(kind=8), pointer :: psidel1(:) => null()
    real(kind=8), pointer :: psidel2(:) => null()
    real(kind=8), pointer :: sincos_angle_a(:) => null()
    real(kind=8), pointer :: sincos_angle_b(:) => null()
    real(kind=8), pointer :: sincos_angle_g(:) => null()
    real(kind=8), pointer :: dplmod(:) => null()
    real(kind=8), pointer :: dplmod1(:) => null()
    real(kind=8), pointer :: dplmod2(:) => null()
    real(kind=8), pointer :: vint(:) => null()
    real(kind=8), pointer :: raide0(:) => null()
    real(kind=8), pointer :: limy(:) => null()
    real(kind=8), pointer :: kcine(:) => null()
    real(kind=8), pointer :: puis(:) => null()
    real(kind=8), pointer :: limu(:) => null()
    character(len=8), pointer :: nofdep(:) => null()
    character(len=8), pointer :: nofvit(:) => null()
!
! --------------------------------------------------------------------------------------------------
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_nl = sd_nl_
!
    call nlget(sd_nl, _INTERNAL_VARS, vr=vint, buffer=buffnl)
    call nlget(sd_nl, _INTERNAL_VARS_INDEX, vi=vindx, buffer=buffnl)
    start = vindx(nl_ind)

!   ddl de rotation ou non
    call nlget(sd_nl, _ECRCIN_KELA, iocc=nl_ind, vr=raide0, buffer=buffnl)
    raide(:) = raide0(:)

    l_rota = ASTER_FALSE
    nbddl = 3
    if (sum(raide(4:6)) .gt. 0.d0) then
        l_rota = ASTER_TRUE
        nbddl = 6
    end if
!
    deploc(:) = 0.d0
!
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode, buffer=buffdtm)
    call dtmget(sd_dtm, _MULTI_AP, kscal=monmot, buffer=buffdtm)
    multi_support = monmot(1:3) .eq. 'OUI'
    if (multi_support) then
        call dtmget(sd_dtm, _CALC_SD, kscal=nomres, buffer=buffdtm)
        call dtmget(sd_dtm, _NB_EXC_T, iscal=nbexci, buffer=buffdtm)
!
        call jeveuo(nomres//'.FDEP', 'L', vk8=nofdep)
        call jeveuo(nomres//'.FVIT', 'L', vk8=nofvit)
!
        AS_ALLOCATE(vr=coedep, size=nbexci)
        AS_ALLOCATE(vr=coevit, size=nbexci)
        do iex = 1, nbexci
            coedep(iex) = 0.d0
            coevit(iex) = 0.d0
            if (nofdep(iex) .ne. ' ') then
                call fointe('F', nofdep(iex), 1, ['INST'], [time], &
                            coedep(iex), ier)
            end if
            if (nofvit(iex) .ne. ' ') then
                call fointe('F', nofvit(iex), 1, ['INST'], [time], &
                            coevit(iex), ier)
            end if
        end do
    end if
!
!
    nbno = 1
    call nlget(sd_nl, _OBST_TYP, iocc=nl_ind, kscal=obst_typ, buffer=buffnl)
    call nlget(sd_nl, _MODAL_DEPL_NO1, iocc=nl_ind, vr=dplmod1, buffer=buffnl)
    if (multi_support) call nlget(sd_nl, _PSI_DELT_NO1, vr=psidel1, buffer=buffnl)
!
    if (obst_typ(1:2) .eq. 'BI') then
        nbno = 2
        call nlget(sd_nl, _MODAL_DEPL_NO2, iocc=nl_ind, vr=dplmod2, buffer=buffnl)
        if (multi_support) call nlget(sd_nl, _PSI_DELT_NO2, vr=psidel2, buffer=buffnl)

        call nlget(sd_nl, _SINCOS_ANGLE_A, iocc=nl_ind, vr=sincos_angle_a, buffer=buffnl)
        call nlget(sd_nl, _SINCOS_ANGLE_B, iocc=nl_ind, vr=sincos_angle_b, buffer=buffnl)
        call nlget(sd_nl, _SINCOS_ANGLE_G, iocc=nl_ind, vr=sincos_angle_g, buffer=buffnl)
        sina = sincos_angle_a(1)
        cosa = sincos_angle_a(2)
        sinb = sincos_angle_b(1)
        cosb = sincos_angle_b(2)
        sing = sincos_angle_g(1)
        cosg = sincos_angle_g(2)
    end if
!
    origob(:) = 0.d0
    do ino = 1, nbno
! Point toward the modal displacement for the concerned node / 1 or 2 /
        dplmod => dplmod1
        if (multi_support) psidel => psidel1
        if (ino .eq. 2) then
            dplmod => dplmod2
            if (multi_support) psidel => psidel2
        end if
! Conversion of generalized displacements/velocities
! back to the physical (global) basis
        if (multi_support) then
            call tophys_ms(dplmod, psidel, coedep, depl, depglo, nbddl)
        else
            call tophys(dplmod, depl, depglo, nbddl)
        end if
! Conversion of these vectors to the local basis if needed
        if (nbno .eq. 2) then
            call gloloc(depglo(1:3), origob, sina, cosa, sinb, &
                        cosb, sing, cosg, deploc(1+(ino-1)*nbddl))
            if (l_rota) then
                call gloloc(depglo(4:6), origob, sina, cosa, sinb, &
                            cosb, sing, cosg, deploc(4+(ino-1)*nbddl))

            end if
        else
            deploc(1:nbddl) = depglo(1:nbddl)
        end if
    end do
    vdeploc(:) = 0.d0
    if (nbno .eq. 1) then
        vdeploc(1:nbddl) = deploc(1:nbddl)
    else
        vdeploc(1:nbddl) = deploc(nbddl+1:2*nbddl)-deploc(1:nbddl)
    end if
!

!   At initialization, the step is set to zero
    if (abs(step) .le. r8prem()) then
        deplocm(:) = 0.d0
        varmoi(:) = 0.d0
        fglob(:) = 0.d0
    else
        deplocm(:) = vint(start-1+1:start-1+12)
        varmoi(:) = vint(start-1+13:start-1+30)
        fglob(:) = vint(start-1+31:start-1+36)
    end if
! Physical (behavior) parameters
    okdire(:) = ASTER_FALSE
    coeflo(:, :) = 0.d0

    call nlget(sd_nl, _ECRCIN_LIMY, iocc=nl_ind, vr=limy, buffer=buffnl)
    call nlget(sd_nl, _ECRCIN_KCIN, iocc=nl_ind, vr=kcine, buffer=buffnl)
    call nlget(sd_nl, _ECRCIN_PUIS, iocc=nl_ind, vr=puis, buffer=buffnl)
    call nlget(sd_nl, _ECRCIN_LIMU, iocc=nl_ind, vr=limu, buffer=buffnl)

!   ordre des paramètres attendu dans dinon3 :
!   'LIMU_DX', 'PUIS_DX', 'KCIN_DX', 'LIMY_DX' ...
    do iddl = 1, nbddl
        if (raide(iddl) .gt. 0.d0) then
            okdire(iddl) = ASTER_TRUE
            coeflo(iddl, 3) = kcine(iddl)
            coeflo(iddl, 4) = limy(iddl)
            coeflo(iddl, 2) = puis(iddl)
            coeflo(iddl, 1) = limu(iddl)
        end if
    end do

    dul(:) = deploc(:)-deplocm(:)

    call dinon3(nbno*nbddl, deplocm, dul, deploc, nbno, &
                nbddl, varmoi, raide, 4, coeflo, &
                okdire, varplu, dflocal)

!   Conversion to the global (physical) reference
    if (nbno .eq. 2) then
        call locglo(dflocal(1:3), sina, cosa, sinb, cosb, &
                    sing, cosg, dfglob(1:3))
        if (l_rota) call locglo(dflocal(4:6), sina, cosa, sinb, cosb, &
                                sing, cosg, dfglob(4:6))
    else
        dfglob(1:nbddl) = dflocal(1:nbddl)
    end if
    fglob(:) = fglob(:)+dfglob(:)
!
!   Force linéaire correspondante
    floc_lin(:) = raide0(:)*vdeploc(:)
    if (nbno .eq. 2) then
        call locglo(floc_lin(1:3), sina, cosa, sinb, cosb, &
                    sing, cosg, fglob_lin(1:3))
        if (l_rota) call locglo(floc_lin(4:6), sina, cosa, sinb, cosb, &
                                sing, cosg, fglob_lin(4:6))
    else
        fglob_lin(1:nbddl) = floc_lin(1:nbddl)
    end if
!
!   Force corrective par rapport au linéaire
    fglob_correc(:) = fglob(:)-fglob_lin(:)

!   Generalized force on nodes
    if (nbno .eq. 2) then
        call togene(dplmod1, fglob_correc, fext, nbddl_=nbddl)
        call togene(dplmod2, fglob_correc, fext, coef=-1.d0, nbddl_=nbddl)
    else
        call togene(dplmod1, fglob_correc, fext, coef=-1.d0, nbddl_=nbddl)
    end if
!
!   Internal variables, storage
    if (step .gt. 0.d0) then
        finish = vindx(nl_ind+1)
        ASSERT((finish-start) .eq. NBVARINT_DECI)
! Local displacement of node 1
        vint(start) = deploc(1)
        vint(start+1) = deploc(2)
        vint(start+2) = deploc(3)
        vint(start+3) = deploc(4)
        vint(start+4) = deploc(5)
        vint(start+5) = deploc(6)
! Local displacement of node 2
        vint(start+6) = deploc(7)
        vint(start+7) = deploc(8)
        vint(start+8) = deploc(9)
        vint(start+9) = deploc(10)
        vint(start+10) = deploc(11)
        vint(start+11) = deploc(12)
! Internal variables of the behavior
        vint(start-1+13:start-1+30) = varplu(:)
! Global forces
        vint(start-1+31:start-1+36) = fglob(:)
    end if
!
    if (multi_support) then
        AS_DEALLOCATE(vr=coedep)
        AS_DEALLOCATE(vr=coevit)
    end if
end subroutine
