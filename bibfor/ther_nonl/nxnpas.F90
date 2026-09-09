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
subroutine nxnpas(sddisc, solver, nume_inst, ds_print, &
                  lnkry, l_evol, l_stat, &
                  para, time_curr, deltat, reasma, &
                  tpsthe)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/diinst.h"
#include "asterfort/nmnkft.h"
#include "asterfort/SetTableColumn.h"
#include "asterfort/nonlinDSPrintInitTimeStep.h"
#include "asterfort/nonlinDSPrintHeadTimeStep.h"
!
    character(len=19), intent(in) :: sddisc, solver
    type(NL_DS_Print), intent(inout) :: ds_print
    integer(kind=8), intent(in) :: nume_inst
    aster_logical, intent(in) :: lnkry, l_evol, l_stat
    real(kind=8), intent(inout) :: para(2)
    real(kind=8), intent(out) :: time_curr, deltat
    aster_logical, intent(out) :: reasma
    real(kind=8), intent(out) :: tpsthe(6)
!
! --------------------------------------------------------------------------------------------------
!
! THER_NON_LINE - Algorithm
!
! Updates for new time step
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  solver           : datastructure for solver parameters
! IO  ds_print         : datastructure for printing parameters
! In  nume_inst        : index of current time step
! In  l_stat           : .true. is stationnary
! In  l_evol           : .true. if transient
! In  lnkry            : .true. if Newton-Krylov solver
! IO  para             : parameters for time
!                            (1) THETA
!                            (2) DELTA
! Out time_curr        : current time
! Out deltat           : increment of time
! Out reasma           : flag to assemble matrix
! Out tpsthe           : parameters for time
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: theta, theta_read, time_prev
!
! --------------------------------------------------------------------------------------------------
!
    theta_read = para(1)
!
! - RECUPERATION DU PAS DE TEMPS ET DES PARAMETRES DE RESOLUTION
!
    if (l_stat) then
        deltat = -1.d150
        theta = 1.d0
        if (l_evol) then
            time_curr = diinst(sddisc, nume_inst)
        else
            time_curr = 0.d0
        end if
    else
        time_curr = diinst(sddisc, nume_inst)
        time_prev = diinst(sddisc, nume_inst-1)
        deltat = time_curr-time_prev
        theta = theta_read
    end if
    para(2) = deltat
!
! - NEWTON-KRYLOV : COPIE DANS LA SD SOLVEUR DE LA PRECISION DE LA
!                     RESOLUTION POUR LA PREDICTION (FORCING-TERM)
    if (lnkry) then
        call nmnkft(solver, sddisc)
    end if
!
! - MATRICE TANGENTE REACTUALISEE POUR UN NOUVEAU DT
!
    reasma = ASTER_TRUE
!
! - Save values
!
    tpsthe(1) = time_curr
    tpsthe(2) = deltat
    tpsthe(3) = theta
    tpsthe(4) = r8vide()
    tpsthe(5) = r8vide()
    tpsthe(6) = r8vide()
!
! - Print management - Initializations for convergence table
!
    call SetTableColumn(ds_print%table_cvg, flag_acti_=ASTER_TRUE)
    call nonlinDSPrintInitTimeStep(ds_print)
!
! - Print management - Print head for new step time
!
    call nonlinDSPrintHeadTimeStep(sddisc, nume_inst, ds_print)
!
end subroutine
