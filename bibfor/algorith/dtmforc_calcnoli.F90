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

subroutine dtmforc_calcnoli(sd_dtm_, sd_nl_, buffdtm, buffnl, &
                            time, dt, depl, vite, fext, &
                            nbmode, nlacc, nlcase, idx_start, idx_end)
    implicit none

!
!
! dtmforc_choc : Calculates the stop/choc forces at the current step (t)
!
!       nl_ind           : nonlinearity index (for sd_nl access)
!       sd_dtm_, buffdtm : dtm data structure and its buffer
!       sd_nl_ , buffnl  : nl  data structure and its buffer
!       time             : current instant t
!       depl, vite       : structural modal displacement and velocity at "t"
!       fext_nl          : projected total non-linear force
!       fext_tgt         : projected tangential component of the fext_nl
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dtmforc_choc.h"
#include "asterfort/dtmforc_galet.h"
#include "asterfort/dtmforc_ants.h"
#include "asterfort/dtmforc_flam.h"
#include "asterfort/dtmforc_deci.h"
#include "asterfort/dtmforc_decr.h"
#include "asterfort/dtmforc_dvis.h"
#include "asterfort/dtmforc_rede.h"
#include "asterfort/dtmforc_revi.h"
#include "asterfort/dtmforc_rotf.h"
#include "asterfort/nlget.h"

!
!   -0.1- Input/output arguments
    integer(kind=8), intent(in)  :: nbmode, nlacc, nlcase, idx_start, idx_end
    character(len=*), intent(in)  :: sd_dtm_
    character(len=*), intent(in)  :: sd_nl_
    integer(kind=8), pointer  :: buffdtm(:)
    integer(kind=8), pointer  :: buffnl(:)
    real(kind=8), intent(in)  :: time
    real(kind=8), pointer  :: depl(:)
    real(kind=8), pointer  :: vite(:)
    real(kind=8), pointer :: fext(:)
    real(kind=8), intent(in)  :: dt
!
!   -0.2- Local variables
    integer(kind=8)           :: i, par_coorno(2)
    integer(kind=8)           ::  nl_ind, nl_type
    real(kind=8)      ::  eps
    character(len=8)  :: sd_dtm, sd_nl
!
    real(kind=8), pointer :: fext_nl(:) => null()
    real(kind=8), pointer :: fext_tgt(:) => null()
!
    data par_coorno/_COOR_NO1, _COOR_NO2/
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_nl = sd_nl_
    eps = r8prem()

    call nlget(sd_nl, _F_TOT_WK, vr=fext_nl, buffer=buffnl)
    call nlget(sd_nl, _F_TAN_WK, vr=fext_tgt, buffer=buffnl)

    do nl_ind = idx_start, idx_end
        call nlget(sd_nl, _NL_TYPE, iocc=nl_ind, iscal=nl_type, buffer=buffnl)
        select case (nl_type)
!
        case (NL_CHOC)
            call dtmforc_choc(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              time, depl, vite, fext_nl, fext_tgt)
!           Special case of choc nonlinearities that can be implicitely treated
            if ((nlcase .eq. 0) .or. (nlacc .eq. 1)) then
                do i = 1, nbmode
                    fext(i) = fext(i)+fext_nl(i)
                end do
            else
                do i = 1, nbmode
                    fext(i) = fext(i)+fext_tgt(i)
                end do
            end if

        case (NL_BUCKLING)
            call dtmforc_flam(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              time, depl, vite, fext)

        case (NL_ANTI_SISMIC)
            call dtmforc_ants(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              time, depl, vite, fext)

        case (NL_DIS_VISC)
            call dtmforc_dvis(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              time, dt, depl, vite, fext)

        case (NL_DIS_ECRO_TRAC)
            call dtmforc_decr(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              time, dt, depl, vite, fext)

        case (NL_DIS_CHOC_ELAS)
            call dtmforc_galet(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                               time, dt, depl, fext)

        case (NL_CRACKED_ROTOR)
            call dtmforc_rotf(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              time, depl, fext)
! !
        case (NL_FX_RELATIONSHIP)
            call dtmforc_rede(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              depl, fext)
! !
        case (NL_FV_RELATIONSHIP)
            call dtmforc_revi(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              vite, fext)
! !
        case (NL_DIS_ECRO_CINE)
            call dtmforc_deci(nl_ind, sd_dtm, sd_nl, buffdtm, buffnl, &
                              time, dt, depl, fext)

            !
        case default
            ASSERT(.false.)

        end select
    end do

end subroutine
