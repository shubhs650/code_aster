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
subroutine getBehaviourAlgo(lPlaneStress, relaComp, &
                            relaCompPY, relaMecaPY, &
                            factorKeyword, iFactorKeyword, &
                            algo_inte, algo_inte_r)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterc/lcalgo.h"
#include "asterc/lctest.h"
#include "asterf_types.h"
#include "asterfort/getvtx.h"
#include "asterfort/utlcal.h"
#include "asterfort/utmess.h"
!
    aster_logical, intent(in) :: lPlaneStress
    character(len=16), intent(in) :: relaComp, relaCompPY, relaMecaPY
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=16), intent(out) :: algo_inte
    real(kind=8), intent(out) :: algo_inte_r
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get algorithm for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  lPlaneStress     : flag for plane stress model
! In  relaComp         : behaviour (RELATION keyword)
! In  relaCompPY       : behaviour (RELATION keyword) - For Python
! In  relaMecaPY       : mechanical part of behaviour - For Python
! In  factorKeyword    : factor keyword to read (COMPORTEMENT)
! In  iFactorKeyword   : index of factor keyword
! Out algo_inte        : algorithm for integration of behaviour
! Out algo_inte_r      : identifier for algorithm for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    algo_inte = ' '
    algo_inte_r = 0.d0

! - Get ALGO_INTE
    call getvtx(factorKeyword, 'ALGO_INTE', iocc=iFactorKeyword, scal=algo_inte, nbret=iret)
    if (iret .eq. 0) then
        call lcalgo(relaCompPY, algo_inte)
    else
        call lctest(relaMecaPY, 'ALGO_INTE', algo_inte, iret)
        if (iret .eq. 0) then
            call utmess('F', 'COMPOR1_45', nk=3, valk=[algo_inte, 'ALGO_INTE', relaComp])
        end if
    end if

! - Get ALGO_INTE - Plane stress
    if (lPlaneStress) then
        if (relaComp .eq. 'VMIS_ECMI_LINE' .or. relaComp .eq. 'VMIS_ECMI_TRAC' .or. &
            relaComp .eq. 'VMIS_ISOT_LINE' .or. relaComp .eq. 'VMIS_ISOT_TRAC') then
            algo_inte = 'SECANTE'
        end if
    end if

! - Convert name of algorithm to identifier
    call utlcal('NOM_VALE', algo_inte, algo_inte_r)
!
end subroutine
