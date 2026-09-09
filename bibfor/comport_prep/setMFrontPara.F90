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
! aslint: disable=W0413
!
subroutine setMFrontPara(prepCrit, iFactorKeyword)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterc/mgis_get_double_mfront_parameter.h"
#include "asterc/mgis_get_integer_mfront_parameter.h"
#include "asterc/mgis_set_double_parameter.h"
#include "asterc/mgis_set_integer_parameter.h"
#include "asterc/mgis_set_outofbounds_policy.h"
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/utmess.h"
!
    type(BehaviourPrep_Crit), pointer :: prepCrit(:)
    integer(kind=8), intent(in) :: iFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Set parameters for MFront
!
! --------------------------------------------------------------------------------------------------
!
! Ptr prepCrit         : pointer to behaviour criteria
! In  iFactorKeyword   : index of factor keyword (for map)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: resi_inte_mfront_rela, valr(2)
    integer(kind=8):: solvBehavType, iveriborne
    character(len=16) :: adrsMGIS
    real(kind=8) :: resi_inte
    integer(kind=8) :: iter_inte_maxi, iter_inte_mfront_maxi, vali(2)
!
! --------------------------------------------------------------------------------------------------
!
    iveriborne = prepCrit(iFactorKeyword)%iveriborne
    adrsMGIS = prepCrit(iFactorKeyword)%prepExte%adrsMGIS
    solvBehavType = prepCrit(iFactorKeyword)%prepExte%solvBehavType

! - Set values
    if (solvBehavType == SOLV_BEHAV_MGIS_OFFI .or. &
        solvBehavType == SOLV_BEHAV_MGIS_PROTO) then
        if (associated(prepCrit(iFactorKeyword)%resi_inte)) then
            resi_inte = prepCrit(iFactorKeyword)%resi_inte
            call mgis_get_double_mfront_parameter(adrsMGIS, "epsilon", resi_inte_mfront_rela)
            if (resi_inte_mfront_rela .ne. 0.d0 .and. &
                resi_inte .gt. resi_inte_mfront_rela) then
                valr(1) = resi_inte
                valr(2) = resi_inte_mfront_rela
                call utmess('A', 'COMPOR6_16', nr=2, valr=valr)
            end if
            call mgis_set_double_parameter(adrsMGIS, "epsilon", resi_inte)
        end if
        if (associated(prepCrit(iFactorKeyword)%iter_inte_maxi)) then
            iter_inte_maxi = prepCrit(iFactorKeyword)%iter_inte_maxi
            call mgis_get_integer_mfront_parameter(adrsMGIS, "iterMax", iter_inte_mfront_maxi)
            if (iter_inte_mfront_maxi .ne. 0 .and. &
                iter_inte_maxi .ne. iter_inte_mfront_maxi) then
                vali(1) = iter_inte_maxi
                vali(2) = iter_inte_mfront_maxi
                call utmess('I', 'COMPOR6_17', ni=2, vali=vali)
            end if
            call mgis_set_integer_parameter(adrsMGIS, "iterMax", iter_inte_maxi)
        end if
        call mgis_set_outofbounds_policy(adrsMGIS, iveriborne)
    end if
!
end subroutine
