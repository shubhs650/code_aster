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
subroutine getMFrontPlaneStress(relaComp, relaCompPY, &
                                factorKeyword, iFactorKeyword, &
                                l_mfront_cp)
!
    implicit none
!
#include "asterc/lctest.h"
#include "asterf_types.h"
#include "asterfort/getvtx.h"
!
    character(len=16), intent(in) :: relaComp, relaCompPY
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    aster_logical, intent(out) :: l_mfront_cp
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get type of plane stress hypothesis for MFront
!
! --------------------------------------------------------------------------------------------------
!
! In  relaComp         : behaviour (RELATION keyword)
! In  relaCompPY       : behaviour (RELATION keyword) - For Python
! In  factorKeyword    : factor keyword to read (COMPORTEMENT)
! In  iFactorKeyword   : index of factor keyword
! Out l_mfront_cp      : .true. if analytical plane stress
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=16) :: answer
!
! --------------------------------------------------------------------------------------------------
!
    l_mfront_cp = ASTER_FALSE
!
    if (relaComp .eq. 'MFRONT') then
! ----- Prototype
        call getvtx(factorKeyword, 'ALGO_CPLAN', iocc=iFactorKeyword, scal=answer)
        l_mfront_cp = answer .eq. 'ANALYTIQUE'
    else
! ----- Official
        call lctest(relaCompPY, 'MODELISATION', 'C_PLAN', iret)
        l_mfront_cp = iret .ne. 0
    end if
!
end subroutine
