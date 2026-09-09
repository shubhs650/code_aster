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
subroutine pmdocr(carcriList)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/carc_chck.h"
#include "asterfort/carc_delete.h"
#include "asterfort/carc_info.h"
#include "asterfort/carc_read_pt.h"
#include "asterfort/setBehaviourParaValue.h"
!
    real(kind=8), intent(out) :: carcriList(CARCRI_SIZE)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics/SIMU_POINT_MAT)
!
! Get list of parameters for integration of constitutive law
!
! --------------------------------------------------------------------------------------------------
!
! Out carcriList       : list of parameters for integration of constitutive law
!
! --------------------------------------------------------------------------------------------------
!
    type(BehaviourPrep_MapCarcri) :: prepMapCarcri
!
! --------------------------------------------------------------------------------------------------
!
    carcriList(1:CARCRI_SIZE) = 0.d0

! - Create carcriList informations objects
    call carc_info(prepMapCarcri)

! - Read informations from command file
    call carc_read_pt(prepMapCarcri)

! - Some checks
    call carc_chck(prepMapCarcri)

! - Set in list
    call setBehaviourParaValue(prepMapCarcri%prepCrit, &
                               prepMapCarcri%parm_theta_thm, prepMapCarcri%parm_alpha_thm, &
                               carcriList_=carcriList(1:CARCRI_SIZE))

! - Cleaning
    call carc_delete(prepMapCarcri)
!
end subroutine
