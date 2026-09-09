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
subroutine thmCheckPorosity(relaMeca, ds_thm)
!
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: relaMeca
    type(THM_DS), intent(in) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Check porosity for some behaviours
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  relaMeca         : relation for mechanical part
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: poro_init, poro_meca, poro_diff, poro_tole
    integer(kind=8) :: propCode(1)
    real(kind=8) :: propVale(1)
!
! --------------------------------------------------------------------------------------------------
!
    poro_init = ds_thm%ds_parainit%poro_init

! - Check
    if (relaMeca .eq. 'CAM_CLAY') then
        poro_tole = 1.D-6
        call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                    ' ', relaMeca, &
                    0, ' ', [0.d0], &
                    1, ['PORO'], propVale, &
                    propCode, 0)
        poro_meca = propVale(1)
        poro_diff = abs(poro_meca-poro_init)
        if (abs(poro_diff) .gt. poro_tole) then
            call utmess('F', 'THM2_60', sk=relaMeca)
        end if
    end if
!
end subroutine
