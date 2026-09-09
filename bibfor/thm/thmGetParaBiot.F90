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
subroutine thmGetParaBiot(ds_thm)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8nnem.h"
#include "asterfort/rcvala.h"
#include "asterfort/THM_type.h"
!
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get Biot resumeters (for porosity evolution) (THM_DIFFU)
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 4
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'BIOT_COEF', 'BIOT_L   ', &
                                                         'BIOT_N   ', 'BIOT_T   '/)
    real(kind=8) :: emmag, phi0
    real(kind=8), parameter :: eps = 1.d-21
!
! --------------------------------------------------------------------------------------------------
!
    propVale(:) = r8nnem()

    call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                ' ', 'THM_DIFFU', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 0, nan='OUI')
    ds_thm%ds_material%biot%coef = propVale(1)
    ds_thm%ds_material%biot%l = propVale(2)
    ds_thm%ds_material%biot%n = propVale(3)
    ds_thm%ds_material%biot%t = propVale(4)

! - Type
    if (propCode(1) .eq. 0) then
        ds_thm%ds_material%biot%type = BIOT_TYPE_ISOT
    else
        if (propCode(4) .eq. 0) then
            ds_thm%ds_material%biot%type = BIOT_TYPE_ORTH
        else
            ds_thm%ds_material%biot%type = BIOT_TYPE_ISTR
        end if
    end if

! - If small storage coefficient
    if (ds_thm%ds_material%hydr%l_emmag) then
        emmag = ds_thm%ds_material%hydr%emmag
        phi0 = ds_thm%ds_parainit%poro_init
        if (emmag .lt. eps) then
            ds_thm%ds_material%biot%coef = phi0
            ds_thm%ds_material%biot%l = phi0
            ds_thm%ds_material%biot%t = phi0
            ds_thm%ds_material%biot%t = phi0
        end if
    end if
!
end subroutine
