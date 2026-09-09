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
subroutine thmEvalGravity(ds_thm, time, gravity)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvala.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: time
    real(kind=8), intent(out) :: gravity(3)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute gravity
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  time             : current time
! Out gravity          : gravity
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: grav_func(1)
    integer(kind=8), parameter :: nbProp = 3
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'PESA_X', 'PESA_Y', 'PESA_Z'/)
!
! --------------------------------------------------------------------------------------------------
!
    gravity = 0.d0
    propVale = 0.d0

! - Get parameters
    call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                ' ', 'THM_DIFFU', &
                0, ' ', [0.0d0], &
                nbProp, propName, propVale, &
                propCode, 0, nan='NON')

! - Get function
    call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                ' ', 'THM_DIFFU', &
                1, 'INST', [time], &
                1, 'PESA_MULT', grav_func, &
                propCode(1), 0, nan='NON')
    if (propCode(1) .eq. 1) then
        grav_func(1) = 1.d0
    end if
!
    gravity(1) = grav_func(1)*propVale(1)
    gravity(2) = grav_func(1)*propVale(2)
    gravity(3) = grav_func(1)*propVale(3)
!
end subroutine
