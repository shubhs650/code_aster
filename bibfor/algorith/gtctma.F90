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

subroutine gtctma(elem_coor, elem_nbnode, elem_code, elem_dime, &
                  ctcoor)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/reerel.h"
#include "MeshTypes_type.h"
!
!

    real(kind=8), intent(in) :: elem_coor(3, MT_NNOMAX2D)
    integer(kind=8), intent(in) :: elem_nbnode
    character(len=8), intent(in) :: elem_code
    integer(kind=8), intent(in) :: elem_dime
    real(kind=8), intent(out) :: ctcoor(3)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Get center of a given contact element
!
! --------------------------------------------------------------------------------------------------
!
! In elem_coor        : coordinates of nodes for current element
! In elem_nbnode      : number of node for current element
! In elem_code        : code of current element
! In elem_dime        : dimension of current element
! Out ctcoor           : coordonate of the center
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: rfcoor(3)
!
! --------------------------------------------------------------------------------------------------
!
!
! - Initialisation
!
    ctcoor(1:3) = 0.d0
    select case (elem_code)
    case ('SE2', 'SE3', 'SE4')
        rfcoor(1) = 0.d0
        rfcoor(2) = 0.d0
        rfcoor(3) = 0.d0
    case ('TR3', 'TR6', 'TR7', 'TR1')
        rfcoor(1) = 1.d0/3.d0
        rfcoor(2) = 1.d0/3.d0
        rfcoor(3) = 0.d0
    case ('QU4', 'QU8', 'QU9', 'Q12')
        rfcoor(1) = 0.d0
        rfcoor(2) = 0.d0
        rfcoor(3) = 0.d0
    case default
        ASSERT(.false.)
    end select
!
! - Compute center
!
    call reerel(elem_code, elem_nbnode, 3, elem_coor, rfcoor, ctcoor)
!
! - Print check
!
    !write(*,*)ctcoor(:)
!
end subroutine
