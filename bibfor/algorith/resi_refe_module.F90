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

module resi_refe_module

    implicit none
    private
    public:: Init, GetRef, Check, RESI_REFE

#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"

    type RESI_REFE
        character(len=16), private:: nomte = ''
        integer(kind=8), private:: maxi_refe = 0
        character(len=8), pointer, private:: names(:) => null()
        real(kind=8), pointer, private:: values(:) => null()
        aster_logical, allocatable, private:: requested_values(:)
    contains
        procedure, pass:: Init
        procedure, pass:: GetRef
        procedure, pass:: Check
    end type RESI_REFE

contains

    subroutine Init(self, nomte)
        implicit none
        class(RESI_REFE), intent(out):: self
        character(len=*), intent(in):: nomte
        ! ------------------------------------------------------------------------------------------
        self%nomte = nomte
        call jevech('PRESICMP', 'L', vk8=self%names)
        call jevech('PRESIREF', 'L', vr=self%values)
        self%maxi_refe = size(self%names)
        ASSERT(self%maxi_refe .ge. 1)
        ASSERT(size(self%values) .eq. self%maxi_refe)

        allocate (self%requested_values(self%maxi_refe))
        self%requested_values = ASTER_FALSE
    end subroutine Init

    function GetRef(self, comp_name_) result(vale_refe)
        implicit none
        class(RESI_REFE), intent(inout):: self
        character(len=*), intent(in):: comp_name_
        real(kind=8) :: vale_refe
        ! ------------------------------------------------------------------------------------------
        integer(kind=8), pointer:: pos(:) => null()
        integer(kind=8):: idx
        character(len=8):: comp_name
        character(len=16) :: kmess(2)
        ! ------------------------------------------------------------------------------------------
        ASSERT(self%maxi_refe .ne. 0)
        comp_name = comp_name_

        pos = findloc(self%names, comp_name)
        ASSERT(size(pos) .eq. 1)
        idx = pos(1)
        ASSERT(idx .ge. 1 .and. idx .le. self%maxi_refe)

        vale_refe = self%values(idx)
        if (vale_refe .eq. r8vide()) then
            kmess(1) = self%nomte
            kmess(2) = comp_name
            call utmess('F', 'MECANONLINE5_55', nk=2, valk=kmess)
        end if

        self%requested_values(idx) = ASTER_TRUE
    end function GetRef

    subroutine Check(self)
        implicit none
        class(RESI_REFE), intent(in):: self
        ! ------------------------------------------------------------------------------------------
        integer(kind=8):: i
        character(len=16) :: kmess(2)
        ! ------------------------------------------------------------------------------------------
        do i = 1, self%maxi_refe
            if (.not. self%requested_values(i) .and. self%values(i) .ne. r8vide()) then
                kmess(1) = self%nomte
                kmess(2) = self%names(i)
                call utmess('A', 'MECANONLINE5_57', nk=2, valk=kmess)
            end if
        end do
    end subroutine Check
end module
