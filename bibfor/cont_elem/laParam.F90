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
subroutine laParam(parameters)
!
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "contact_module.h"
#include "jeveux.h"
!
    type(ContactParameters), intent(inout) :: parameters
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Get parameters from ContactComputaction.cxx/getData()
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jcont
!
    call jevech('PCONFR', 'L', jcont)
!
! - Reso Geom
!
    parameters%reso_geom = RESO_GEOM_NEWT
!
! - Contact
!
    parameters%algo_cont = nint(zr(jcont+23))
    parameters%type_cont = nint(zr(jcont+24))
    parameters%vari_cont = nint(zr(jcont+25))
    parameters%jac_type = nint(zr(jcont+26))
    parameters%inte_type = nint(zr(jcont+27))
!
    select case (parameters%vari_cont)
    case (CONT_VARI_NONE)
        parameters%vari_cont_coef = 0.d0
    case (CONT_VARI_RAPI)
        parameters%vari_cont_coef = 0.d0
    case (CONT_VARI_ROBU)
        parameters%vari_cont_coef = -1.d0
    case (CONT_VARI_CLAS)
        parameters%vari_cont_coef = -1.d0
    case (CONT_VARI_SYME)
        parameters%vari_cont_coef = 1.d0
    case default
        ASSERT(ASTER_FALSE)
    end select
!
! - Friction
!
    parameters%l_fric = (zr(jcont+30) > 0.5d0)
    parameters%algo_fric = nint(zr(jcont+31))
    parameters%type_fric = nint(zr(jcont+32))
    parameters%threshold_given = zr(jcont+34)
!
! - Other
!
    parameters%proj_tole = zr(jcont+40)
    parameters%cont_init = nint(zr(jcont+41))
    parameters%E = zr(jcont+45)
!
end subroutine
