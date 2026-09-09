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
subroutine te0589(option, nomte)
!
    use pipeElem_module
!
    implicit none
!
#include "asterfort/pipeElem_type.h"
#include "asterfort/tutemp.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: TUYAU_*
!
! Option: CHAR_MECA_TEMP_R
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbNode, nbFourier, nbDof
!
! --------------------------------------------------------------------------------------------------
!
    call pipeGetDime(nomte, 'RIGI', &
                     nbNode, nbFourier, nbDof)

! - Compute options
    call tutemp(nbNode, nbDof, nbFourier)
!
end subroutine
