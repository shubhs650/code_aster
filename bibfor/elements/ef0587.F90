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
subroutine ef0587(nomte)
!
    use pipeElem_module
    implicit none
!
#include "asterfort/pipeElem_type.h"
#include "asterfort/tuefgeElno.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: TUYAU
!
! Option: EFGE_ELNO (non-linear case)
!
! --------------------------------------------------------------------------------------------------
!
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    aster_logical, parameter :: lLine = ASTER_FALSE
    integer(kind=8) :: nbFourier, nbDof, nbNode
!
! --------------------------------------------------------------------------------------------------
!
    call pipeGetDime(nomte, fami, &
                     nbNode, nbFourier, nbDof)

! - Compute option
    call tuefgeElno(lLine, nbNode, nbDof, nbFourier)
!
end subroutine
