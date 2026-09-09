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

subroutine te0131(nomopt, nomte)
    implicit none
!
#include "asterf_types.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/readVector.h"
#include "asterfort/writeVector.h"
#include "FE_module.h"
!
    character(len=16) :: nomte, nomopt
!
! ------------------------------------------------------
! -- HHO_DEPL_MECA: copy field only
!-------------------------------------------------------
    integer(kind=8) :: ndim, nno, nsize
    real(kind=8) :: field(MAX_BV_CG)
!
    call elrefe_info(fami="RIGI", ndim=ndim, nno=nno)
    nsize = nno*ndim
!
    call readVector("PDEPLPR", nsize, field)
    call writeVector("PDEPL_R", nsize, field)
!
end subroutine
