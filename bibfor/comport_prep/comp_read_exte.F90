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
subroutine comp_read_exte(factorKeyword, iFactorKeyword, &
                          librNameUMAT, subrNameUMAT, nbVariUMAT)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
!
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=255), intent(out) :: librNameUMAT
    character(len=255), intent(out) :: subrNameUMAT
    integer(kind=8), intent(out) :: nbVariUMAT
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get parameters for external programs (UMAT)
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read (COMPORTEMENT)
! In  iFactorKeyword   : index of factor keyword
! Out librNameUMAT     : name of library
! Out subrNameUMAT     : name of behaviour in library
! Out nbVariUMAT       : number of internal variables
!
! --------------------------------------------------------------------------------------------------
!
    librNameUMAT = ' '
    subrNameUMAT = ' '
    nbVariUMAT = 0
    ASSERT(iFactorKeyword .ne. 0)
    call getvtx(factorKeyword, 'LIBRAIRIE', iocc=iFactorKeyword, scal=librNameUMAT)
    call getvtx(factorKeyword, 'NOM_ROUTINE', iocc=iFactorKeyword, scal=subrNameUMAT)
    call getvis(factorKeyword, 'NB_VARI', iocc=iFactorKeyword, scal=nbVariUMAT)
!
end subroutine
