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
subroutine comp_read_mfront(factorKeyword, iFactorKeyword, adrsMGIS)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/jeveuo.h"
!
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=16), intent(out) :: adrsMGIS
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get parameters for external programs (MFRONT)
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read (COMPORTEMENT)
! In  iFactorKeyword   : index of factor keyword
! Out adrsMGIS         : MGIS address
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: compMFront
    character(len=16), pointer :: compMFrontAddr(:) => null()
    integer(kind=8) :: nbret
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(iFactorKeyword .ne. 0)
    call getvid(factorKeyword, "COMPOR_MFRONT", iFactorKeyword, scal=compMFront, nbret=nbret)
    ASSERT(nbRet .eq. 1)
    call jeveuo(compMFront//'.ADDR', 'L', vk16=compMFrontAddr)
    adrsMGIS = compMFrontAddr(1)
!
end subroutine
