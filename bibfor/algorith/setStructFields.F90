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
subroutine setStructFields(caraElemZ, nbFieldInMax, lchin, lpain, nbFieldIn)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/mecara.h"
!
    character(len=*), intent(in) :: caraElemZ
    integer(kind=8), intent(in) :: nbFieldInMax
    character(len=*), intent(inout) :: lpain(nbFieldInMax)
    character(len=*), intent(inout) :: lchin(nbFieldInMax)
    integer(kind=8), intent(inout) :: nbFieldIn
!
! --------------------------------------------------------------------------------------------------
!
! Add fields for structural elements
!
! --------------------------------------------------------------------------------------------------
!
! In  caraElem         : name of elementary characteristics (field)
! In  nbFieldInMax     : maximum number of input fields
! IO  lpain            : list of input parameters
! IO  lchin            : list of input fields
! IO  nbFieldIn        : number of input fields
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbFieldAdd
    character(len=24) :: chcara(18)
!
! --------------------------------------------------------------------------------------------------
!
    call mecara(caraElemZ, chcara)
    nbFieldAdd = 14
    ASSERT(nbFieldIn+nbFieldAdd .le. nbFieldInMax)
    lpain(nbFieldIn+1) = 'PCADISK'
    lchin(nbFieldIn+1) = chcara(2)
    lpain(nbFieldIn+2) = 'PCADISM'
    lchin(nbFieldIn+2) = chcara(3)
    lpain(nbFieldIn+3) = 'PCADISA'
    lchin(nbFieldIn+3) = chcara(4)
    lpain(nbFieldIn+4) = 'PCAGEPO'
    lchin(nbFieldIn+4) = chcara(5)
    lpain(nbFieldIn+5) = 'PCAGNPO'
    lchin(nbFieldIn+5) = chcara(6)
    lpain(nbFieldIn+6) = 'PCASECT'
    lchin(nbFieldIn+6) = chcara(8)
    lpain(nbFieldIn+7) = 'PCAARPO'
    lchin(nbFieldIn+7) = chcara(9)
    lpain(nbFieldIn+8) = 'PCACABL'
    lchin(nbFieldIn+8) = chcara(10)
    lpain(nbFieldIn+9) = 'PCAGNBA'
    lchin(nbFieldIn+9) = chcara(11)
    lpain(nbFieldIn+10) = 'PCAPOUF'
    lchin(nbFieldIn+10) = chcara(13)
    lpain(nbFieldIn+11) = 'PVENTCX'
    lchin(nbFieldIn+11) = chcara(14)
    lpain(nbFieldIn+12) = 'PCINFDI'
    lchin(nbFieldIn+12) = chcara(15)
    lpain(nbFieldIn+13) = 'PNBSP_I'
    lchin(nbFieldIn+13) = chcara(16)
    lpain(nbFieldIn+14) = 'PFIBRES'
    lchin(nbFieldIn+14) = chcara(17)
    nbFieldIn = nbFieldIn+nbFieldAdd
!
end subroutine
