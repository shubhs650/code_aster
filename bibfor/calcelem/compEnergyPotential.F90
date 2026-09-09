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
subroutine compEnergyPotential(optionZ, &
                               modelZ, materCodeZ, caraElemZ, comporZ, &
                               chdispZ, chharmZ, chgeomZ, &
                               chtimeZ, chvarcZ, chvrefZ, &
                               l_temp, chtempZ, &
                               ligrelZ, jvBaseZ, epotElemZ, &
                               codret)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/detrsd.h"
#include "asterfort/exisd.h"
#include "asterfort/meceuc.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: optionZ
    character(len=*), intent(in) :: modelZ, materCodeZ, caraElemZ, comporZ
    character(len=*), intent(in) :: chdispZ, chharmZ, chgeomZ, chtimeZ
    character(len=*), intent(in) :: chvarcZ, chvrefZ
    aster_logical, intent(in) :: l_temp
    character(len=*), intent(in) :: chtempZ
    character(len=*), intent(in) :: ligrelZ, jvBaseZ, epotElemZ
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Fields computation
!
! Utility to compute ETHE_ELEM and EPOT_ELEM
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    character(len=1) :: jvBase
    character(len=8) :: model, caraElem
    character(len=24) :: chdisp, chelem, chtemp
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    chdisp = chdispZ
    chelem = epotElemZ
    chtemp = chtempZ
    jvBase = jvBaseZ
    model = modelZ
    caraElem = caraElemZ
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '
    codret = 0

! - Add input fields
    if (l_temp) then
        lpain(1) = 'PTEMPER'
        lchin(1) = chtemp
    else
        lpain(1) = 'PDEPLAR'
        lchin(1) = chdisp
    end if
    lpain(2) = 'PCOMPOR'
    lchin(2) = comporZ
    lpain(3) = 'PGEOMER'
    lchin(3) = chgeomZ
    lpain(4) = 'PHARMON'
    lchin(4) = chharmZ
    lpain(5) = 'PMATERC'
    lchin(5) = materCodeZ
    lpain(6) = 'PVARCRR'
    lchin(6) = chvrefZ
    lpain(7) = 'PVARCPR'
    lchin(7) = chvarcZ
    lpain(8) = 'PINSTR'
    lchin(8) = chtimeZ
    nbFieldIn = 8

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Set output field
    lchout(1) = chelem
    lpaout(1) = 'PENERDR'

! - Computation (with preparation for COMPLEX fields)
    call meceuc('C', optionZ, caraElem, ligrelZ, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, jvBase)
    call exisd('CHAMP_GD', lchout(1), iret)
    if (iret .eq. 0) then
        codret = 1
        call utmess('A', 'CALCCHAMP_89', sk=optionZ)
    end if

! - Clean
    call detrsd('CHAM_ELEM_S', chelem)
!
end subroutine
