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
subroutine compEnergyKinetic(lModal, &
                             modelZ, materCodeZ, caraElemZ, &
                             dispZ, viteZ, chFreqZ, chgeomZ, &
                             chmasdZ, chvarcZ, &
                             ligrelZ, jvBaseZ, ecinElemZ, &
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
    aster_logical, intent(in) :: lModal
    character(len=*), intent(in) :: modelZ, materCodeZ, caraElemZ
    character(len=*), intent(in) :: dispZ, viteZ, chgeomZ, chmasdZ
    character(len=*), intent(in) :: chFreqZ, chvarcZ
    character(len=*), intent(in) :: ligrelZ, jvBaseZ, ecinElemZ
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Fields computation
!
! Utility to compute ECIN_ELEM
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'ECIN_ELEM'
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    character(len=1) :: jvBase
    character(len=8) :: model, caraElem
    character(len=24) :: chdisp, chvite, chfreq, chelem
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    caraElem = caraElemZ
    chdisp = dispZ
    chvite = viteZ
    chfreq = chFreqZ
    chelem = ecinElemZ
    codret = 0
    jvBase = jvBaseZ
    model = modelZ
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '

! - Add input fields
    lpain(1) = 'POMEGA2'
    lchin(1) = chFreqZ
    lpain(2) = 'PMASDIA'
    lchin(2) = chmasdZ
    lpain(3) = 'PGEOMER'
    lchin(3) = chgeomZ
    lpain(4) = 'PMATERC'
    lchin(4) = materCodeZ
    lpain(5) = 'PVARCPR'
    lchin(5) = chvarcZ
    lpain(6) = 'PDEPLAR'
    if (lModal) then
        lchin(6) = chdisp
    else
        lchin(6) = chvite
    end if
    nbFieldIn = 6

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Set output field
    lchout(1) = chelem
    lpaout(1) = 'PENERCR'

! - Computation (with preparation for COMPLEX fields)
    call meceuc('C', option, caraElem, ligrelZ, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, jvBase)
    call exisd('CHAMP_GD', lchout(1), iret)
    if (iret .eq. 0) then
        codret = 1
        call utmess('A', 'CALCCHAMP_89', sk=option)
    end if

! - Clean
    call detrsd('CHAM_ELEM_S', chelem)
!
end subroutine
