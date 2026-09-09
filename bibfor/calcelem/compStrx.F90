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
subroutine compStrx(modelZ, materCodeZ, caraElemZ, comporZ, &
                    dispZ, chgeomZ, &
                    chvarcZ, chvrefZ, &
                    lPoux, loadPres, coefMultR, &
                    ligrelZ, jvBaseZ, strxZ, codret)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/detrsd.h"
#include "asterfort/exisd.h"
#include "asterfort/exixfe.h"
#include "asterfort/jedetc.h"
#include "asterfort/jeexin.h"
#include "asterfort/meceuc.h"
#include "asterfort/mechpo.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: modelZ, materCodeZ, caraElemZ, comporZ
    character(len=*), intent(in) :: dispZ, chgeomZ
    character(len=*), intent(in) :: chvarcZ, chvrefZ
    aster_logical, intent(in) :: lPoux
    character(len=*), intent(in) :: loadPres
    real(kind=8), intent(in) :: coefMultR
    character(len=*), intent(in) :: ligrelZ, strxZ, jvBaseZ
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Fields computation
!
! Utility to compute STRX_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'STRX_ELGA'
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    character(len=1) :: jvBase
    character(len=8) :: model, caraElem
    character(len=24) :: chdisp, chelem, chdynr, suropt
    integer(kind=8) :: iret, nbFieldAdd
    aster_logical ::  lXFEM
    character(len=1), parameter :: coefType = "R"
    complex(kind=8), parameter :: coefMultC = (0.d0, 0.d0)
!
! --------------------------------------------------------------------------------------------------
!
    caraElem = caraElemZ
    chdisp = dispZ
    chelem = strxZ
    codret = 0
    jvBase = jvBaseZ
    model = modelZ
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '
    chdynr = ' '
    suropt = ' '

! - XFEM
    call exixfe(model, iret)
    lXFEM = iret .ne. 0
    if (lXFEM) then
        codret = 1
        call utmess('A', 'CALCCHAMP_7')
        goto 99
    end if

! - Add input fields
    lpain(1) = 'PDEPLAR'
    lchin(1) = chdisp
    lpain(2) = 'PCOMPOR'
    lchin(2) = comporZ
    lpain(3) = 'PGEOMER'
    lchin(3) = chgeomZ
    lpain(4) = 'PMATERC'
    lchin(4) = materCodeZ
    lpain(5) = 'PVARCRR'
    lchin(5) = chvrefZ
    lpain(6) = 'PVARCPR'
    lchin(6) = chvarcZ
    nbFieldIn = 6

! - Add field for beams
    if (lPoux) then
        call mechpo('&&MECHPO', loadPres, model, chdisp, chdynr, &
                    suropt, lpain(nbFieldIn+1), lchin(nbFieldIn+1), nbFieldAdd, &
                    coefType, coefMultR, coefMultC)
        nbFieldIn = nbFieldIn+nbFieldAdd
    end if

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Set output field
    lchout(1) = chelem
    lpaout(1) = 'PSTRXRR'

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
    if (lPoux) then
        call jedetc('V', '&&MECHPO', 1)
    end if
!
99  continue
!
end subroutine
