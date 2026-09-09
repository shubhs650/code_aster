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
subroutine compStress(modelZ, modelLigrelZ, &
                      materCodeZ, caraElemZ, comporZ, &
                      chdispZ, chgeom, &
                      chtime, chharm, &
                      chvarc, chvref, &
                      jvBaseZ, chelemZ, codret)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cesvar.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/exisd.h"
#include "asterfort/exixfe.h"
#include "asterfort/jeexin.h"
#include "asterfort/meceuc.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
#include "asterfort/xajcin.h"
!
    character(len=*), intent(in) :: modelZ, modelLigrelZ
    character(len=*), intent(in) :: materCodeZ, caraElemZ, comporZ
    character(len=*), intent(in) :: chdispZ, chgeom
    character(len=*), intent(in) :: chtime, chharm
    character(len=*), intent(in) :: chvarc, chvref
    character(len=*), intent(in) :: chelemZ, jvBaseZ
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Fields computation
!
! Utility to compute SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'SIEF_ELGA'
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    character(len=1) :: jvBase
    character(len=8) :: model, caraElem

    character(len=19), parameter :: canbsp = '&&MECALC.NBSP'
    character(len=24) :: chdisp, chelem
    integer(kind=8) :: iret
    aster_logical :: lXFEM
!
! --------------------------------------------------------------------------------------------------
!
    chdisp = chdispZ
    chelem = chelemZ
    jvBase = jvBaseZ
    model = modelZ
    caraElem = caraElemZ
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '
    codret = 0
    call exixfe(model, iret)
    lXFEM = iret .ne. 0

! - Add output field
    lchout(1) = chelem
    lpaout(1) = 'PCONTRR'

! - Preparation for dynamic fields ('sub-points')
    call exisd('CHAM_ELEM_S', canbsp, iret)
    if (iret .ne. 1) then
        call cesvar(caraElem, ' ', modelLigrelZ, canbsp)
    end if
    call copisd('CHAM_ELEM_S', 'V', canbsp, chelem)

! - Add input fields
    lpain(1) = 'PDEPLAR'
    lchin(1) = chdisp
    lpain(2) = 'PABSCUR'
    lchin(2) = chgeom(1:8)//'.ABSC_CURV'
    lpain(3) = 'PCOMPOR'
    lchin(3) = comporZ
    lpain(4) = 'PGEOMER'
    lchin(4) = chgeom
    lpain(5) = 'PHARMON'
    lchin(5) = chharm
    lpain(6) = 'PMATERC'
    lchin(6) = materCodeZ
    lpain(7) = 'PVARCRR'
    lchin(7) = chvref
    lpain(8) = 'PVARCPR'
    lchin(8) = chvarc
    lpain(9) = 'PINSTR'
    lchin(9) = chtime
    nbFieldIn = 9

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add XFEM fields
    if (lXFEM) then
        call xajcin(model, 'REFE_FORC_NODA', nbFieldInMax, lchin, lpain, &
                    nbFieldIn)
    end if

! - Computation (with preparation for COMPLEX fields)
    call meceuc('C', option, caraElem, modelLigrelZ, &
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
