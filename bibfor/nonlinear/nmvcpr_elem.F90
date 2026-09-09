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
subroutine nmvcpr_elem(modelZ, materFieldZ, materCodeZ, caraElemZ, &
                       poum, hval_incr, &
                       varcRefeZ, comporZ, &
                       vectElemZ)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/detrsd.h"
#include "asterfort/maveElemCreate.h"
#include "asterfort/nmchex.h"
#include "asterfort/varcCalcComp.h"
#include "asterfort/varcCalcMeta.h"
#include "asterfort/varcCalcPrep.h"
#include "asterfort/varcDetect.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, materCodeZ
    character(len=1), intent(in) :: poum
    character(len=*), intent(in) :: varcRefeZ, comporZ
    character(len=19), intent(in) :: hval_incr(*)
    character(len=*), intent(in) :: vectElemZ
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Command variables - Elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  model          : name of model
! In  materField     : name of material characteristics (field)
! In  materCode      : name of coded material
! In  caraElem       : name of elementary characteristics (field)
! In  poum           :  '-' or '+' for command variables evaluation
! In  hval_incr      : hat-variable for incremental values
! In  varcRefe       : name of reference command variables vector
! In  compor         : name of comportment definition (field)
! In  vectElem      : elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: jvBase = "V"
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
    aster_logical :: l_temp, l_hydr, l_ptot
    aster_logical :: l_sech, l_epsa, l_meta
    character(len=19) :: sigmPrev, variPrev, varcPrev, varcCurr
    integer(kind=8) :: nbFieldIn, nbFieldOut
    character(len=24) :: multComp
    character(len=24), parameter :: chsith = '&&VECTME.CHSITH'
!
! --------------------------------------------------------------------------------------------------
!
    multComp = comporZ

! - Get fields from hat-variables - Begin of time step
    call nmchex(hval_incr, 'VALINC', 'SIGMOI', sigmPrev)
    call nmchex(hval_incr, 'VALINC', 'VARMOI', variPrev)
    call nmchex(hval_incr, 'VALINC', 'COMMOI', varcPrev)
    call nmchex(hval_incr, 'VALINC', 'COMPLU', varcCurr)

! - Detect external state variables
    call varcDetect(materFieldZ, l_temp, l_hydr, l_ptot, l_sech, l_epsa, l_meta)

! - Prepare elementary vectors
    call detrsd('VECT_ELEM', vectElemZ)
    call maveElemCreate(jvBase, vectElemZ, modelZ)

! - Preparation
    call varcCalcPrep(modelZ, caraElemZ, materCodeZ, &
                      poum, &
                      l_temp, l_meta, &
                      varcRefeZ, varcPrev, varcCurr, &
                      comporZ, multComp, chsith, &
                      sigmPrev, variPrev, &
                      nbFieldInMax, nbFieldOutMax, &
                      nbFieldIn, nbFieldOut, &
                      lpain, lchin, &
                      lpaout, lchout)

! - Calls to CALCUL
    call varcCalcComp(modelZ, chsith, &
                      l_temp, l_hydr, l_ptot, &
                      l_sech, l_epsa, &
                      nbFieldIn, nbFieldOut, &
                      lpain, lchin, &
                      lpaout, lchout, &
                      vectElemZ)

! - Call to CALCUL special for metallurgy (non-incremental)
    if (poum .eq. '+' .and. l_meta) then
        call varcCalcMeta(modelZ, &
                          nbFieldIn, nbFieldOut, &
                          lpain, lchin, &
                          lpaout, lchout, &
                          vectElemZ)
    end if
!
end subroutine
