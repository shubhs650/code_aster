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
subroutine nmvcfo(poum, model, materField, materCode, caraElem, compor, &
                  varcRefe, hval_incr, vectElem)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmvarc_prep.h"
#include "asterfort/nmvccc.h"
#include "asterfort/nmvcd2.h"
#include "asterfort/reajre.h"
#include "asterfort/vemare.h"
!
    character(len=1), intent(in) :: poum
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: materField, materCode, caraElem
    character(len=24), intent(in) :: varcRefe
    character(len=24), intent(in) :: compor
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: vectElem
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Command variables - Vector for reference (residual evaluation)
!
! --------------------------------------------------------------------------------------------------
!
! In  poum      : type of computation
!                      '-' - Previous step
!                      '+' - Current step
! In  model          : name of model
! In  materField     : name of material characteristics (field)
! In  materCode      : name of coded material
! In  caraElem       : name of elementary characteristics (field)
! In  varcRefe       : name of reference command variables vector
! In  compor         : name of comportment definition (field)
! In  hval_incr      : hat-variable for incremental values
! In  vectElem      : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
    character(len=8) :: lpaout(nbFieldOutMax), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOutMax), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn, nbFieldout
!
    aster_logical :: exis_temp, exis_hydr, exis_ptot, exis_sech, exis_epsa
    aster_logical :: exis_meta_zirc, exis_meta_acier, exis_meta, calc_meta
    character(len=19) :: sigmPrev, variPrev, varcPrev, varcCurr
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!

! - Get fields from hat-variables - Begin of time step
    call nmchex(hval_incr, 'VALINC', 'SIGMOI', sigmPrev)
    call nmchex(hval_incr, 'VALINC', 'VARMOI', variPrev)
    call nmchex(hval_incr, 'VALINC', 'COMMOI', varcPrev)
    call nmchex(hval_incr, 'VALINC', 'COMPLU', varcCurr)

! - Get state of external state variables
    call nmvcd2('HYDR', materField, exis_hydr)
    call nmvcd2('PTOT', materField, exis_ptot)
    call nmvcd2('SECH', materField, exis_sech)
    call nmvcd2('EPSA', materField, exis_epsa)
    call nmvcd2('M_ZIRC', materField, exis_meta_zirc)
    call nmvcd2('M_ACIER', materField, exis_meta_acier)
    call nmvcd2('TEMP', materField, exis_temp)
    exis_meta = exis_temp .and. (exis_meta_zirc .or. exis_meta_acier)
    calc_meta = ASTER_FALSE
    if (exis_meta .and. poum .eq. '+') then
        calc_meta = ASTER_TRUE
    end if

! - Prepare elementary vectors
    call jeexin(vectElem(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        call vemare('V', vectElem, model)
    end if
    call jedetr(vectElem(1:19)//'.RELR')
    call reajre(vectElem, ' ', 'V')

! - Fields preparation of elementary vectors
    call nmvarc_prep(poum, model, caraElem, materCode, varcRefe, &
                     compor, exis_temp, &
                     nbFieldInMax, nbFieldIn, lpain, lchin, &
                     nbFieldOutMax, nbFieldout, lpaout, lchout, &
                     sigmPrev, variPrev, varcPrev, varcCurr)

! - Computation of elementary vectors
    call nmvccc(model, &
                nbFieldInMax, nbFieldIn, lpain, lchin, &
                nbFieldOutMax, nbFieldout, lpaout, lchout, &
                exis_temp, exis_hydr, exis_ptot, &
                exis_sech, exis_epsa, calc_meta, &
                vectElem)
!
end subroutine
