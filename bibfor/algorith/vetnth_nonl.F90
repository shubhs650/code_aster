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

subroutine vetnth_nonl(model, caraElem, materCode, time, comporTher, &
                       tempIter, hydrPrev, &
                       varcPrev, varcCurr, &
                       jvBase, vectElemLine, vectElemNLin)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/megeom.h"
#include "asterfort/reajre.h"
#include "asterfort/vemare.h"
!
    character(len=8), intent(in) :: model, caraElem
    character(len=24), intent(in) :: materCode, time, comporTher
    character(len=24), intent(in) :: tempIter, hydrPrev
    character(len=19), intent(in) :: varcPrev, varcCurr
    character(len=1), intent(in) :: jvBase
    character(len=24), intent(in) :: vectElemLine, vectElemNLin
!
! --------------------------------------------------------------------------------------------------
!
! Thermic - Residuals
!
! Evolution for non-linear (CHAR_THER_EVOLNI)
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  caraElem         : name of elementary characteristics (field)
! In  time             : time (<CARTE>)
! In  comporTher       : name of <CARTE> COMPOR
! In  tempIter         : temperature field at current Newton iteration
! In  hydrPrev         : previous hydration
! In  varcCurr         : command variable for current time
! In  varcPrev         : command variable for previous time
! In  vectElemLine     : name of vect_elem result (linear part)
! In  vectElemNLin     : name of vect_elem result (non linear part)
! In  jvBase           : JEVEUX jvBase for object
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
    character(len=16), parameter :: option = 'CHAR_THER_EVOLNI'
    integer(kind=8) :: iret, nbFieldIn
    character(len=8) :: newnom

    character(len=24) :: modelLigrel
    character(len=19) :: resuElemLine, resuElemNLin
    character(len=24) :: chgeom
!
! --------------------------------------------------------------------------------------------------
!
    newnom = '.0000000'
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Prepare VECT_ELEM
    call jeexin(vectElemLine(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        call vemare(jvBase, vectElemLine, model)
    else
        call jedetr(vectElemLine(1:19)//'.RELR')
    end if
    call jeexin(vectElemNLin(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        call vemare(jvBase, vectElemNLin, model)
    else
        call jedetr(vectElemNLin(1:19)//'.RELR')
    end if

! - Get geometry field
    call megeom(model, chgeom)

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PTEMPER'
    lchin(2) = tempIter(1:19)
    lpain(3) = 'PMATERC'
    lchin(3) = materCode(1:19)
    lpain(4) = 'PINSTR'
    lchin(4) = time(1:19)
    lpain(5) = 'PVARCPR'
    lchin(5) = varcCurr(1:19)
    lpain(6) = 'PHYDRPM'
    lchin(6) = hydrPrev(1:19)
    lpain(7) = 'PCOMPOR'
    lchin(7) = comporTher(1:19)
    lpain(8) = 'PVARCMR'
    lchin(8) = varcPrev(1:19)
    nbFieldIn = 8

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Generate new RESU_ELEM name
    resuElemNLin = vectElemNLin(1:8)//'.0000000'
    newnom = resuElemNLin(10:16)
    call gcnco2(newnom)
    resuElemNLin(10:16) = newnom(2:8)

! - Output fields
    lpaout(1) = 'PVECTTI'
    lchout(1) = resuElemNLin
    call corich('E', lchout(1), ichin_=-1)

! - Generate new RESU_ELEM name
    resuElemLine = vectElemLine(1:8)//'.0000000'
    newnom = resuElemLine(10:16)
    call gcnco2(newnom)
    resuElemLine(10:16) = newnom(2:8)

! - Set output fields
    lpaout(2) = 'PVECTTR'
    lchout(2) = resuElemLine
    call corich('E', lchout(2), ichin_=-1)

! - Compute
    call calcul('S', option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')

! - Add RESU_ELEM in VECT_ELEM
    call reajre(vectElemNLin, lchout(1), jvBase)
    call reajre(vectElemLine, lchout(2), jvBase)
!
end subroutine
