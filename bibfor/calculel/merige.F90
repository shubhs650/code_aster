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
subroutine merige(modelZ, caraElemZ, sigm, strx, materElem, &
                  jvBase, numeHarm, disp_, materCode_)
!
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
#include "asterfort/setStructFields.h"
#include "asterfort/xajcin.h"
!
    integer(kind=8) :: numeHarm
    character(len=1) :: jvBase
    character(len=*) :: sigm, strx
    character(len=19) :: materElem
    character(len=*), intent(in) :: modelZ
    character(len=*), intent(in) :: caraElemZ
    character(len=*), optional, intent(in) :: disp_
    character(len=*), optional, intent(in) :: materCode_
!
! --------------------------------------------------------------------------------------------------
!
! Elementary matrix for RIGI_GEOM
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    character(len=16), parameter :: option = 'RIGI_GEOM'
    character(len=24) :: modelLigrel, chgeom, chharm
    character(len=8) :: model, caraElem
    integer(kind=8) :: icode, ier
    aster_logical :: lXFEM

! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    model = modelZ
    caraElem = caraElemZ
    if (model(1:1) .eq. ' ') then
        call utmess('F', 'CALCULEL2_82')
    end if
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '

! - Prepare flags
    call exixfe(model, ier)
    lXFEM = ier .ne. 0

! - Preparation of input fields
    call mecham(option, model, numeHarm, &
                chgeom, chharm, icode)

! - Prepare MATER_ELEM
    call detrsd('MATR_ELEM', materElem)
    call memare(jvBase, materElem, model, option)

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PHARMON'
    lchin(2) = chharm(1:19)
    lpain(3) = 'PCONTRR'
    lchin(3) = sigm
    lpain(4) = 'PSTRXRR'
    lchin(4) = strx
    lpain(5) = 'PEFFORR'
    lchin(5) = sigm
    nbFieldIn = 5
    if (present(disp_)) then
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PDEPLPR'
        lchin(nbFieldIn) = disp_

    end if
    if (present(materCode_)) then
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PMATERC'
        lchin(nbFieldIn) = materCode_
    end if

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add input XFEM fields if required
    if (lxfem) then
        call xajcin(model, option, nbFieldInMax, lchin, lpain, nbFieldIn)
    end if

! - Add output field
    lpaout(1) = 'PMATUUR'
    lchout(1) = materElem(1:15)//'.ME001'

    call calcul('S', option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')
    call reajre(materElem, lchout(1), jvBase)
!
    call jedema()
end subroutine
