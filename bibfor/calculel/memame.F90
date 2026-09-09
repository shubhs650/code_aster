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
subroutine memame(optionz, modelz, materFieldZ, materCodeZ, caraElemz, time, &
                  comporMultz, matrElemz, jvBaseZ, listElemCalcz)
!
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/redetr.h"
#include "asterfort/setStructFields.h"
#include "asterfort/vrcins.h"
#include "asterfort/xajcin.h"
!
    character(len=*), intent(in) :: optionz
    character(len=*), intent(in) :: modelz, materFieldZ, materCodeZ, caraElemz
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: comporMultz, matrElemz
    character(len=*), intent(in) :: jvBaseZ, listElemCalcz
!
! --------------------------------------------------------------------------------------------------
!
! Elementary matrix MASS_*
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option
! In  model            : name of the model
! In  materField       : name of material characteristics (field)
! In  materCode        : name of coded material
! In  caraElem         : name of elementary characteristics (field)
! In  time             : current time
! In  comporMult       : name of comportment definition for PMF (field)
! In  jvBase             : JEVEUX jvBase to create matrElem
! In  matrElem         : elementary matrix
! In  listElemCalc     : list of elements (LIGREL) where matrElem is computed
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    integer(kind=8) :: nbFieldIn, nbFieldOut
    character(len=2) :: codret
    integer(kind=8) :: iret
    integer(kind=8), parameter :: numeHarm = 0
    character(len=16) :: option
    character(len=24), parameter :: chvarc = '&&MERIME.CHVARC'
    character(len=24) :: comporMult, listElemCalc
    character(len=24) :: chgeom, chharm
    character(len=1) :: jvBase
    character(len=8) :: model, caraElem
    character(len=24) :: materField, materCode
    character(len=19) :: matrElem
    integer(kind=8) :: nbSubstruct
    aster_logical :: lxfem, hasFiniteElement
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    option = optionz
    model = modelz
    caraElem = caraElemz
    materField = materFieldZ
    materCode = materCodeZ
    matrElem = matrElemz
    comporMult = comporMultz
    jvBase = jvBaseZ
    listElemCalc = listElemCalcz
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '

! - Prepare flags
    call exixfe(model, iret)
    lxfem = iret .ne. 0
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nbSubstruct)

! - Preparation of input fields
    call mecham(option, model, numeHarm, &
                chgeom, chharm, iret)
    hasFiniteElement = iret .eq. 0

! - Field for external state variables
    call vrcins(model, materField, caraElem, time, chvarc, codret)

! - Prepare RESU_ELEM objects
    call jeexin(matrElem(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        call memare(jvBase, matrElem, model, option, to_aster_logical(nbSubstruct > 0))
    else
        call jedetr(matrElem(1:19)//'.RELR')
    end if

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = materCode(1:19)
    lpain(3) = 'PABSCUR'
    lchin(3) = chgeom(1:8)//'.ABSC_CURV'
    lpain(4) = 'PVARCPR'
    lchin(4) = chvarc(1:19)
    lpain(5) = 'PCOMPOR'
    lchin(5) = comporMult(1:19)
    nbFieldIn = 5

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add input XFEM fields if required
    if (lxfem) then
        call xajcin(model, option, nbFieldInMax, lchin, lpain, nbFieldIn)
    end if
!
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)
!
! - Output fields
    lpaout(1) = 'PMATUUR'
    lchout(1) = matrElem(1:15)//'.M01'
    lpaout(2) = 'PMATUNS'
    lchout(2) = matrElem(1:15)//'.M02'
    if (option .eq. 'MASS_MECA') then
        nbFieldOut = 2
    else
        nbFieldOut = 1
    end if

! - Mass
    if (hasFiniteElement) then
! ----- Compute
        call calcul('S', &
                    option, listElemCalc, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')

! ----- Save RESU_ELEM
        call reajre(matrElem, lchout(1), jvBase)
        if (nbFieldOut .eq. 2) then
            call reajre(matrElem, lchout(2), jvBase)
        end if

    end if

! - Clean
    call redetr(matrElem)
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
