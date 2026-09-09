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

subroutine nmvarc_prep(poum, model, caraElem, materCode, varcRefe, &
                       compor, exis_temp, &
                       nbFieldInMax, nbFieldIn, lpain, lchin, &
                       nbFieldOutMax, nbFieldOut, lpaout, lchout, &
                       sigmPrev, variPrev, varcPrev, varcCurr)
!
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/alchml.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/nmvcex.h"
#include "asterfort/setStructFields.h"
#include "asterfort/xajcin.h"
!
    character(len=1), intent(in) :: poum
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: materCode
    character(len=24), intent(in) :: varcRefe
    character(len=24), intent(in) :: caraElem
    character(len=24), intent(in) :: compor
    aster_logical, intent(in) :: exis_temp
    integer(kind=8), intent(in) :: nbFieldInMax
    character(len=8), intent(inout) :: lpain(nbFieldInMax)
    character(len=19), intent(inout) :: lchin(nbFieldInMax)
    integer(kind=8), intent(out) :: nbFieldIn
    integer(kind=8), intent(in) :: nbFieldOutMax
    character(len=8), intent(inout) :: lpaout(nbFieldOutMax)
    character(len=19), intent(inout) :: lchout(nbFieldOutMax)
    integer(kind=8), intent(out) :: nbFieldOut
    character(len=19), intent(in) :: sigmPrev, variPrev
    character(len=19), intent(in) :: varcPrev, varcCurr
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Command variables - Fields preparation
!
! --------------------------------------------------------------------------------------------------
!
! In  poum           : type of computation
!                      '-' - Previous step
!                      '+' - Current step
! In  model          : name of model
! In  materCode      : name of coded material
! In  caraElem       : name of elementary characteristics (field)
! In  varcRefe       : name of reference command variables vector
! In  compor         : name of comportment definition (field)
! In  exis_temp      : .true. if temperature variable command exists
! In  nbFieldInMax         : maximum number of input fields
! IO  lpain          : list of input parameters
! IO  lchin          : list of input fields
! IO  nbFieldIn           : number of input fields
! In  nbFieldOutMax        : maximum number of output fields
! IO  lpaout         : list of output parameters
! IO  lchout         : list of output fields
! IO  nbFieldOut          : number of output fields
! In  sigmPrev      : stress at previous step
! In  variPrev      : internal variables at previous step
! In  varcPrev      : command variables at previous step
! In  varcCurr      : command variables at current step
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: chsith = '&&NMVCPR.CHSITH'
    integer(kind=8), parameter :: numeHarm = 0
    aster_logical :: lxfem
    integer(kind=8) :: iret
    character(len=19) :: varcAllPrev, varcAllCurr, timeCurr, timePrev
    character(len=24) :: chgeom, modelLigrel, chharm, varcAllRefe
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "
    call exixfe(model, iret)
    lxfem = (iret .eq. 1)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! - Get external state variables
    call nmvcex('TOUT', varcRefe, varcAllRefe)
    call nmvcex('TOUT', varcCurr, varcAllCurr)
    call nmvcex('TOUT', varcPrev, varcAllPrev)

! - Get field for time
    call nmvcex('INST', varcPrev, timePrev)
    call nmvcex('INST', varcCurr, timeCurr)

! - Get geometry field
    call megeom(model, chgeom)

! - Create field for Fourier mode
    call meharm(model, numeHarm, chharm)

! - Add input fields
    lpain(1) = 'PVARCRR'
    lchin(1) = varcAllRefe(1:19)
    lpain(2) = 'PGEOMER'
    lchin(2) = chgeom(1:19)
    lpain(3) = 'PMATERC'
    lchin(3) = materCode(1:19)
    lpain(4) = 'PCONTMR'
    lchin(4) = sigmPrev
    lpain(5) = 'PVARIPR'
    lchin(5) = variPrev
    lpain(6) = 'PCOMPOR'
    lchin(6) = compor(1:19)
    lpain(7) = 'PHARMON'
    lchin(7) = chharm(1:19)
    nbFieldIn = 7

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add specific input fields
    if (poum .eq. '-') then
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PINSTR'
        lchin(nbFieldIn) = timePrev
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVARCPR'
        lchin(nbFieldIn) = varcAllPrev

    else if (poum .eq. '+') then
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PINSTR'
        lchin(nbFieldIn) = timeCurr
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVARCPR'
        lchin(nbFieldIn) = varcAllCurr
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVARCMR'
        lchin(nbFieldIn) = varcAllPrev

    end if

! - Add XFEM input fields
    if (lxfem .and. exis_temp) then
        call xajcin(model, 'CHAR_MECA_TEMP_R', nbFieldInMax, lchin, lpain, nbFieldIn)
    end if

! - Add HHO fields
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Set output fields
    lpaout(1) = 'PVECTUR'
    nbFieldOut = 1

! - Add XFEM output fields
    if (lxfem .and. exis_temp) then
        call detrsd('CHAM_ELEM', chsith)
        call alchml(modelLigrel, 'SIEF_ELGA', 'PCONTRR', 'V', chsith, iret, ' ')
        lpaout(2) = 'PCONTRT'
        lchout(2) = chsith(1:19)
        nbFieldOut = nbFieldOut+1
    end if
!
    call jedema()
end subroutine
