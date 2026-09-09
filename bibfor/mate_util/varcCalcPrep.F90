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
! aslint: disable=W1504
!
subroutine varcCalcPrep(modelZ, caraElemZ, materCodeZ, &
                        poum, &
                        l_temp, l_meta, &
                        varcRefeZ, varcPrevZ, varcCurrZ, &
                        comporZ, multCompZ, chsithz, &
                        sigmz, variz, &
                        nbFieldInMax, nbFieldOutMax, &
                        nbFieldIn, nbFieldOut, &
                        lpain, lchin, &
                        lpaout, lchout)
!
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/nmvcex.h"
#include "asterfort/setStructFields.h"
#include "asterfort/xajcin.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ, materCodeZ
    aster_logical, intent(in) :: l_temp, l_meta
    character(len=1), intent(in) :: poum
    character(len=*), intent(in) :: varcRefeZ, varcPrevZ, varcCurrZ
    character(len=*), intent(in) :: comporZ, multCompZ, chsithz
    character(len=*), intent(in) :: sigmz, variz
    integer(kind=8), intent(in) :: nbFieldInMax, nbFieldOutMax
    integer(kind=8), intent(out) :: nbFieldIn, nbFieldOut
    character(len=8), intent(out) :: lpaout(nbFieldOutMax), lpain(nbFieldInMax)
    character(len=19), intent(out) :: lchout(nbFieldOutMax), lchin(nbFieldInMax)
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! Preparation to compute elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  caraElem         : name of elementary characteristics (field)
! In  materCode        : name of coded material
! In  poum             :  '-' or '+' for command variables evaluation
! In  l_temp           : for temperature
! In  l_meta           : for metallurgy
! In  varcRefe         : name of reference command variables vector
! In  varcPrev         : command variables at previous step
! In  varcCurr         : command variables at current step
! In  compor           : name of comportment definition (field)
! In  mult_comp        : multi-comportment (DEFI_COMPOR for PMF)
! In  chsith           : commande variable for temperature in XFEM
! In  sigm             : stress
! In  vari             : internal variables
! In  nbFieldInMax     : maximum number of input fields
! In  nbFieldOutMax    : maximum number of output fields
! Out nbFieldIn        : effective number of input fields
! Out nbFieldOut       : effective number of output fields
! In  vect_elem        : name of elementary vectors
! Out lpain            : list of input parameters
! Out lchin            : list of input fields
! Out lpaout           : list of output parameters
! Out lchout           : list of output fields
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: numeHarm = 0
    integer(kind=8) :: iret
    aster_logical :: lXFEM
    character(len=8) :: model
    character(len=24) :: caraElem, materCode
    character(len=19) :: modelLigrel
    character(len=24) :: chgeom, chharm
    character(len=24) :: varcAllRefe, varcAllPrev, varcAllCurr, timeCurr, timePrev
!
! --------------------------------------------------------------------------------------------------
!
    model = modelZ
    caraElem = caraElemZ
    materCode = materCodeZ

! - Initializations
    nbFieldIn = 0
    nbFieldOut = 0
    lpaout = ' '
    lpain = ' '
    lchout = ' '
    lchin = ' '
    call exixfe(model, iret)
    lXFEM = iret .ne. 0
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! - Get fields for external state variables
    call nmvcex('TOUT', varcRefeZ, varcAllRefe)
    if (poum .eq. '-') then
        call nmvcex('TOUT', varcPrevZ, varcAllPrev)
    elseif (poum .eq. '+') then
        call nmvcex('TOUT', varcPrevZ, varcAllPrev)
        call nmvcex('TOUT', varcCurrZ, varcAllCurr)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get fields for time
    if (poum .eq. '-') then
        call nmvcex('INST', varcPrevZ, timePrev)
    elseif (poum .eq. '+') then
        call nmvcex('INST', varcCurrZ, timeCurr)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get geometry field
    call megeom(model, chgeom)

! - Create field for Fourier
    call meharm(model, numeHarm, chharm)

! - Add input fields
    lpain(1) = 'PVARCRR'
    lchin(1) = varcAllRefe(1:19)
    lpain(2) = 'PGEOMER'
    lchin(2) = chgeom(1:19)
    lpain(3) = 'PMATERC'
    lchin(3) = materCode(1:19)
    lpain(4) = 'PCONTMR'
    lchin(4) = sigmz(1:19)
    lpain(5) = 'PVARIPR'
    lchin(5) = variz(1:19)
    lpain(6) = 'PHARMON'
    lchin(6) = chharm(1:19)
    nbFieldIn = 6

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add behaviour field => only for metallurgy (non-linear)
    nbFieldIn = nbFieldIn+1
    lpain(nbFieldIn) = 'PCOMPOR'
    if (l_meta) then
        lchin(nbFieldIn) = comporZ(1:19)
    else
        lchin(nbFieldIn) = multCompZ(1:19)
    end if

! - Computation of elementary vectors - Previous
    if (poum .eq. '-') then
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PINSTR'
        lchin(nbFieldIn) = timePrev(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVARCPR'
        lchin(nbFieldIn) = varcAllPrev(1:19)

    elseif (poum .eq. '+') then
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PINSTR'
        lchin(nbFieldIn) = timeCurr(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVARCPR'
        lchin(nbFieldIn) = varcAllCurr(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVARCMR'
        lchin(nbFieldIn) = varcAllPrev(1:19)

    else
        ASSERT(ASTER_FALSE)
    end if

! - Add XFEM input fields
    if (lXFEM .and. l_temp) then
        call xajcin(model, 'CHAR_MECA_TEMP_R', nbFieldInMax, lchin, lpain, nbFieldIn)
    end if

! - Add HHO fields
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add output fields
    lpaout(1) = 'PVECTUR'
    nbFieldOut = 1

! - Add XFEM output field
    if (lXFEM .and. l_temp) then
        call detrsd('CHAM_ELEM', chsithz)
        call alchml(modelLigrel, 'SIEF_ELGA', 'PCONTRR', 'V', chsithz, iret, ' ')
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PCONTRT'
    end if
!
end subroutine
