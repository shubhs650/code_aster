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
! ==================================================================================================
!
! Module to compute External State Variables
!
! ==================================================================================================
!
module ExternalStateVariableComp_module
! ==================================================================================================
! ==================================================================================================
    use coorSyst_module, only: setOrieFields
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: varcCompElem
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/reajre.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
#include "asterfort/varcDetect.h"
#include "asterfort/vemare.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! varcCompElem
!
! Compute elementary vector for external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine varcCompElem(line, &
                            numeHarm, modelZ, caraElemZ, materFieldZ, materCodeZ, &
                            chtimeZ, varcRefeZ, varcZ, &
                            jvBase, vectElemZ, &
                            lCumul_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: line
        integer(kind=8), intent(in) :: numeHarm
        character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, materCodeZ
        character(len=*), intent(in) :: chtimeZ, varcRefeZ, varcZ
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(in) :: vectElemZ
        aster_logical, optional, intent(in) :: lCumul_
! ----- Locals
        integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
        character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
        character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOut)
        character(len=8) :: newnom
        character(len=16) :: option
        character(len=24) :: modelLigrel, resuElem
        character(len=24) :: chgeom, chharm
        aster_logical :: lTemp, lHydr, lPtot, lSech, lEpsa, lMeta, lCumul
        integer(kind=8) :: nbFieldIn
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(line)

! ----- Initializations
        call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)
        lpain = " "
        lpaout = " "
        lchin = " "
        lchout = " "
        lCumul = ASTER_FALSE
        if (present(lCumul_)) then
            lCumul = lCumul_
        end if

! ----- Linear: only for temperature for the moment
        call varcDetect(materFieldZ, lTemp, lHydr, lPtot, lSech, lEpsa, lMeta)
        if (line) then
            if (lHydr .or. lPtot .or. lSech .or. lEpsa .or. lMeta) then
                call utmess('F', 'SOUSTRUC_18')
            end if
        end if

! ----- Get geometry field
        call megeom(modelZ, chgeom)

! ----- Create field for Fourier
        call meharm(modelZ, numeHarm, chharm)

! ----- Add input fields
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PMATERC'
        lchin(2) = materCodeZ
        lpain(3) = 'PHARMON'
        lchin(3) = chharm
        lpain(4) = 'PINSTR'
        lchin(4) = chtimeZ
        lpain(5) = 'PVARCRR'
        lchin(5) = varcRefeZ
        lpain(6) = 'PVARCPR'
        lchin(6) = varcZ
        nbFieldIn = 6

! ----- Add fields for structural elements
        call setStructFields(caraElemZ, nbFieldInMax, lchin, lpain, nbFieldIn)

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElemZ)

! ----- Set output field
        lpaout(1) = 'PVECTUR'

! ----- Allocate result
        if (.not. lCumul) then
            call detrsd('VECT_ELEM', vectElemZ)
            call vemare(jvBase, vectElemZ, modelZ)
            call reajre(vectElemZ, ' ', jvBase)
        end if
        newnom = '.0000000'
        resuElem = vectElemZ(1:8)//'.0000000'

! ----- Compute
        if (lTemp) then
            call gcnco2(newnom)
            resuElem(10:16) = newnom(2:8)
            call corich('E', resuElem, ichin_=-1)
            lchout(1) = resuElem
            option = 'CHAR_MECA_TEMP_R'
            call calcul('C', option, modelLigrel, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        jvBase, &
                        'OUI')
            call reajre(vectElemZ, resuElem, jvBase)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module ExternalStateVariableComp_module
