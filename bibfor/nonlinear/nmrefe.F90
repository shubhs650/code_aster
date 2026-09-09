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
subroutine nmrefe(model, compor, materCode, caraElem, nume_dof, &
                  ds_conv, hatValinc, hatVeelem, hatVeasse)
!
    use NonLin_Datastructure_type
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assmiv.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/nmchex.h"
#include "asterfort/reajre.h"
#include "asterfort/vemare.h"
#include "asterfort/setStructFields.h"
#include "asterfort/xajcin.h"
!
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: compor
    character(len=24), intent(in) :: materCode
    character(len=24), intent(in) :: caraElem
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_Conv), intent(in) :: ds_conv
    character(len=19), intent(in) :: hatValinc(*)
    character(len=19), intent(in) :: hatVeelem(*)
    character(len=19), intent(in) :: hatVeasse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Computation
!
! Compute reference vector for RESI_REFE_RELA
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  compor           : name of comportment definition (field)
! In  materCode        : name of coded material
! In  caraElem         : name of elementary characteristics (field)
! In  nume_dof         : name of numbering (NUME_DDL)
! In  ds_conv          : datastructure for convergence management
! In  hatValinc        : hat variable for algorithm fields
! In  hatVeelem        : hat variable for elementary vectors
! In  hatVeasse        : hat variable for vectors
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    character(len=16), parameter :: option = 'REFE_FORC_NODA'
    character(len=19) :: vectElem, vectAsse, dispPrev
    character(len=19) :: modelLigrel
    character(len=19), parameter :: resuElem = '&&NMREFE.VEREFE'
    character(len=24) :: chgeom
    integer(kind=8) :: nbFieldIn, ier
    aster_logical :: lXFEM
!
! --------------------------------------------------------------------------------------------------
!

! - Initializations
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call exixfe(model, ier)
    lXFEM = ier .ne. 0

! - Get names of fields
    call nmchex(hatValinc, 'VALINC', 'DEPMOI', dispPrev)
    call nmchex(hatVeelem, 'VEELEM', 'CNREFE', vectElem)
    call nmchex(hatVeasse, 'VEASSE', 'CNREFE', vectAsse)

! - Get geometry field
    call megeom(model, chgeom)

! - Preparation of VECT_ELEM
    call detrsd('VECT_ELEM', vectElem)
    call vemare('V', vectElem, model)

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PRESICMP'
    lchin(2) = ds_conv%cresicmp
    lpain(3) = 'PRESIREF'
    lchin(3) = ds_conv%cresiref
    lpain(4) = 'PCOMPOR'
    lchin(4) = compor(1:19)
    lpain(5) = 'PMATERC'
    lchin(5) = materCode(1:19)
    lpain(6) = 'PDEPLMR'
    lchin(6) = dispPrev
    nbFieldIn = 6

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

! - Add HHO field
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Set output field
    lpaout(1) = 'PVECTUR'
    lchout(1) = resuElem

! - Computation
    call calcul('S', option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')

! - Copying output field
    call reajre(vectElem, lchout(1), 'V')

! - Assembly
    call assmiv('V', vectAsse, 1, vectElem, [1.d0], &
                nume_dof, 1)
!
end subroutine
