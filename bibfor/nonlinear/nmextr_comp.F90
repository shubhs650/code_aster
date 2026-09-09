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
subroutine nmextr_comp(field, field_disc, field_type, &
                       meshZ, modelZ, caraElemZ, &
                       ds_material, ds_constitutive, &
                       dispCurrZ, strxCurrZ, varcCurrZ, time, ligrelZ_)
!
    use NonLin_Datastructure_type
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/mecact.h"
#include "asterfort/setStructFields.h"
!
    character(len=19), intent(in) :: field
    character(len=24), intent(in) :: field_type
    character(len=4), intent(in) :: field_disc
    character(len=*), intent(in) :: modelZ, meshZ, caraElemZ
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=*), intent(in) :: dispCurrZ, strxCurrZ, varcCurrZ
    real(kind=8), intent(in) :: time
    character(len=*), optional, intent(in) :: ligrelZ_
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Field extraction datastructure
!
! Compute fields when not a default in nonlinear operator
!
! ONLY EPSI_ELGA !
!
! --------------------------------------------------------------------------------------------------
!
! In  field            : name of field
! In  field_disc       : localization of field (discretization: NOEU or ELGA)
! In  field_type       : type of field (name in results datastructure)
! In  model            : name of model
! In  mesh             : name of mesh
! In  caraElem         : name of datastructure for elementary parameters (CARTE)
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  dispCurr         : current displacements
! In  varcCurr         : command variable for current time
! In  time             : current time
! In  strxCurr         : fibers information for current time
! In  ligrel           : current LIGREL (if not present: on all model)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
!
    character(len=16), parameter :: option = 'EPSI_ELGA'
    character(len=24), parameter :: chtime = '&&NMEXTR_COMP.CHTIME'
    integer(kind=8), parameter :: numeHarm = 0
    character(len=24) :: chgeom, chharm
    character(len=19) :: ligrel, modelLigrel
    integer(kind=8) :: nbFieldIn
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(field_type .eq. 'EPSI_ELGA')
    ASSERT(field_disc .eq. 'ELGA')

! - Initializations
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)
    if (present(ligrelZ_)) then
        ligrel = ligrelZ_
    else
        ligrel = modelLigrel
    end if

! - Create time field
    call mecact('V', chtime, 'MAILLA', meshZ, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=time)

! - Get geometry field
    call megeom(modelZ, chgeom)

! - Create Fourier field
    call meharm(modelZ, numeHarm, chharm)

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PDEPLAR'
    lchin(2) = dispCurrZ(1:19)
    lpain(3) = 'PMATERC'
    lchin(3) = ds_material%mateco(1:19)
    lpain(4) = 'PINSTR'
    lchin(4) = chtime(1:19)
    lpain(5) = 'PVARCPR'
    lchin(5) = varcCurrZ(1:19)
    lpain(6) = 'PVARCRR'
    lchin(6) = ds_material%varc_refe(1:19)
    lpain(7) = 'PCOMPOR'
    lchin(7) = ds_constitutive%compor(1:19)
    lpain(8) = 'PHARMON'
    lchin(8) = chharm(1:19)
    lpain(9) = 'PSTRXMR'
    lchin(9) = strxCurrZ(1:19)
    nbFieldIn = 9

! - Add fields for structural elements
    call setStructFields(caraElemZ, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElemZ)

! - Set output field
    lpaout(1) = 'PDEFOPG'
    lchout(1) = field

! - Computation
    call calcul('S', option, ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')
!
end subroutine
