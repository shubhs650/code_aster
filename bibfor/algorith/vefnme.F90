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
subroutine vefnme(optionZ, modelZ, materCode, caraElem, &
                  compor, nh, ligrelInZ, &
                  varcz, sigmz, strxz, dispz, &
                  jvBase, vectElemZ)
!
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/maveElemCreate.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/reajre.h"
#include "asterfort/setStructFields.h"
#include "asterfort/xajcin.h"
!
    character(len=*), intent(in) :: optionZ, modelZ
    character(len=24), intent(in) :: materCode, caraElem
    character(len=19), intent(in) :: compor
    integer(kind=8), intent(in) :: nh
    character(len=*), intent(in) :: ligrelInZ
    character(len=*), intent(in) :: sigmz, varcz, strxz, dispz
    character(len=1), intent(in) :: jvBase
    character(len=*), intent(in) :: vectElemZ
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Option: FORC_NODA
!         FONL_NOEU
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : NOM DU MODELE (NECESSAIRE SI SIGMA EST UNE CARTE)
! IN  SIGMA  : NOM DU CHAM_ELEM (OU DE LA CARTE) DE CONTRAINTES
! IN  CARA   : NOM DU CARA_ELEM
! IN  DEPMOI : NOM DU CHAM_NO DE DEPLACEMENTS PRECEDENTS
! IN  DEPDEL : NOM DU CHAM_NO D'INCREMENT DEPLACEMENTS
! IN  MATCOD : NOM DU MATERIAU CODE
! IN  COMPOR : NOM DE LA CARTE DE COMPORTEMENT
! IN  NH     : NUMERO D'HARMONIQUE DE FOURIER
! IN  PARTPS : INSTANT PRECEDENT ET ACTUEL
! IN  CARCRI : CARTE DES CRITERES ET DE THETA
! IN  CHVARC : NOM DU CHAMP DE VARIABLE DE COMMANDE
! IN  LIGREZ : (SOUS-)LIGREL DE MODELE POUR CALCUL REDUIT
!                  SI ' ', ON PREND LE LIGREL DU MODELE
! OUT VECELZ : VECT_ELEM RESULTAT.
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    character(len=19), parameter :: chharm = '&&VEFNME.NUME_HARM'
    aster_logical :: lXFEM
    character(len=8) :: mesh, newnom, model
    character(len=16) :: option
    character(len=19) :: vectElem, resuElem
    character(len=19) :: ligrelCalc, ligrelIn
    character(len=19) :: chgeom
    integer(kind=8) :: iret, nbFieldIn
    character(len=19) :: sigm, varc, strx, disp
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "
    model = modelZ
    sigm = sigmz
    varc = varcz
    strx = strxz
    disp = dispz
    option = optionZ
    call exixfe(model, iret)
    lXFEM = (iret .eq. 1)

! - Get FED to compute
    ligrelCalc = " "
    ligrelIn = ligrelInZ
    if (ligrelIn .eq. ' ') then
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrelCalc)
    else
        ligrelCalc = ligrelIn
    end if

! - Get mesh
    if (disp .ne. ' ') then
        call dismoi('NOM_MAILLA', disp, 'CHAM_NO', repk=mesh)
    else if (sigm .ne. ' ') then
        call dismoi('NOM_MAILLA', sigm, 'CHAM_ELEM', repk=mesh)
    else
        ASSERT(ASTER_FALSE)
    end if
    chgeom = mesh(1:8)//'.COORDO'

! - Create field for Fourier mode
    call mecact('V', chharm, 'MAILLA', mesh, 'HARMON', &
                ncmp=1, nomcmp='NH', si=nh)

! - Suppress old vectElem result
    vectElem = vectElemZ
    call maveElemCreate(jvBase, vectElem, model)

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = materCode(1:19)
    lpain(3) = 'PCOMPOR'
    lchin(3) = compor
    lpain(4) = 'PSIEFR'
    lchin(4) = sigm
    lpain(5) = 'PDEPLAR'
    lchin(5) = disp
    lpain(6) = 'PHARMON'
    lchin(6) = chharm
    lpain(7) = 'PVARCPR'
    lchin(7) = varc
    lpain(8) = 'PSTRXMR'
    lchin(8) = strx
    nbFieldIn = 8

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add XFEM fields
    if (lXFEM) then
        call xajcin(model, option, nbFieldInMax, lchin, lpain, nbFieldIn)
    end if

! - Add HHO field
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Set output field
    newnom = '.0000000'
    resuElem = vectElem(1:8)//'.0000000'
    call gcnco2(newnom)
    resuElem(10:16) = newnom(2:8)
    call corich('E', resuElem, ichin_=-1)
    lpaout(1) = 'PVECTUR'
    lchout(1) = resuElem

! - Computation
    call calcul('S', option, ligrelCalc, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')

! - Copying output field
    call reajre(vectElem, resuElem, jvBase)

! - Clean
    call detrsd('CHAMP_GD', chharm)
!
    call jedema()
end subroutine
