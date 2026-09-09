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
subroutine vefnme_cplx(optionZ, jvBase, &
                       model, materCode, caraElem, &
                       comporZ, timePrev, timeCurr, nh, ligrelInZ, varcZ, &
                       sigmPrev, sigmZ, strxz, deplz, vectElemZ)
!
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/exixfe.h"
#include "asterfort/gcnco2.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/reajre.h"
#include "asterfort/sepach.h"
#include "asterfort/setStructFields.h"
#include "asterfort/vemare.h"
#include "asterfort/xajcin.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: optionZ
    character(len=1), intent(in) :: jvBase
    character(len=8), intent(in) :: model
    real(kind=8), intent(in) :: timePrev, timeCurr
    character(len=8), intent(in) :: caraElem
    character(len=24), intent(in) :: materCode
    character(len=*), intent(in) :: ligrelInZ
    integer(kind=8), intent(in) :: nh
    character(len=*), intent(in) :: comporZ, sigmZ, sigmPrev, varcZ, strxz, deplz
    character(len=*), intent(inout) :: vectElemZ(2)
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
    character(len=19), parameter :: chharm = '&&VEFNME.NUME_HARM'
    character(len=19), parameter :: chtimePrev = '&&VEFNME.CH_INSTAM'
    character(len=19), parameter :: chtimeCurr = '&&VEFNME.CH_INSTAP'
    integer(kind=8), parameter :: nbFieldOut = 1
    integer(kind=8), parameter :: nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchoutReal(nbFieldOut), lchoutCplx(2)
    character(len=19) :: lchin(nbFieldInMax), lchinReal(nbFieldInMax), lchinImag(nbFieldInMax)
    character(len=19) :: chdecr(nbFieldInMax), chdeci(nbFieldInMax)
!
    integer(kind=8) :: fieldIsCplx(nbFieldInMax)
    aster_logical :: lXFEM
    character(len=8) :: mesh, newnom, physQuanName
    character(len=19) :: ligrelCalc, ligrelIn
    character(len=19) :: chgeom, vectElemR, vectElemI, compor, fieldIn
    character(len=16) :: option
    integer(kind=8) :: iret, iexi, iFieldIn, nbFieldIn
    character(len=19) :: sigm, varc, strx, depl
    character(len=19) :: fieldInReal, fieldInCplx
    aster_logical :: lcmplx, lsspt
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lpain = " "
    lpaout = " "
    lchin = " "
    lchinReal = " "
    lchinImag = " "
    lchoutReal = " "
    lchoutCplx = " "
    sigm = sigmZ
    varc = varcZ
    strx = strxz
    depl = deplz
    compor = comporZ
    newnom = '.0000000'
    option = optionZ
    if (optionZ .ne. 'FONL_NOEU') then
        option = 'FORC_NODA'
    end if
    call exixfe(model, iret)
    lXFEM = (iret .eq. 1)

! - Get mesh
    if (depl .ne. ' ') then
        call dismoi('NOM_MAILLA', depl, 'CHAM_NO', repk=mesh)
    else if (sigm .ne. ' ') then
        call dismoi('NOM_MAILLA', sigm, 'CHAM_ELEM', repk=mesh)
    else
        ASSERT(ASTER_FALSE)
    end if
    chgeom = mesh(1:8)//'.COORDO'

! - Select vectElem
    vectElemR = vectElemZ(1)
    if (vectElemR .eq. ' ') then
        vectElemR = '&&VEFNME'
    end if
    vectElemI = vectElemZ(2)
    if (vectElemI .eq. ' ') then
        vectElemI = '&&VEFNMI'
    end if

! - Get FED to compute
    ligrelCalc = " "
    ligrelIn = ligrelInZ
    if (ligrelIn .eq. ' ') then
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrelCalc)
    else
        ligrelCalc = ligrelIn
    end if
!
! - <CARTE> for structural elements
!
    !   call mecara(caraElem, chcara)
!
! - Create field for Fourier mode
    call mecact('V', chharm, 'MAILLA', mesh, 'HARMON', &
                ncmp=1, nomcmp='NH', si=nh)

! - Create fields for time
    call mecact('V', chtimePrev, 'MAILLA', mesh, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=timePrev)
    call mecact('V', chtimeCurr, 'MAILLA', mesh, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=timeCurr)

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = materCode(1:19)
    lpain(3) = 'PCOMPOR'
    lchin(3) = compor
    lpain(4) = 'PSIEFR'
    lchin(4) = sigmZ
    lpain(5) = 'PDEPLAR'
    lchin(5) = depl
    lpain(6) = 'PCONTGM'
    lchin(6) = sigmPrev
    lpain(7) = 'PHARMON'
    lchin(7) = chharm
    lpain(8) = 'PINSTMR'
    lchin(8) = chtimePrev
    lpain(9) = 'PINSTPR'
    lchin(9) = chtimeCurr
    lpain(10) = 'PVARCPR'
    lchin(10) = varc
    lpain(11) = 'PSTRXMR'
    lchin(11) = strx
    nbFieldIn = 11

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

! - Detect complex case
    lcmplx = ASTER_FALSE
    chdeci = " "
    chdecr = " "
    do iFieldIn = 1, nbFieldIn
        fieldIsCplx(iFieldIn) = 0
        fieldIn = lchin(iFieldIn)
        if (fieldIn .eq. ' ') cycle
        ASSERT(fieldIn .ne. ' ')
        call exisd('CHAMP', fieldIn, iexi)
        if (iexi .eq. 0) cycle
        call dismoi('NOM_GD', fieldIn, 'CHAMP', repk=physQuanName)
        if (physQuanName(5:6) .eq. '_C') then
            lcmplx = ASTER_TRUE
            fieldIsCplx(iFieldIn) = 1
            fieldInReal = '&&VEFNME.CHXX.R'
            fieldInCplx = '&&VEFNME.CHXX.I'
            call codent(iFieldIn, 'D0', fieldInReal(12:13))
            call codent(iFieldIn, 'D0', fieldInCplx(12:13))
            call sepach(caraElem, fieldIn, 'V', fieldInReal, fieldInCplx)
            chdecr(iFieldIn) = fieldInReal
            chdeci(iFieldIn) = fieldInCplx
        end if
    end do

! - Set output fields
    lpaout(1) = 'PVECTUR'
    if (lcmplx) then
        call gcnco2(newnom)
        lchoutCplx(1) = vectElemR(1:8)//newnom
        call corich('E', lchoutCplx(1), ichin_=-1)
        call gcnco2(newnom)
        lchoutCplx(2) = vectElemI(1:8)//newnom
        call corich('E', lchoutCplx(2), ichin_=-1)
    else
        call gcnco2(newnom)
        lchoutReal(1) = vectElemR(1:8)//newnom
        call corich('E', lchoutReal(1), ichin_=-1)
    end if
    call exisd('CHAM_ELEM_S', lchoutReal(1), iexi)
    lsspt = (iexi .ne. 0)

! - Suppress old vectElem result
    call detrsd('VECT_ELEM', vectElemR)
    if (lcmplx) then
        call detrsd('VECT_ELEM', vectElemI)
    end if
    call vemare(jvBase, vectElemR, model)
    if (lcmplx) then
        do iFieldIn = 1, nbFieldIn
            if (fieldIsCplx(iFieldIn) .eq. 0) then
                lchinReal(iFieldIn) = lchin(iFieldIn)
                lchinImag(iFieldIn) = lchin(iFieldIn)
            else
                lchinReal(iFieldIn) = chdecr(iFieldIn)
                lchinImag(iFieldIn) = chdeci(iFieldIn)
            end if
        end do
        call vemare(jvBase, vectElemI, model)
    end if

! - APPEL A CALCUL
    if (lcmplx) then
        if (lsspt) then
            call copisd('CHAM_ELEM_S', 'V', lchoutReal(1), lchoutCplx(1))
        end if
        call calcul('S', option, ligrelCalc, &
                    nbFieldIn, lchinReal, lpain, &
                    nbFieldOut, lchoutCplx(1), lpaout, &
                    jvBase, 'OUI')
        call reajre(vectElemR, lchoutCplx(1), jvBase)
        vectElemZ(1) = vectElemR//'.RELR'
        if (lsspt) then
            call copisd('CHAM_ELEM_S', 'V', lchoutReal(1), lchoutCplx(2))
        end if
        call calcul('S', option, ligrelCalc, &
                    nbFieldIn, lchinImag, lpain, &
                    nbFieldOut, lchoutCplx(2), lpaout, &
                    jvBase, 'OUI')
        call reajre(vectElemI, lchoutCplx(2), jvBase)
        vectElemZ(2) = vectElemI//'.RELR'
    else
        call calcul('S', option, ligrelCalc, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchoutReal, lpaout, &
                    jvBase, 'OUI')
        call reajre(vectElemR, lchoutReal(1), jvBase)
        vectElemZ(1) = vectElemR//'.RELR'
    end if

! - Cleaning
    do iFieldIn = 1, nbFieldIn
        if (fieldIsCplx(iFieldIn) .ne. 0) then
            call detrsd('CHAMP', chdecr(iFieldIn))
            call detrsd('CHAMP', chdeci(iFieldIn))
        end if
    end do

    call detrsd('CHAMP_GD', chharm)
    call jedema()
end subroutine
