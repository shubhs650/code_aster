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
subroutine pascou(materField, materCode, caraElem, sddyna, sddisc)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/megeom.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "jeveux.h"
!
    character(len=24), intent(in) :: materField, materCode, caraElem
    character(len=19), intent(in) :: sddyna, sddisc
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE DYNA_NON_LINE (UTILITAIRE)
!
! EVALUATION DU PAS DE TEMPS DE COURANT POUR LE MODELE
!
! --------------------------------------------------------------------------------------------------
!
! IN  MATE   : CHAMP MATERIAU
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! IN  SDDYNA : SD DEDIEE A LA DYNAMIQUE (CF NDLECT)
! IN  SDDISC : SD DISCRETISATION
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    integer(kind=8), parameter :: numeInstInit = 0
    integer(kind=8) :: ibid, jcesd, jcesl, n1, i
    integer(kind=8) :: nbCell, iCell, iad, nbinst, cellNume
    real(kind=8) :: dtcou, valeur, phi, timeInit
    aster_logical :: booneg, boopos
    character(len=2) :: codret
    character(len=8) :: model, stopCFL, mesh
    character(len=19), parameter :: chams = '&&OP0070.CHAMS'
    character(len=19), parameter :: chvarc = '&&PASCOU.CH_VARC_R'
    character(len=24) :: chgeom, modelLigrel
    real(kind=8), pointer :: ditr(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - INITIALISATIONS
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '
    call getvid(' ', 'MODELE', scal=model, nbret=ibid)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NB_MA_MAILLA', modelLigrel, 'LIGREL', repi=nbCell)

! - CHAMP DES VARIABLES DE COMMANDE
    timeInit = diinst(sddisc, numeInstInit)
    call vrcins(model, materField, caraElem, timeInit, chvarc, codret)

! - Add input fields
    call megeom(model, chgeom)
    lpain(1) = 'PMATERC'
    lchin(1) = materCode(1:19)
    lpain(2) = 'PGEOMER'
    lchin(2) = chgeom(1:19)
    lpain(3) = 'PVARCPR'
    lchin(3) = chvarc(1:19)
    nbFieldIn = 3

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add output field
    lpaout(1) = 'PCOURAN'
    lchout(1) = '&&OP0070.PASCOU'

! - Compute
    call calcul('S', 'PAS_COURANT', modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')

!  PASSAGE D'UN CHAM_ELEM EN UN CHAM_ELEM_S
    call celces(lchout(1), 'V', chams)
    call jeveuo(chams//'.CESD', 'L', jcesd)
    call jeveuo(chams//'.CESL', 'L', jcesl)
    call jeveuo(chams//'.CESV', 'L', vr=cesv)

! A L'ISSUE DE LA BOUCLE :
! BOONEG=TRUE SI L'ON N'A PAS PU CALCULER DTCOU POUR AU MOINS UN ELMNT
! BOOPOS=TRUE SI L'ON A CALCULE DTCOU POUR AU MOINS UN ELEMENT
    dtcou = -1.d0
    booneg = ASTER_FALSE
    boopos = ASTER_FALSE
    cellNume = 1
    do iCell = 1, nbCell
        call cesexi('C', jcesd, jcesl, iCell, 1, 1, 1, iad)
        if (iad .gt. 0) then
            valeur = cesv(iad)
        else if (iad .eq. 0) then
            cycle
        end if
        if (valeur .lt. 0) then
            booneg = ASTER_TRUE
        else
            boopos = ASTER_TRUE
            if (dtcou .gt. 0) then
                if (valeur .le. dtcou) then
                    dtcou = valeur
                    cellNume = iCell
                end if
            else
                dtcou = valeur
            end if
        end if
    end do
!
    call getvtx('SCHEMA_TEMPS', 'STOP_CFL', iocc=1, scal=stopCFL, nbret=n1)
!
! BOOPOS=TRUE SI L'ON A CALCULE DTCOU POUR AU MOINS UN ELEMENT
    if (boopos) then
        if (booneg) then
            call utmess('A', 'DYNAMIQUE_3')
        end if

! ----- Output
        if (ndynlo(sddyna, 'DIFF_CENT')) then
            dtcou = dtcou/(2.d0)
            call utmess('I', 'DYNAMIQUE_5', si=cellNume, sr=dtcou)
        else
            if (ndynlo(sddyna, 'TCHAMWA')) then
                phi = ndynre(sddyna, 'PHI')
                dtcou = dtcou/(phi*2.d0)
                call utmess('I', 'DYNAMIQUE_6', si=cellNume, sr=dtcou)
            else
                call utmess('F', 'DYNAMIQUE_1')
            end if
        end if

!       VERIFICATION DE LA CONFORMITE DE LA LISTE D'INSTANTS
        call utdidt('L', sddisc, 'LIST', 'NBINST', vali_=nbinst)
        call jeveuo(sddisc//'.DITR', 'L', vr=ditr)
        do i = 1, nbinst-1
            if (ditr(i+1)-ditr(i) .gt. dtcou) then
                if (stopCFL(1:3) .eq. 'OUI') then
                    call utmess('F', 'DYNAMIQUE_2')
                else
                    call utmess('A', 'DYNAMIQUE_2')
                end if
            end if
        end do
    else if (stopCFL(1:3) .eq. 'OUI') then
        call utmess('F', 'DYNAMIQUE_4')
    else if (stopCFL(1:3) .eq. 'NON') then
        call utmess('A', 'DYNAMIQUE_4')
    end if
!
    call jedema()
!
end subroutine
