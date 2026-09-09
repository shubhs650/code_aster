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

! --------------------------------------------------------------------
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your optionZ) any later version.
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
subroutine memam2(optionZ, &
                  modelZ, materFieldZ, materCodeZ, caraElemZ, &
                  compor, time, chacceZ, &
                  vectElemZ, jvBaseZ, ligrelZ)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: optionZ, modelZ, materFieldZ, materCodeZ, caraElemZ
    character(len=24), intent(in) :: compor
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: chacceZ, vectElemZ, jvBaseZ, ligrelZ
!
! --------------------------------------------------------------------------------------------------
!
!     CALCULE LES VECTEURS ELEMENTAIRES ( MASSE_MECA * CHACCE )
!
! --------------------------------------------------------------------------------------------------
!
! IN  : OPTION : OPTION DE CALCUL
! IN  : MODELE : NOM DU MODELE (OBLIGATOIRE)
! IN  : MATE   : CARTE DE MATERIAUX
! IN  : CARA   : CHAMP DE CARAC_ELEM
! IN  : TIME   : INSTANT DE CALCUL
! IN  : CHACCE : CHAMP D'ACCELERATION
! IN  : VECEL  : NOM DU VECT_ELEM RESULTAT
! IN  : BASEZ  : NOM DE LA BASE
! IN  : LIGREZ  : (SOUS-)LIGREL DE MODELE POUR CALCUL REDUIT
!                  SI ' ', ON PREND LE LIGREL DU MODELE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    character(len=19), parameter :: chvarc = '&&MEMAM2.VARC'
    integer(kind=8), parameter :: numeHarm = 0
    character(len=1) :: jvBase
    character(len=2) :: codret
    character(len=8) :: newnom
    character(len=24) :: ligrel, chgeom, chharm, vectElem, resuElem
    integer(kind=8) :: icode, iret
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    newnom = '.0000000'
    vectElem = vectElemZ
    jvBase = jvBaseZ
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '
    if (modelZ(1:1) .eq. ' ') then
        call utmess('F', 'CALCULEL2_82')
    end if
    ligrel = ligrelZ
    if (ligrel .eq. ' ') then
        call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=ligrel)
    end if

! - Preparation of input fields
    call mecham('MASS_MECA', modelZ, numeHarm, &
                chgeom, chharm, icode)

! - Get external state variable
    call vrcins(modelZ, materFieldZ, ' ', time, chvarc, codret)

! - Prepare MATR_ELEM
    call memare(jvBase, vectElemZ, modelZ, optionZ, ASTER_TRUE)
!
    call jeexin(vectElem(1:19)//'.RELR', iret)
    if (iret .gt. 0) call jedetr(vectElem(1:19)//'.RELR')
    if (icode .eq. 1) then
        goto 10
    end if

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = materCodeZ
    lpain(3) = 'PVARCPR'
    lchin(3) = chvarc
    lpain(4) = 'PACCELR'
    lchin(4) = chacceZ
    lpain(5) = 'PCOMPOR'
    lchin(5) = compor(1:19)
    nbFieldIn = 5

! - Add fields for structural elements
    call setStructFields(caraElemZ, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElemZ)

! - Add output field
    lpaout(1) = 'PVECTUR'
    resuElem = '&&MEMAM2.???????'
    call gcnco2(newnom)
    resuElem(10:16) = newnom(2:8)
    lchout(1) = resuElem(1:19)

    call corich('E', resuElem, ichin_=-1)
    call calcul('S', optionZ, ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')
!
    call reajre(vectElem, lchout(1), jvBase)
!
10  continue
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
