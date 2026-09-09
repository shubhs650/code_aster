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
subroutine vechmx(model, listLoad, iLoad, nbLoadIndx, listLoadIndxJv, &
                  nbFieldInMax, lpain, lchin, nbFieldIn, vectElem)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/corich.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisllc.h"
#include "asterfort/lisltc.h"
#include "asterfort/lisopt.h"
#include "asterfort/reajre.h"
!
    character(len=8), intent(in) :: model
    character(len=19), intent(in) :: listLoad
    integer(kind=8), intent(in) :: iLoad, nbLoadIndx
    character(len=24), intent(in) :: listLoadIndxJv
    integer(kind=8), intent(in) :: nbFieldInMax
    character(len=8), intent(inout) :: lpain(nbFieldInMax)
    character(len=19), intent(inout) :: lchin(nbFieldInMax)
    integer(kind=8), intent(inout) :: nbFieldIn
    character(len=19), intent(in) :: vectElem
!
! --------------------------------------------------------------------------------------------------
!
! CALCUL DES VECTEURS ELEMENTAIRES DES CHARGEMENTS MECANIQUES
! DE NEUMANN (VOIR DEFINITION DANS LISDEF)
!
! CALCUL EFFECTIF - BOUCLE SUR LES TYPES DE CHARGEMENT
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMO   : NOM DU MODELE
! IN  LISCHA : SD LISTE DES CHARGES
! IN  ICHAR  : INDICE DE LA CHARGE
! IN  NOMLIS : LISTE DES INDEX DES CHARGES
! IN  NBCH   : LONGUEUR DE NOMLIS
! IN  NBIN_MAXI   : NOMBRE MAXI DE CHAMPS D'ENTREE
! IN  LPAIN  : LISTE DES PARAMETRES IN
! IN  LCHIN  : LISTE DES CHAMPS IN
! IN  LASTIN : NOMBRE EFFECTIF DE CHAMPS IN
! OUT VECELE : VECT_ELEM RESULTAT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOut = 1
    character(len=8) :: lpaout(nbFieldOut)
    character(len=19) :: lchout(nbFieldOut)
    integer(kind=8) :: iLoadIndx, iret, loadInx, nbFieldInMod
    character(len=16) :: option
    character(len=8) :: parain, paraou, newnom
    character(len=8) :: loadType
    character(len=19) :: carte, ligrelCalc
    character(len=13) :: loadPreObject
    integer(kind=8), pointer :: listLoadInx(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    newnom = '.0000000'

! - PREFIXE DE L'OBJET DE LA CHARGE
    call lisllc(listLoad, iLoad, loadPreObject)

! - TYPE DE LA CHARGE
    call lisltc(listLoad, iLoad, loadType)

! - Generate name of output field
    call codent(iLoad, 'D0', newnom(2:8))
    lchout(1) = '&&VECHMX.'//newnom(2:8)
    call corich('E', lchout(1), ichin_=iLoad)

! - LISTE DES INDEX DES CHARGES
    call jeveuo(listLoadIndxJv, 'L', vi=listLoadInx)

! - CALCUL
    do iLoadIndx = 1, nbLoadIndx
        loadInx = listLoadInx(iLoadIndx)
        call lisopt(loadPreObject, model, loadType, loadInx, option, &
                    parain, paraou, carte, ligrelCalc)
        call jeexin(carte(1:19)//'.DESC', iret)
        if (iret .ne. 0) then

! --------- Input field
            nbFieldInMod = nbFieldIn+1
            lchin(nbFieldInMod) = carte
            lpain(nbFieldInMod) = parain
            ASSERT(nbFieldInMod .le. nbFieldInMax)

! --------- Output field
            lpaout(1) = paraou
            call calcul('S', option, ligrelCalc, &
                        nbFieldInMod, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        'V', 'OUI')

! --------- RESU_ELEM DANS LE VECT_ELEM
            call exisd('CHAMP_GD', lchout(1), iret)
            ASSERT(iret .gt. 0)
            call reajre(vectElem, lchout(1), 'V')
        end if
    end do
!
    call jedema()
end subroutine
