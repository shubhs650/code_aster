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
subroutine carelo(model, caraElem, jvBase, &
                  chrel1, chrel2, chrel3)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/megeom.h"
!
    character(len=8), intent(in) :: caraElem, model
    character(len=1), intent(in) :: jvBase
    character(len=19), intent(in) :: chrel1, chrel2, chrel3
!
! --------------------------------------------------------------------------------------------------
!
!       CALCULER LES REPERES LOCAUX DES ELEMENTS
!
! --------------------------------------------------------------------------------------------------
!
!     IN MODELE  : MODELE
!     IN CARELE  : CARA_ELEM
!     IN BASE    : G/V
!     OUT CHREL1 : 1ER  VECTEUR DU REPERE LOCAL
!     OUT CHREL2 : 2EME VECTEUR DU REPERE LOCAL
!     OUT CHREL3 : 3EME VECTEUR DU REPERE LOCAL
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 3
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    character(len=19) :: modelLigrel, chgeom
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! - Get geometry field
    call megeom(model, chgeom)

! - Set input fields
    lchin(1) = chgeom
    lpain(1) = 'PGEOMER'
    nbFieldIn = 1

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Set output fields
    lchout(1) = chrel1
    lpaout(1) = 'PREPLO1'
    lchout(2) = chrel2
    lpaout(2) = 'PREPLO2'
    lchout(3) = chrel3
    lpaout(3) = 'PREPLO3'

! - Compute
    call calcul('C', 'REPERE_LOCAL', modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'NON')
!
    call jedema()
end subroutine
