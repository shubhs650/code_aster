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
subroutine modirepcham(fieldOut, fieldIn)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/calcul.h"
#include "asterfort/cesvar.h"
#include "asterfort/checkConsistencyLigrel.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecoor.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: fieldOut, fieldIn
!
! --------------------------------------------------------------------------------------------------
!
!     COMMANDE : MODI_REPERE / CHAM_GD
!
!   in
!       fieldIn  : Nom du champ en entrée
!   out
!       fieldOut  : Nom du champ en sortie
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    integer(kind=8) :: ifm, niv, nret, iret
    character(len=8) :: mesh, caraElem, caraElemMesh, caraElemModel
    character(len=16) :: repere
    character(len=19), parameter :: chpass = '&&REPCHA.MATPASS'
    character(len=24) :: ligrel, option
    character(len=24) :: chgeom
    aster_logical :: lreuse, lret
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "
    lreuse = (fieldIn .eq. fieldOut)

!   Définition du repère utilisé
    call getvtx(' ', 'REPERE', scal=repere, nbret=nret)
    if (nret .eq. 0 .and. .not. lreuse) then
        call utmess('F', 'MODELISA3_2')
    else if (repere .ne. 'GLOBAL_UTIL' .and. .not. lreuse) then
        call utmess('F', 'MODELISA3_3')
    end if

!   Lecture du concept CARA_ELEM
    call getvid(' ', 'CARA_ELEM', scal=caraElem, nbret=nret)
    if (nret .eq. 0 .and. .not. lreuse) then
        call utmess('F', 'MODELISA3_7')
    end if

!   Informations sur le champ en entrée.
    call dismoi('NOM_OPTION', fieldIn, 'CHAMP', repk=option)
    call dismoi('NOM_LIGREL', fieldIn, 'CHAMP', repk=ligrel)
    call dismoi('NOM_MAILLA', fieldIn, 'CHAMP', repk=mesh)
    if (option .ne. 'INI_SP_RIGI') then
        call utmess('F', 'MODELISA3_1')
    end if
!
! --------------------------------------------------------------------------------------------------
!   Vérification que CARCOQUE existe
    call exisd('CARTE', caraElem//'.CARCOQUE', iret)
    if (iret .eq. 0) then
        call utmess('F', 'MODELISA3_4')
    end if
!   Vérification que CANBSP existe
    call exisd('CHAM_ELEM', caraElem//'.CANBSP', iret)
    if (iret .eq. 0) then
        call utmess('F', 'MODELISA3_4')
    end if

!   Nom du mesh sous-jacent à la carte. Le même que celui du champ.
    call dismoi('NOM_MAILLA', caraElem, 'CARA_ELEM', repk=caraElemMesh)
    if (mesh .ne. caraElemMesh) then
        call utmess('F', 'MODELISA3_5')
    end if

!   Nom du modèle sous-jacent à CARA_ELEM. Le même que celui du champ.
    call dismoi('NOM_MODELE', caraElem, 'CARA_ELEM', repk=caraElemModel)
    call checkConsistencyLigrel(caraElemModel, ligrel, lret)
    if (.not. lret) then
        call utmess('F', 'MODELISA3_6')
    end if
!
! --------------------------------------------------------------------------------------------------
!   Matrice de passage du repère global vers le repère utilisateur

! - Set input fields
    call mecoor(ligrel, chgeom)
    nbFieldIn = 1
    lpain(nbFieldIn) = 'PGEOMER'
    lchin(nbFieldIn) = chgeom(1:19)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)
! - Set output field
    lchout(1) = chpass
    lpaout(1) = 'PMATPASS'

! - Compute
    call calcul('C', 'REPERE_LOCAL', ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'NON')

! --------------------------------------------------------------------------------------------------
!   Changement de repère

    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Set input fields
    nbFieldIn = 1
    lchin(nbFieldIn) = chpass
    lpain(nbFieldIn) = 'PMATPASS'
    nbFieldIn = nbFieldIn+1
    lchin(nbFieldIn) = fieldIn
    lpain(nbFieldIn) = 'PSIEFR'

! - Set output field
    lchout(1) = fieldOut
    lpaout(1) = 'PCONTPR'
    call cesvar(caraElem, ' ', ligrel, lchout(1))

! - Compute
    call calcul('C', 'MODI_REPERE', ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'G', 'NON')
!
! --------------------------------------------------------------------------------------------------
!   Ménage
    call jedetr(chpass)
!
    call jedema()
end subroutine
