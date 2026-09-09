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
subroutine ccchno(option, numeStore, resultIn, resultOut, fieldNameOut, &
                  lRestCell, nbRestCell, restCellJv, &
                  mesh, model, caraElem, jvBase, &
                  ligrel, ligrelHasBeenChanged, codret, &
                  fieldNameIn, ideb_, ifin_, vcham_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ccvrrl.h"
#include "asterfort/celces.h"
#include "asterfort/cescns.h"
#include "asterfort/cesred.h"
#include "asterfort/cnscno.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/inigrl.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    integer(kind=8), intent(in) :: numeStore
    character(len=8), intent(in) :: resultIn, resultOut
    character(len=24), intent(out) :: fieldNameOut
    aster_logical, intent(in) :: lRestCell
    integer(kind=8), intent(in) :: nbRestCell
    character(len=24), intent(in) :: restCellJv
    character(len=8), intent(in) :: mesh, model, caraElem
    character(len=1), intent(in) :: jvBase
    character(len=24), intent(in) :: ligrel
    integer(kind=8), intent(out) :: codret
    aster_logical, intent(in) :: ligrelHasBeenChanged
    character(len=19), intent(in) :: fieldNameIn
    integer(kind=8), optional, intent(in) :: ideb_, ifin_
    character(len=24), optional, intent(in) :: vcham_
!
! --------------------------------------------------------------------------------------------------
!
!  CALC_CHAMP
!
!  Compute option for nodal field
!
! --------------------------------------------------------------------------------------------------
!
! IN  :
!   OPTION  K16  NOM DE L'OPTION A CALCULER
!   NUMORD  I    NUMERO D'ORDRE COURANT
!   RESUIN  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT IN
!   RESUOU  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT OUT
!   MESMAI  K24  NOM DU VECTEUR CONTENANT LES MAILLES SUR LESQUELLES
!                LE CALCUL EST DEMANDE
!   NOMAIL  K8   NOM DU MAILLAGE SUR LEQUEL LE CALCUL EST REALISE
!   MODELE  K8   NOM DU MODELE
!   CARAEL  K8   NOM DU CARAEL
!   BASOPT  K1   BASE SUR LAQUELLE DOIT ETRE CREE LE CHAMP DE SORTIE
!   LIGREL  K24  NOM DU LIGREL
!   OPTIONNEL: NOCHOU K19 NOM DU CHAMP INITIAL
!   OPTIONNEL: IDEB/IFIN   I   INDICES DEBUT/FIN CHAM_NOS SIMULTANES
!   OPTIONNEL: VCHAM  K24 VECTEUR DES NOMS DES CHAM_NOS SIMULTANES
!
! IN/OUT :
!   LICHOU  K8*  LISTE DES CHAMPS OUT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nmaxob = 30
    integer(kind=8) :: adobj(nmaxob)
    character(len=24) :: noobj(nmaxob)
    integer(kind=8) :: ier, ngr, igr, iret
    integer(kind=8) :: nbobj, nbsp, nb, jvRestCell
    aster_logical :: ldist, lRDM
    character(len=8) :: k8b, answer
    character(len=16) :: optionElno
    character(len=19) :: valk(4)
    character(len=19), parameter :: fieldElnoS = '&&CALCOP.CHAMS0'
    character(len=19) :: fieldElno
    character(len=19), parameter :: fieldNodeS = '&&CALCOP.CHAMS1'
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    optionElno = option(1:5)//'ELNO'
    fieldNameOut = fieldNameIn

! - For time "distribution)
    if (present(ideb_) .and. present(ifin_) .and. present(vcham_)) then
        ldist = .True.
        nb = ifin_-ideb_+1
    else
        ldist = .False.
        nb = 0
    end if

! - Try to get ELNO field
    call rsexch(' ', resultIn, optionElno, numeStore, fieldElno, ier)
    if (ier .ne. 0) then
        call rsexch(' ', resultOut, optionElno, numeStore, fieldElno, ier)
    end if

! - Does ELNO field exist ?
    call exisd('CHAMP_GD', fieldElno, ier)
    if (ier .eq. 0) then
        fieldNameOut = ' '
        valk(1) = optionElno
        valk(2) = option
        call utmess('A', 'CALCCHAMP_2', nk=4, valk=valk, si=numeStore)
        if (ldist) call utmess('F', 'PREPOST_18')
        goto 999
    end if

! - Convert ELNO field and reduce it on lsit of cells
    call celces(fieldElno, 'V', fieldElnoS)
    if (lRestCell) then
        call jeveuo(restCellJv, 'L', jvRestCell)
        call cesred(fieldElnoS, nbRestCell, zi(jvRestCell), 0, [k8b], 'V', fieldElnoS)
    end if

! - Check continuity of local coordinate system
    answer = 'NON'
    call dismoi('EXI_RDM', ligrel, 'LIGREL', repk=answer)
    lRDM = answer(1:3) .eq. 'OUI'

!   cette vérification ne doit être faite que dans le cas ou le modèle contient des éléments
!   de structure et que pour certains champs qui sont en repère local
    if (lRDM .and. &
        ((option(1:4) .eq. 'EPSI') .or. (option(1:4) .eq. 'SIGM') .or. &
         (option(1:4) .eq. 'SIEF') .or. &
         (option(1:4) .eq. 'DEGE') .or. (option(1:4) .eq. 'EFGE'))) then
        if (ligrelHasBeenChanged) then
!           pour les coques 3d certaines initialisations sont nécessaires pour pouvoir utiliser
!           les routines de changement de repère propres aux coques 3d
            call dismoi('EXI_COQ3D', ligrel, 'LIGREL', repk=answer)
            if (answer .eq. 'OUI' .and. ligrelHasBeenChanged) then
                call jelira(ligrel(1:19)//'.LIEL', 'NUTIOC', ngr)
                do igr = 1, ngr
                    call inigrl(ligrel, igr, nmaxob, adobj, noobj, nbobj)
                end do
            end if
        end if
        if (caraElem .ne. ' ') then
            call ccvrrl(mesh, model, caraElem, &
                        lRestCell, nbRestCell, restCellJv, &
                        fieldElnoS, codret)
        else
            call utmess('A', 'CALCULEL4_2', sk=option)
        end if
    end if

! - Convert ELNO field
    call cescns(fieldElnoS, ' ', 'V', fieldNodeS, 'A', codret)
    if (ldist .and. (nb .ge. 2)) then
        call cnscno(fieldNodeS, ' ', 'NON', jvBase, fieldNameIn, ' ', iret, &
                    nbz=nb, vchamz=vcham_)
    else
        call cnscno(fieldNodeS, ' ', 'NON', jvBase, fieldNameIn, ' ', iret)
    end if

!   VERIFICATION POUR LES CHAMPS A SOUS-POINT
    call dismoi('MXNBSP', fieldElno, 'CHAM_ELEM', repi=nbsp)
    if ((nbsp .gt. 1) .and. (iret .eq. 1)) then
        valk(1) = optionElno
        valk(2) = option
        call utmess('F', 'CALCULEL4_16', nk=2, valk=valk)
    end if
!
    call detrsd('CHAM_ELEM_S', fieldElnoS)
    call detrsd('CHAM_NO_S', fieldNodeS)
!
999 continue
    call jedema()
!
end subroutine
