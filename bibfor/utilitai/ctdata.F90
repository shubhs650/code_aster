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

subroutine ctdata(mesnoe, mesmai, nkcha, fieldDisc, toucmp, &
                  nkcmp, nkvari, nbcmp, chpgs, chpsu, mesh, &
                  nbNode, nbCell, nbField, physQuanScal)
!
! --------------------------------------------------------------------------------------------------
!
!                 OPERATEUR CREA_TABLE , MOT-CLE FACTEUR RESU
!
!          RECUPERER LES DONNEES UTILES POUR CONSTRUIRE LA TABLE
!              (COMPOSANTES,NOEUDS,MAILLES,...)
!
! --------------------------------------------------------------------------------------------------
!
!        IN     : NKCHA  (K24) : OBJET DES NOMS DE CHAMP
!                 NBVAL (I)    : NOMBRE DE VALEURS D'ACCES
!        IN/OUT : MESNOE (K24) : OBJET DES NOMS DE NOEUD
!                 MESMAI (K24) : OBJET DES NOMS DE MAILLE
!                 NKCMP  (K24) : OBJET DES NOMS DE COMPOSANTES  (NOM_CMP)
!                 NKVARI (K24) : OBJET DES NOMS DE VAR. INTERNES (NOM_VARI)
!                 NCHSPG (K24) : NOM DU CHAM_ELEM_S DES COORDONNES DES
!                                POINTS DE GAUSS (REMPLI SI TYCH='ELGA')
!        OUT    : TYCH   (K4)  : TYPE DE CHAMP (=NOEU,ELXX,CART)
!                 TOUCMP (L)   : INDIQUE SI TOUT_CMP EST RENSEIGNE
!                 NBCMP  (I)   : NOMBRE DE COMPOSANTES LORSQUE
!                                NOM_CMP EST RENSEIGNE, 0 SINON
!                 NDIM   (I)   : DIMENSION GEOMETRIQUE (=2 OU 3)
!                 NOMA   (K8)  : NOM DU MAILLAGE
!                 NBNO   (I)   : NOMBRE DE NOEUDS UTILISATEUR
!                 NBMA   (I)   : NOMBRE DE MAILLES UTILISATEUR
!                 TSCA  (K1)  : TYPE DE LA GRANDEUR (REEL)
!
! --------------------------------------------------------------------------------------------------
!
    use coorSyst_module, only: setOrieFields
    use MGIS_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cesvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/exlim2.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecoor.h"
#include "asterfort/reliem.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/rsGetOneBehaviourFromResult.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
#include "asterfort/varinonu.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8) :: nbcmp, nbNode, nbCell, nbField
    character(len=1) :: physQuanScal
    character(len=4) :: fieldDisc
    character(len=8) :: mesh
    character(len=24) :: mesnoe, mesmai, nkcha, nkvari, nkcmp
    character(len=19) :: chpgs, chpsu
    aster_logical :: toucmp
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    integer(kind=8) :: nbFieldIn, nbFieldOut
    integer(kind=8) :: iField, iret, jlno, n1, jlma, n2, n3, ierr, iNode, iCell
    integer(kind=8) :: igrel, nbVari, nbRet
    integer(kind=8), pointer :: repe(:) => null()
    character(len=8) :: physQuanName, caraElem
    character(len=8) :: typmcl(4), result
    character(len=16) :: motcle(4), fieldName
    character(len=19) :: ligrel, ligrelField, cel19, compor
    character(len=24) :: chgeom, field
    aster_logical :: hasCaraElem
    integer(kind=8), pointer :: listStore(:) => null()
    integer(kind=8) :: nbStore
    character(len=16), pointer :: variName(:) => null()
    character(len=8), pointer :: cmpName(:) => null()
    character(len=24), pointer :: listField(:) => null()
    aster_logical :: lFieldUser, lResultUser
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    fieldDisc = ' '
    ligrel = '&&CTDATA.LIGREL'
    physQuanScal = ' '
    result = ' '
    hasCaraElem = ASTER_FALSE
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '

! - Get result or field from user ?
    call getvid('RESU', 'RESULTAT', iocc=1, scal=result, nbret=nbRet)
    lResultUser = nbRet .ne. 0
    call getvid('RESU', 'CHAM_GD', iocc=1, nbval=0, nbret=nbRet)
    lFieldUser = nbRet .ne. 0

!   DETERMINATION DU TYPE DE CHAMP
    call jeveuo(nkcha, 'L', vk24=listField)

    do iField = 1, nbField
        field = listField(iField)
        if (field .ne. '&&CHAMP_INEXISTANT') then
! --------- Parameters of field
            call dismoi('TYPE_CHAMP', field, 'CHAMP', repk=fieldDisc)
            call dismoi('NOM_GD', field, 'CHAMP', repk=physQuanName)
            call dismoi('NOM_MAILLA', field, 'CHAMP', repk=mesh)
            call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nbNode)
            call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCell)
            call dismoi('TYPE_SCA', physQuanName, 'GRANDEUR', repk=physQuanScal)
            if (physQuanScal .ne. 'R') then
                call utmess('F', 'TABLE0_42')
            end if
            if (fieldDisc(1:2) .eq. 'EL') then
                call dismoi('NOM_LIGREL', field, 'CHAMP', repk=ligrelField)
                call jeveuo(ligrelField//'.REPE', 'L', vi=repe)
            end if
            if (fieldDisc .eq. 'ELGA') then
                if (lResultUser) then
                    call dismoi('CARA_ELEM', field, 'RESULTAT', repk=caraElem, arret='C', ier=iret)
                    if (iret .eq. 0) then
                        hasCaraElem = ASTER_TRUE
                    end if
                else if (lFieldUser) then
                    call getvid('RESU', 'CARA_ELEM', iocc=1, scal=caraElem, nbret=nbRet)
                    if (nbRet .ne. 0) then
                        hasCaraElem = ASTER_TRUE
                    end if
                end if
            end if
            goto 61
        end if
    end do
61  continue
!

!   RECUPERATION DES NOEUDS,MAILLES
    if (fieldDisc .eq. 'NOEU') then
        motcle(1) = 'NOEUD'
        motcle(2) = 'GROUP_NO'
        motcle(3) = 'MAILLE'
        motcle(4) = 'GROUP_MA'
        typmcl(1) = 'NOEUD'
        typmcl(2) = 'GROUP_NO'
        typmcl(3) = 'MAILLE'
        typmcl(4) = 'GROUP_MA'
        call getvtx('RESU', 'TOUT', iocc=1, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call wkvect(mesnoe, 'V V I', nbNode, jlno)
            do iNode = 1, nbNode
                zi(jlno+iNode-1) = iNode
            end do
        else
            call reliem(' ', mesh, 'NU_NOEUD', 'RESU', 1, &
                        4, motcle, typmcl, mesnoe, nbNode)
            call jeveuo(mesnoe, 'L', jlno)
        end if
        nbCell = 0
!
    else if (fieldDisc(1:2) .eq. 'EL' .or. fieldDisc .eq. 'CART') then
!       VERIFICATIONS
        call getvtx('RESU', 'NOEUD', iocc=1, nbval=0, nbret=n1)
        call getvtx('RESU', 'GROUP_NO', iocc=1, nbval=0, nbret=n2)
        n3 = -n1-n2
        if (n3 .ne. 0) then
            call utmess('F', 'TABLE0_41')
        end if
        motcle(1) = 'MAILLE'
        motcle(2) = 'GROUP_MA'
        typmcl(1) = 'MAILLE'
        typmcl(2) = 'GROUP_MA'
        call getvtx('RESU', 'TOUT', iocc=1, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call wkvect(mesmai, 'V V I', nbCell, jlma)
            if (fieldDisc .eq. 'CART') then
                do iCell = 1, nbCell
                    zi(jlma+iCell-1) = iCell
                end do
            else
!               on ne garde que les mailles du ligrel :
                do iCell = 1, nbCell
                    igrel = repe(1+2*(iCell-1))
                    if (igrel .gt. 0) zi(jlma+iCell-1) = iCell
                end do
            end if
        else
            call reliem(' ', mesh, 'NU_MAILLE', 'RESU', 1, &
                        2, motcle, typmcl, mesmai, nbCell)
        end if
        nbNode = 0
!
        if (fieldDisc .eq. 'ELGA') then
!           calcul de ligrel
            call jeveuo(mesmai, 'L', jlma)
            call jelira(mesmai, 'LONMAX', nbCell)
            call exlim2(zi(jlma), nbCell, ligrelField, 'V', ligrel)
            call mecoor(ligrelField, chgeom)

! --------- Add input field
            lchin(1) = chgeom(1:19)
            lpain(1) = 'PGEOMER'
            nbFieldIn = 1

            if (hasCaraElem) then
! ------------- Add fields for structural elements
                call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! ------------- Add fields for orientation
                call setOrieFields(nbFieldInMax, lpain, lchin, &
                                   nbFieldIn, caraElem)
            end if

! --------- Add output fields
            if (hasCaraElem) then
                lchout(1) = '&&CTDATA.PGCOOR'
                lpaout(1) = 'PCOORPG'
                lchout(2) = '&&CTDATA.SUCOOR'
                lpaout(2) = 'PCOORSU'
                nbFieldOut = 2
!               Champ ELGA aux sous-points
                call cesvar(caraElem, ' ', ligrel, lchout(1))
            else
                lchout(1) = '&&CTDATA.PGCOOR'
                lpaout(1) = 'PCOORPG'
                nbFieldOut = 1
                chpsu = ' '
            end if
!
            call calcul('S', 'COOR_ELGA', ligrel, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        'V', 'OUI')
            call celces(lchout(1), 'V', chpgs)
            if (nbFieldOut .eq. 2) then
!               Si c'est un élément sans sous-point le champ n'est pas calculé
                cel19 = lchout(2) (1:19)
                call exisd('CHAM_ELEM', cel19, ierr)
!               Si le champ existe (ierr=1), il est transformé en CES
                if (ierr .eq. 1) then
                    call celces(lchout(2), 'V', chpsu)
                else
                    chpsu = ' '
                end if
            end if
            call detrsd('LIGREL', ligrel)
        end if
    end if
!
!   RECUPERATION DES COMPOSANTES
    call getvtx('RESU', 'TOUT_CMP', iocc=1, nbval=0, nbret=n1)
    if (n1 .ne. 0) then
        nbcmp = 0
        toucmp = ASTER_TRUE
        call wkvect(nkcmp, 'V V K8', 1, vk8=cmpName)
        cmpName(1) = ' '
    else
        toucmp = .false.
        call getvtx('RESU', 'NOM_CMP', iocc=1, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            nbcmp = -n1
            call wkvect(nkcmp, 'V V K8', nbcmp, vk8=cmpName)
            call getvtx('RESU', 'NOM_CMP', iocc=1, nbval=nbcmp, vect=cmpName)
        else
! --------- Get internal state variables
            call getvtx('RESU', 'NOM_VARI', iocc=1, nbval=0, nbret=nbVari)
            nbVari = -nbVari
            ASSERT(nbVari .gt. 0)
            call wkvect(nkvari, 'V V K16', nbVari, vk16=variName)
            call getvtx('RESU', 'NOM_VARI', iocc=1, nbval=nbVari, vect=variName)
            nbcmp = nbVari
            call wkvect(nkcmp, 'V V K8', nbCell*nbcmp, vk8=cmpName)
            if (result .eq. ' ') then
                call utmess('F', 'EXTRACTION_24')
            end if
            call getvtx('RESU', 'NOM_CHAM', iocc=1, scal=fieldName)
            if (fieldName(1:7) .ne. 'VARI_EL') then
                call utmess('F', 'EXTRACTION_25', sk=fieldName)
            end if
            ASSERT(nbCell .gt. 0)

! --------- Get list of storing index
            call rs_get_liststore(result, nbStore)
            if (nbStore .ne. 0) then
                AS_ALLOCATE(vi=listStore, size=nbStore)
                call rs_get_liststore(result, nbStore, listStore)
            end if

! --------- Get behaviour (only one !)
            call rsGetOneBehaviourFromResult(result, nbStore, listStore, compor)
            if (compor .eq. '#SANS') then
                call utmess('F', 'RESULT1_5')
            end if
            if (compor .eq. '#PLUSIEURS') then
                call utmess('F', 'RESULT1_6')
            end if
            AS_DEALLOCATE(vi=listStore)

! --------- Get name of internal state variables
            if (hasMFront(compor)) then
                call utmess('F', "COMPOR6_6")
            end if
            call varinonu(ligrelField, compor, &
                          nbCell, zi(jlma), &
                          nbVari, variName, cmpName)
        end if
    end if
!
    call jedema()
!
end subroutine
