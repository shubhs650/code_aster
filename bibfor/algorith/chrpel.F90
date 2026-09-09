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
! aslint: disable=W1501
!
subroutine chrpel(fieldOutZ, repereZ, fieldNameZ, iOccField, fieldDimeZ, &
                  model, caraElem, ligrelCalcZ, lModelVariable)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterc/getfac.h"
#include "asterc/r8dgrd.h"
#include "asterf_types.h"
#include "asterfort/angvxy.h"
#include "asterfort/assach.h"
#include "asterfort/assert.h"
#include "asterfort/calc_coor_elga.h"
#include "asterfort/calcul.h"
#include "asterfort/carelo.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/cesvar.h"
#include "asterfort/chrgd.h"
#include "asterfort/chrpan.h"
#include "asterfort/copisd.h"
#include "asterfort/cylrep.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getelem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/matrot.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/normev.h"
#include "asterfort/selectComp.h"
#include "asterfort/sepach.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: fieldOutZ, repereZ, fieldDimeZ
    integer(kind=8), intent(in) :: iOccField
    character(len=*), intent(in) :: fieldNameZ
    character(len=8), intent(in) :: model, caraElem
    character(len=*), intent(in) :: ligrelCalcZ
    aster_logical, intent(in) :: lModelVariable
!
! --------------------------------------------------------------------------------------------------
!
!                       CHANGEMENT DE REPERE DANS LE CAS D'UN CHAM_ELEM
!
! --------------------------------------------------------------------------------------------------
!
!       fieldOutZ      : nom du champ a traiter (champ out)
!       repereZ      : type de repereZ (utilisateur ou cylindrique
!                         ou coque ou coque_util_intr ou coque_intr_util
!                         ou coque_util_cyl)
!       fieldNameZ    : nom de type de champ
!       iOccField       : numero d'occurrence
!       ligrelCalc    : nom du FED servant à restreindre le calcul
!       fieldDimeZ  : type du champ :'tens' 'vect' ou 'coque'
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyw = "AFFE"
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    character(len=8) :: paoutc
    integer(kind=8) :: ii, jj, kk, ino, iad, ipt, isp
    integer(kind=8) :: jcesd, jcesv, jcesl, nbpt
    integer(kind=8) :: ilcnx1, nbsp, iCell
    integer(kind=8) :: nbRet, fieldNbCell, iret
    integer(kind=8) :: ndim, nbCell, cellNume
    integer(kind=8) :: iret0, iret1, nncp
    integer(kind=8) :: ierk, ipaxe, ipaxe2
    integer(kind=8) :: nbpg, nodeNume, ipg, cellNbNode
    integer(kind=8) :: type_pt, ndim_type
    integer(kind=8) :: iocc, nocc
    integer(kind=8) :: iexist, jcesd_gauss, jcesl_gauss, icoo, selectNbCmp
    integer(kind=8), parameter :: type_unknown = 0, type_noeud = 1, type_gauss = 2
!   nb max de points (noeuds|gauss) par élément
    integer(kind=8), parameter :: nptmax = 30
    integer(kind=8), dimension(6) :: permvec
    real(kind=8) :: valr, xnormr, tmp
    real(kind=8), dimension(3) :: xbary, angnot
    real(kind=8), dimension(3) :: orig, axez, vectx, vecty
    real(kind=8), dimension(9) :: valecarte
    real(kind=8), dimension(3, 3) :: pgl, pgcyl, pgu, pglelem
    real(kind=8), dimension(3, nptmax), target :: xno, xpg
    character(len=3) :: physQuanScal
    integer(kind=8), parameter :: selectNbCmpMax = 8
    character(len=8) :: selectCmpName(selectNbCmpMax)
    integer(kind=8), parameter :: nbCmp = 9
    character(len=8), parameter :: cmpName(nbCmp) = (/'ALPHA', 'BETA ', 'REP  ', &
                                                      'AXE_X', 'AXE_Y', 'AXE_Z', &
                                                      'O_X  ', 'O_Y  ', 'O_Z  '/)
    character(len=16) :: option, fieldDime
    character(len=8) :: mesh, answer, physQuanName, fieldDisc
    character(len=19) :: ligrelField, ligrelCalc
    character(len=19), parameter :: canbsp = '&&CHRPEL.NBSP'
    character(len=19), parameter :: chams0 = '&&CHRPEL.CHAMS0', chams1 = '&&CHRPEL.CHAMS1'
    character(len=19) :: carte, chr, chi, ch1, ch2
    character(len=19), parameter :: changl = '&&CHRPEL.ANGL'
    character(len=19) :: celgauss, cesgauss
    character(len=24) :: chgeom
    character(len=24), parameter :: listCellJv = '&&CHRPEL.MES_MAILLES'
    character(len=24) :: valk(3)
    integer(kind=8) :: jcesvrepso(3), jcesdrepso(3), adressev(3)
    integer(kind=8) :: nbptii, nbspii, ncmpii
    character(len=19) :: chrel(3), chres(3)
    integer(kind=8), pointer :: listCell(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    real(kind=8), pointer :: coo_gauss(:) => null()
    real(kind=8), dimension(:, :), pointer :: xpt => null()
    character(len=8), pointer :: cesk(:) => null()
    aster_logical :: exi_cmp, exi_local, lAllCell
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    ipaxe = 0
    fieldDime = fieldDimeZ
    ligrelCalc = ligrelCalcZ
    call dismoi('NOM_LIGREL', fieldOutZ, 'CHAM_ELEM', repk=ligrelField)
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "

!   DEFINITION ET CREATION DU CHAM_ELEM SIMPLE CHAMS1 A PARTIR DU CHAM_ELEM CHAMP1
    call celces(fieldOutZ, 'V', chams0)

! - Select components in field
    call selectComp(chams0, fieldNameZ, fieldDime, selectNbCmp, selectCmpName, ndim_type)
    ASSERT(selectNbCmp .le. selectNbCmpMax)
    call cesred(chams0, 0, [0], selectNbCmp, selectCmpName, 'V', chams1)
    call detrsd('CHAM_ELEM_S', chams0)
    call jeveuo(chams1//'.CESK', 'L', vk8=cesk)
    call jeveuo(chams1//'.CESD', 'L', jcesd)
    mesh = cesk(1)
    physQuanName = cesk(2)
    call dismoi('TYPE_SCA', physQuanName, 'GRANDEUR', repk=physQuanScal)

!   ON EXCLUT LES MOT-CLES 'NOEUD' ET 'GROUP_NO'
    call getvtx(factorKeyw, 'NOEUD', iocc=iOccField, nbval=0, nbret=iret0)
    call jeexin(mesh//'.GROUPENO', ierk)
    if (ierk .ne. 0) then
        call getvtx(factorKeyw, 'GROUP_NO', iocc=iOccField, nbval=0, nbret=iret1)
    else
        iret1 = 0
    end if
    if ((iret0 .lt. 0) .or. (iret1 .lt. 0)) then
        valk(1) = 'NOEUD ou GROUP_NO'
        valk(2) = fieldNameZ
        valk(3) = ' '
        call utmess('F', 'ALGORITH12_42', nk=3, valk=valk)
    end if

!   nombre total de mailles du champ
    fieldNbCell = zi(jcesd-1+1)
!
    ndim = 3
    call dismoi('Z_CST', mesh, 'MAILLAGE', repk=answer)
    if (answer .eq. 'OUI') then
        ndim = 2
    end if
    if (ndim .gt. ndim_type) then
        call utmess('F', 'ALGORITH12_45', sk=fieldDime)
    else if (ndim .lt. ndim_type) then
        ndim = 3
        call utmess('A', 'ALGORITH12_44', sk=fieldDime)
    end if
!
    call jeexin(mesh//'.CONNEX', iret)
    ASSERT(iret .ne. 0)
    call jeveuo(mesh//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    call jeveuo(chams1//'.CESV', 'E', jcesv)
    call jeveuo(chams1//'.CESL', 'L', jcesl)
!
!   Si le champ est exprimé dans le repère local des éléments
!       Construction du champ des repère locaux
    exi_local = ASTER_FALSE
    if (fieldDime .eq. '1D_GENE') then
        chrel(1) = '&&CHRPEL.REPLO_1'
        chrel(2) = '&&CHRPEL.REPLO_2'
        chrel(3) = '&&CHRPEL.REPLO_3'
        chres(1) = '&&CHRPEL.REPSO_1'
        chres(2) = '&&CHRPEL.REPSO_2'
        chres(3) = '&&CHRPEL.REPSO_3'
        call carelo(model, caraElem, 'V', chrel(1), chrel(2), chrel(3))
        exi_local = ASTER_TRUE
!
        do ii = 1, 3
            call celces(chrel(ii), 'V', chres(ii))
            call jeveuo(chres(ii)//'.CESV', 'L', jcesvrepso(ii))
            call jeveuo(chres(ii)//'.CESD', 'L', jcesdrepso(ii))
            call detrsd('CHAM_ELEM', chrel(ii))
        end do
    end if
!   Le mot-clé AFFE définit les caractéristiques du nouveau repère
!   On peut définir un repère variable en définissant ces paramètres par mailles/groupes de mailles
    call getfac(factorKeyw, nocc)
!   Boucle sur les occurrences de AFFE
    do iocc = 1, nocc
! ----- Get list of cells from user
        call getelem(mesh, factorKeyw, iocc, ' ', listCellJv, nbCell)
        lAllCell = ASTER_FALSE
        if (nbCell .gt. 0) then
            call jeveuo(listCellJv, 'L', vi=listCell)
        else
            lAllCell = ASTER_TRUE
            nbCell = fieldNbCell
        end if

        axez = 0.d0; orig = 0.d0; angnot = 0.d0
!       Changement de repère "UTILISATEUR"
        if (repereZ(1:11) .eq. 'UTILISATEUR') then
!           SI LE NOUVEAU REPERE EST DONNE VIA DES VECTEURS
            call getvr8(factorKeyw, 'VECT_X', iocc=iocc, nbval=3, &
                        vect=vectx, nbret=nbRet)
            if (nbRet .ne. 0) then
                call getvr8(factorKeyw, 'VECT_Y', iocc=iocc, nbval=3, &
                            vect=vecty, nbret=nbRet)
                if (ndim .ne. 3) then
                    call utmess('F', 'ALGORITH2_4')
                end if
                call angvxy(vectx, vecty, angnot)
            else
                if (ndim .eq. 3) then
                    call getvr8(factorKeyw, 'ANGL_NAUT', iocc=iocc, nbval=3, &
                                vect=angnot, nbret=nbRet)
                    if (nbRet .ne. 3) then
                        call utmess('F', 'ALGORITH2_7')
                    end if
                else
                    call getvr8(factorKeyw, 'ANGL_NAUT', iocc=iocc, &
                                scal=angnot(1), nbret=nbRet)
                    if (nbRet .ne. 1) then
                        valr = angnot(1)
                        call utmess('A', 'ALGORITH12_43', sr=angnot(1))
                    end if
                end if
                angnot = angnot*r8dgrd()
            end if
!
!           Matrice de passage du repère global vers le repère utilisateur
            call matrot(angnot, pgl)
!           Matrot retourne la transposée de la matrice de passage : on transpose pour avoir
!           la matrice de passage
            pgu = transpose(pgl)
!
!           Appliquer le changement de repère pour les mailles sélectionnées
!
            do iCell = 1, nbCell
                if (lAllCell) then
                    cellNume = iCell
                else
                    cellNume = listCell(iCell)
                end if

!               Si le champ est dans le repère local de l'élément on va chercher la matrice
!               de passage du repére local au global
                if (exi_local) then
!                   Vérification du nombre : de point, de sous points, des composantes
!                   Récupération de l'adresse des valeurs des composantes
                    do ii = 1, 3
                        nbptii = zi(jcesdrepso(ii)-1+5+4*(cellNume-1)+1)
                        nbspii = zi(jcesdrepso(ii)-1+5+4*(cellNume-1)+2)
                        ncmpii = zi(jcesdrepso(ii)-1+5+4*(cellNume-1)+3)
                        ASSERT((nbptii .eq. 1) .and. (nbspii .eq. 1) .and. (ncmpii .eq. 3))
                        adressev(ii) = jcesvrepso(ii)-1+zi(jcesdrepso(ii)-1+5+4*(cellNume-1)+4)
                    end do
!                   La matrice de changement de repère lié à l'élément (optenu par matrot)
                    do ii = 1, 3
                        pglelem(1, ii) = zr(adressev(1)+ii)
                        pglelem(2, ii) = zr(adressev(2)+ii)
                        pglelem(3, ii) = zr(adressev(3)+ii)
                    end do
!                   Passage dans le repère Global           Fglob = transpose(pglelem) . Floc
!                   Passage dans le repère utilisateur      Futil = pgl . Fglob
!                   Au Final                                Futil = pgl . transpose(pglelem) . Floc
                    do ii = 1, 3
                        do jj = 1, 3
                            tmp = 0.0
                            do kk = 1, 3
                                tmp = tmp+pgl(ii, kk)*pglelem(jj, kk)
                            end do
                            pgu(ii, jj) = tmp
                        end do
                    end do
                end if
                nbpt = zi(jcesd-1+5+4*(cellNume-1)+1)
                nbsp = zi(jcesd-1+5+4*(cellNume-1)+2)
                cipt1: do ipt = 1, nbpt
                    do isp = 1, nbsp
                        exi_cmp = ASTER_FALSE
                        do ii = 1, selectNbCmp
                            call cesexi('C', jcesd, jcesl, cellNume, ipt, &
                                        isp, ii, iad)
                            if (iad .gt. 0) then
                                exi_cmp = ASTER_TRUE
                            end if
                        end do
                        if (exi_cmp) then
                            call chrgd(selectNbCmp, jcesd, jcesl, jcesv, cellNume, &
                                       ipt, isp, fieldDime, physQuanScal, pgu)
                        else
                            cycle cipt1
                        end if
                    end do
                end do cipt1
            end do
!
!       Changement de repère "CYLINDRIQUE"
        else if (repereZ(1:11) .eq. 'CYLINDRIQUE') then
!
            if (fieldDime .eq. 'VECTR_3D') then
                call utmess('F', 'ALGORITH2_31')
            end if
!
            call dismoi('TYPE_CHAMP', fieldOutZ, 'CHAMP', repk=fieldDisc, arret='C', ier=iret)
            if (ndim .eq. 3) then
                call getvr8(factorKeyw, 'ORIGINE', iocc=iocc, nbval=3, &
                            vect=orig, nbret=nbRet)
                if (nbRet .ne. 3) then
                    call utmess('F', 'ALGORITH2_8')
                end if
                call getvr8(factorKeyw, 'AXE_Z', iocc=iocc, nbval=3, &
                            vect=axez, nbret=nbRet)
                if (nbRet .eq. 0) then
                    call utmess('F', 'ALGORITH2_9')
                end if
            else
                call getvr8(factorKeyw, 'ORIGINE', iocc=iocc, nbval=2, &
                            vect=orig, nbret=nbRet)
                if (nbRet .ne. 2) then
                    call utmess('A', 'ALGORITH2_10')
                end if
                call getvr8(factorKeyw, 'AXE_Z', iocc=iocc, nbval=0, nbret=nbRet)
                if (nbRet .ne. 0) then
                    call utmess('A', 'ALGORITH2_11')
                end if
                axez(1) = 0.0d0
                axez(2) = 0.0d0
                axez(3) = 1.0d0
            end if
            xnormr = 0.0d0
            call normev(axez, xnormr)
            call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=vale)
!
!           Permutation des composantes en dimension 2
!           Initialisation à l'identité
            permvec(:) = (/(ii, ii=1, 6)/)
            if (ndim == 2) then
                select case (fieldDime(1:4))
                case ('TENS')
                    permvec(4) = 5
                case ('VECT')
                    permvec(2) = 3
                    permvec(3) = 2
                end select
            end if
!
!           Localisation du champ : noeuds/pts de Gauss
            type_pt = type_unknown
            if (fieldDisc(1:4) == 'VECT') then
                type_pt = type_noeud
            end if
            if (fieldDisc(1:4) == 'ELNO') then
                type_pt = type_noeud
            else if (fieldDisc(1:4) == 'ELGA') then
                type_pt = type_gauss
            end if
            ASSERT(type_pt /= type_unknown)
!
!           Si le champ est un champ 'ELGA', on a besoin des
!           coordonnées des points de Gauss dans chaque élément
            if (type_pt == type_gauss) then
                if (lModelVariable) then
                    call utmess('F', 'RESULT4_90')
                end if
!               On utilise calc_coor_elga qui retourne un champ par élément
!               contenant les coordonnées des points de Gauss
                celgauss = '&&CHRPEL.CEL_GAUSS'
                call exisd('CHAMP', celgauss, iexist)
                if (iexist .eq. 0) then
                    call megeom(model, chgeom)
                    call calc_coor_elga(model, ligrelField, chgeom, celgauss, caraElem)
                end if
!               On transforme ce champ en champ simple
                cesgauss = '&&CHRPEL.CES_GAUSS'
                call celces(celgauss, 'V', cesgauss)
                call jeveuo(cesgauss//'.CESD', 'L', jcesd_gauss)
                call jeveuo(cesgauss//'.CESL', 'L', jcesl_gauss)
                call jeveuo(cesgauss//'.CESV', 'L', vr=coo_gauss)
            end if
!
!           Boucle sur les mailles à transformer
            do iCell = 1, nbCell
!               Récupération de cellNume : indice de la maille courante
!                               cellNbNode : nombre de noeuds de la maille courante
                if (lAllCell) then
                    cellNume = iCell
                else
                    cellNume = listCell(iCell)
                end if
                cellNbNode = zi(ilcnx1+cellNume)-zi(ilcnx1-1+cellNume)
!               Quelques caractéristiques du champ simple à transformer sur cette maille :
!               nbpg : nombre de points de Gauss
                nbpg = zi(jcesd-1+5+4*(cellNume-1)+1)
!               nbsp : nombre de sous-points
                nbsp = zi(jcesd-1+5+4*(cellNume-1)+2)
!               selectNbCmp : nombre de composantes
                selectNbCmp = zi(jcesd-1+5+4*(cellNume-1)+3)
!
!               Coordonnées des noeuds de la maille courante
                xno = 0.d0
                do ino = 1, cellNbNode
                    nodeNume = connex(zi(ilcnx1+cellNume-1)+ino-1)
                    xno(1, ino) = vale(1+3*(nodeNume-1)-1+1)
                    xno(2, ino) = vale(1+3*(nodeNume-1)-1+2)
                    if (ndim == 3) then
                        xno(3, ino) = vale(1+3*(nodeNume-1)-1+3)
                    end if
                end do
!
                select case (type_pt)
                case (type_noeud)
!                       Noeuds de la maille
                    nbpt = cellNbNode
                    xpt => xno(:, :)
                case (type_gauss)
!                       Points de Gauss de la maille, dont il faut récupérer les coordonnées
                    do ipg = 1, nbpg
                        do icoo = 1, 3
                            call cesexi('S', jcesd_gauss, jcesl_gauss, cellNume, ipg, &
                                        1, icoo, iad)
                            xpg(icoo, ipg) = coo_gauss(iad)
                        end do
                    end do
                    nbpt = nbpg
                    xpt => xpg(:, :)
                case default
                    nbpt = 0
                    ASSERT(ASTER_FALSE)
                end select
!
!               Boucle sur les points (cette partie est commune aux champs ELNO et ELGA)
                cipt2: do ipt = 1, nbpt
!                   Calcul de la matrice de passage vers le repère cylindrique
                    call cylrep(ndim, xpt(:, ipt), axez, orig, pgcyl, &
                                ipaxe)
!                   Si le point x appartient à l'axe du repère cylindrique
                    if (ipaxe > 0) then
                        call utmess('A', 'ALGORITH2_13')
!                       Calcul de la matrice de passage au centre de gravité de l'élément
                        xbary(:) = sum(xno(:, 1:cellNbNode), dim=2)
                        xbary(:) = xbary(:)/cellNbNode
                        ipaxe2 = 0
                        call cylrep(ndim, xbary, axez, orig, pgcyl, &
                                    ipaxe2)
!                       Si le centre de gravité de l'élément est aussi sur l'axe, on s'arrête
                        if (ipaxe2 > 0) then
                            call utmess('F', 'ALGORITH2_13')
                        end if
                    end if
!                   Boucle sur les sous-points
                    do isp = 1, nbsp
                        exi_cmp = ASTER_TRUE
                        do ii = 1, selectNbCmp
!                           la composante ii du champ existe-t-elle?
                            exi_cmp = ASTER_FALSE
                            call cesexi('S', jcesd, jcesl, cellNume, ipt, &
                                        isp, ii, iad)
                            if (iad .gt. 0) then
                                exi_cmp = ASTER_TRUE
                            end if
                        end do
                        if (exi_cmp) then
                            call chrgd(selectNbCmp, jcesd, jcesl, jcesv, cellNume, &
                                       ipt, isp, fieldDime, physQuanScal, pgcyl, &
                                       permvec)
                        else
                            cycle cipt2
                        end if
                    end do
                end do cipt2
            end do
!
            if (ipaxe .ne. 0) then
                call utmess('A', 'ALGORITH17_22', si=ipaxe)
            end if
        end if
!
        call jeexin(listCellJv, iret)
        if (iret .ne. 0) call jedetr(listCellJv)
    end do
!
    if ((repereZ(1:11) .eq. 'CYLINDRIQUE') .or. (repereZ(1:11) .eq. 'UTILISATEUR')) then
!       Champ simple -> Cham_elem
        call dismoi('NOM_OPTION', fieldOutZ, 'CHAM_ELEM', repk=option)
        call cescel(chams1, ligrelField, option, ' ', 'OUI', &
                    nncp, 'G', fieldOutZ, 'F', nbRet)
        call detrsd('CHAM_ELEM_S', chams1)
        if (exi_local) then
            do ii = 1, 3
                call detrsd('CHAM_ELEM_S', chres(ii))
            end do
        end if
    end if
!
! --------------------------------------------------------------------------------------------------
!   Changement de repère sur une coque
    if ((repereZ(1:5) .eq. 'COQUE') .or. (repereZ(1:15) .eq. 'COQUE_INTR_UTIL') .or. &
        (repereZ(1:15) .eq. 'COQUE_UTIL_INTR') .or. (repereZ(1:14) .eq. 'COQUE_UTIL_CYL')) then
!        Verifier ligrelCalc
        if (ligrelCalcZ(1:1) .eq. ' ') then
            ligrelCalc = ligrelField
        end if
!       Pour l'instant on ne traite pas le cas de plusieurs occurrences du mot-clé AFFE
        if (nocc /= 1) then
            call utmess('F', 'ALGORITH17_23', sk=repereZ, si=nocc)
        end if
!
        call megeom(model, chgeom)
!
        if ((fieldDime(1:10) .eq. 'COQUE_GENE') .and. (repereZ(1:14) .eq. 'COQUE_UTIL_CYL')) then
            call utmess('F', 'ELEMENTS5_55', nk=2, valk=(/'COQUE_UTIL_CYL', 'COQUE_GENE    '/))
        end if
        if (fieldDime(1:10) .eq. 'COQUE_GENE') then
            option = 'REPE_GENE'
!           Nb de paramètres en entrée de l'option
            nbFieldIn = 4
        else if (fieldDime(1:7) .eq. 'TENS_3D') then
            option = 'REPE_TENS'
            nbFieldIn = 5
        else
            call utmess('F', 'ELEMENTS5_53', sk=fieldDime)
        end if
!
!       GENERATION D UN CHAMP D'ANGLES (CARTE CONSTANTE)
        carte = '&&CHRPEL.ANGL_REP'
        valecarte(:) = 0.0d0
!
        if (repereZ .eq. 'COQUE_INTR_UTIL') then
            valecarte(3) = 1.d0
        else if (repereZ .eq. 'COQUE_UTIL_INTR') then
            valecarte(3) = 2.d0
        else if (repereZ .eq. 'COQUE_UTIL_CYL') then
            valecarte(3) = 3.d0
        end if
        call mecact('V', carte, 'MODELE', model, 'CAORIE_R', &
                    ncmp=nbCmp, lnomcmp=cmpName, vr=valecarte)
!
!       CREATION D UN CHAM_ELEM D'ANGLES EN LISANT LES ANGL_REP
        call chrpan(model, carte, option, changl)
!
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom(1:19)
        lpain(2) = 'PANGREP'
        lchin(2) = changl
        nbFieldIn = 2

! ----- Add fields for structural elements
        call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElem)

        nbFieldIn = nbFieldIn+1
        lchin(nbFieldIn) = fieldOutZ
        call dismoi('NOM_GD', fieldOutZ, 'CHAMP', repk=physQuanName)
        call dismoi('TYPE_SCA', physQuanName, 'GRANDEUR', repk=physQuanScal)
!
        if (fieldDime .eq. 'COQUE_GENE') then
            if (fieldNameZ .eq. 'EFGE_ELGA') then
                lpain(nbFieldIn) = 'PEFGAIN'
                lpaout(1) = 'PEFGAOUT'
                if (physQuanScal .eq. 'C') then
                    paoutc = 'PEFGAOUC'
                end if
            else if ((fieldNameZ .eq. 'EFGE_ELNO') .or. (fieldNameZ .eq. 'EGRU_ELNO')) then
                lpain(nbFieldIn) = 'PEFNOIN'
                lpaout(1) = 'PEFNOOUT'
                if (physQuanScal .eq. 'C') then
                    paoutc = 'PEFNOOUC'
                end if
            else if (fieldNameZ .eq. 'DEGE_ELGA') then
                lpain(nbFieldIn) = 'PDGGAIN'
                lpaout(1) = 'PDGGAOUT'
                if (physQuanScal .eq. 'C') then
                    paoutc = 'PDGGAOUC'
                end if
            else if (fieldNameZ .eq. 'DEGE_ELNO') then
                lpain(nbFieldIn) = 'PDGNOIN'
                lpaout(1) = 'PDGNOOUT'
                if (physQuanScal .eq. 'C') then
                    paoutc = 'PDGNOOUC'
                end if
            else if (fieldNameZ .eq. 'SIEF_ELGA') then
                lpain(nbFieldIn) = 'PEFGAIN'
                lpaout(1) = 'PEFGAOUT'
                if (physQuanScal .eq. 'C') then
                    paoutc = 'PEFGAOUC'
                end if
            else
                call utmess('F', 'ELEMENTS5_51', sk=fieldNameZ)
            end if
        else if (fieldDime .eq. 'TENS_3D') then
            if (fieldNameZ .eq. 'SIGM_ELGA') then
                lpain(nbFieldIn) = 'PCOGAIN'
                lpaout(1) = 'PCOGAOUT'
            else if (fieldNameZ .eq. 'SIGM_ELNO') then
                lpain(nbFieldIn) = 'PCONOIN'
                lpaout(1) = 'PCONOOUT'
            else if (fieldNameZ .eq. 'EPSI_ELGA') then
                lpain(nbFieldIn) = 'PDEGAIN'
                lpaout(1) = 'PDEGAOUT'
            else if (fieldNameZ .eq. 'EPSI_ELNO') then
                lpain(nbFieldIn) = 'PDENOIN'
                lpaout(1) = 'PDENOOUT'
            else
                call utmess('F', 'ELEMENTS5_52', sk=fieldNameZ)
            end if
            !
        end if
        call exisd('CHAM_ELEM_S', canbsp, iret1)
        if (iret1 .ne. 1) then
            call dismoi('MXNBSP', fieldOutZ, 'CHAM_ELEM', repi=nbsp)
!
!           SI LE CHAMP A DEJA ETE EXTRAIT IL FAUT APPELER CESVAR AVEC CE CHAMP
            if (nbsp .eq. 1) then
                call cesvar(fieldOutZ(1:19), ' ', ligrelCalc, canbsp)
            else
                call cesvar(caraElem, ' ', ligrelCalc, canbsp)
            end if
        end if
        lchout(1) = chams1
        call copisd('CHAM_ELEM_S', 'V', canbsp, lchout(1))
!
        if (physQuanScal .eq. 'C') then
            chr = '&&CHRPEL.CHR'
            chi = '&&CHRPEL.CHI'
            ch1 = '&&CHRPEL.CH1'
            ch2 = '&&CHRPEL.CH2'
            call sepach(caraElem, lchin(nbFieldIn), 'V', chr, chi)
            lchin(nbFieldIn) = chr
            call calcul('S', option, ligrelCalc, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, ch1, lpaout, &
                        'V', 'OUI')
            lchin(nbFieldIn) = chi
            call calcul('S', option, ligrelCalc, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, ch2, lpaout, &
                        'V', 'OUI')
            call assach(ch1, ch2, 'V', lchout(1), parout=paoutc)
            call detrsd('CHAMP', chr)
            call detrsd('CHAMP', chi)
            call detrsd('CHAMP', ch1)
            call detrsd('CHAMP', ch2)
        elseif (physQuanScal .eq. 'R') then
            call calcul('S', option, ligrelCalc, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        'V', 'OUI')
        else
            ASSERT(ASTER_FALSE)
        end if
        call detrsd('CHAM_ELEM_S', lchout(1))
        call copisd('CHAMP_GD', 'G', lchout(1), fieldOutZ)
    end if
!
    call detrsd('CHAM_ELEM_S', canbsp)
!
    call jedema()
!
end subroutine chrpel
