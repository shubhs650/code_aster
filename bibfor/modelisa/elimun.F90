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
subroutine elimun(mesh, model, zoneKeyword, nbUnilZone, &
                  nbgdcuJv, compcuJv, noponoJv, nolinoJv, &
                  lisnoeJv, poinoeJv, &
                  nbNodeUnil)
!
    implicit none
!
#include "asterfort/exiscp.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/palino.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: mesh, model
    character(len=16), intent(in) :: zoneKeyword
    integer(kind=8), intent(in) :: nbUnilZone
    character(len=24), intent(in) :: nbgdcuJv, compcuJv, noponoJv, nolinoJv
    character(len=24), intent(in) :: lisnoeJv, poinoeJv
    integer(kind=8), intent(inout) :: nbNodeUnil
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Clean list of nodes for LIAISON_UNILATER
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  zoneKeyword      : keyword to read information on zone
! In  nbUnilZone       : number of zones with LIAISON_UNILATER
! IN  NZOCU  : NOMBRE DE ZONES DE LIAISON_UNILATERALE
! IN  NBGDCU : NOM JEVEUX DE LA SD INFOS POINTEURS GRANDEURS
! IN  COMPCU : NOM JEVEUX DE LA SD CONTENANT LES GRANDEURS DU MEMBRE
!              DE GAUCHE
! IN  NOPONO : NOM DE L'OBJET CONTENANT LE VECTEUR D'INDIRECTION
! IN  NOLINO : NOM DE L'OBJET CONTENANT LA LISTE DES NOEUDS
! IN  POINOE : NOM DE L'OBJET CONTENANT LE VECTEUR D'INDIRECTION
!               DES NOEUDS APRES NETTOYAGE
! IN  LISNOE : NOM DE L'OBJET CONTENANT LES NOEUDS APRES NETTOYAGE
! IN  NBNOE  : NOMBRE DE NOEUDS DANS LA LISTE RESULTANTE
!                VAUT NBNOE = NBTOT-NBSUP
! I/O NNOCO  : NOMBRE DE TOTAL DE NOEUDS POUR TOUTES LES OCCURRENCES
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: k8bla = " "
    character(len=24), parameter :: nodeElimJv = '&&ELIMUN.ELIM'
    integer(kind=8) :: jdebut, jdecal, jdecat
    character(len=8) :: cmpName
    integer(kind=8) :: iNode, jNode, iCmp, iUnilZone, iUserElim
    integer(kind=8) :: nodeNume1, nodeNume2
    integer(kind=8) :: nbNode, nbNodeSupp, nbCmpDoesntExist, nbCmp, ntNodeSupp
    integer(kind=8) :: nbNodeDouble, nbNodeElim, nbNodeWithoutCmp, nbUserElim
    integer(kind=8) :: cmpExist(1)
    integer(kind=8), pointer :: nolino(:) => null(), nopono(:) => null(), nbgdcu(:) => null()
    integer(kind=8), pointer :: poinoe(:) => null(), lisnoe(:) => null(), nodeElim(:) => null()
    character(len=8), pointer :: compcu(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    ntNodeSupp = 0

! - Access to datastructures
    call jeveuo(nolinoJv, 'E', vi=nolino)
    call jeveuo(noponoJv, 'L', vi=nopono)
    call jeveuo(nbgdcuJv, 'L', vi=nbgdcu)
    call jeveuo(compcuJv, 'L', vk8=compcu)

! - Create object
    call wkvect(poinoeJv, 'V V I', nbUnilZone+1, vi=poinoe)
    poinoe(1) = 1
!
    do iUnilZone = 1, nbUnilZone
! ----- Access to nodes of zone
        nbNode = nopono(1+iUnilZone)-nopono(iUnilZone)
        jdebut = nopono(iUnilZone)

! ----- ELIMINATION DES PURS DOUBLONS
        nbNodeDouble = 0
        do iNode = 1, nbNode
            nodeNume1 = nolino(1-2+jdebut+iNode)
            if (nodeNume1 .ne. 0) then
                do jNode = iNode+1, nbNode
                    nodeNume2 = nolino(1-2+jdebut+jNode)
                    if ((nodeNume1 .eq. nodeNume2) .and. (nodeNume2 .ne. 0)) then
                        nolino(1-2+jdebut+jNode) = 0
                        nbNodeDouble = nbNodeDouble+1
                    end if
                end do
            end if
        end do

! ----- RECUPERATION INFOS SANS_NOEUD, SANS_GROUP_NO
        call palino(mesh, zoneKeyword, 'SANS_GROUP_NO', 'SANS_NOEUD', iUnilZone, &
                    nodeElimJv)
        call jeveuo(nodeElimJv, 'L', vi=nodeElim)
        nbUserElim = nodeElim(1)

! ----- ELIMINATION DES SANS_GROUP_NO, SANS_NOEUD
        nbNodeElim = 0
        do iUserElim = 1, nbUserElim
            nodeNume1 = nodeElim(iUserElim)
            if (nodeNume1 .ne. 0) then
                do jNode = 1, nbNode
                    nodeNume2 = nolino(1-2+jdebut+jNode)
                    if ((nodeNume1 .eq. nodeNume2) .and. (nodeNume2 .ne. 0)) then
                        nolino(1-2+jdebut+jNode) = 0
                        nbNodeElim = nbNodeElim+1
                    end if
                end do
            end if
        end do

! ----- ELIMINATION DES NOEUDS NE COMPORTANT AUCUNE DES GRANDEURS
        nbCmp = nbgdcu(1+iUnilZone)-nbgdcu(iUnilZone)
        jdecat = nbgdcu(iUnilZone)
        nbNodeWithoutCmp = 0
        do iNode = 1, nbNode
            nodeNume1 = nolino(1-2+jdebut+iNode)
            if (nodeNume1 .ne. 0) then
                nbCmpDoesntExist = 0
                do iCmp = 1, nbCmp
                    cmpName = compcu(jdecat+iCmp-1)
                    call exiscp(cmpName, k8bla, model, 1, 'NUM', &
                                k8bla, [nodeNume1], cmpExist)
                    if (cmpExist(1) .eq. 0) then
                        nbCmpDoesntExist = nbCmpDoesntExist+1
                    end if
                end do
                if (nbCmpDoesntExist .eq. nbCmp) then
                    nolino(1-2+jdebut+iNode) = 0
                    nbNodeWithoutCmp = nbNodeWithoutCmp+1
                end if
            end if
        end do

! ----- NOMBRE DE NOEUDS A SUPPRIMER
        nbNodeSupp = nbNodeDouble+nbNodeElim+nbNodeWithoutCmp
        ntNodeSupp = ntNodeSupp+nbNodeSupp

! ----- MAJ VECTEUR POINTEUR INDIRECT (POINOE)
        poinoe(1+iUnilZone) = poinoe(iUnilZone)+nbNode-nbNodeSupp
        if (nbNode .eq. nbNodeSupp) then
            call utmess('F', 'UNILATER_48')
        end if
    end do

! - CREATION DU VECTEUR RESULTANT
    nbNodeUnil = nbNodeUnil-ntNodeSupp
    call wkvect(lisnoeJv, 'V V I', nbNodeUnil, vi=lisnoe)

! - ELIMINATION EFFECTIVE DES NOEUDS
    jdecal = 0
    do iUnilZone = 1, nbUnilZone
        nbNode = nopono(1+iUnilZone)-nopono(iUnilZone)
        jdebut = nopono(iUnilZone)
        do iNode = 1, nbNode
            nodeNume1 = nolino(1-2+jdebut+iNode)
            if (nodeNume1 .ne. 0) then
                lisnoe(1+jdecal) = nodeNume1
                jdecal = jdecal+1
            end if
        end do
    end do
!
    call jedetr(nodeElimJv)
!
    call jedema()
!
end subroutine
