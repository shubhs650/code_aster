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
subroutine rapo2d(numeDofZ, iocc, listRelaZ, loadZ)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/afrela.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/assvec.h"
#include "asterfort/calcul.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim1.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprel.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/malin1.h"
#include "asterfort/mecact.h"
#include "asterfort/mesomm.h"
#include "asterfort/reajre.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/vemare.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: numeDofZ
    integer(kind=8), intent(in) :: iocc
    character(len=*), intent(in) :: listRelaZ, loadZ
!
! --------------------------------------------------------------------------------------------------
!
! LIAISON_ELEM
!
! For Beam/solid (2D)
!
! --------------------------------------------------------------------------------------------------
!
! In  numeDof          : name of numbering object (NUME_DDL)
! In  iocc             : index of factor keyword
! In  listRela         : name of object for linear relations
! In  load             : load
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = "LIAISON_ELEM"
    integer(kind=8), parameter :: nbCmpMaxi = 320
    character(len=4), parameter :: valeType = "REEL"
    character(len=8), parameter :: physQuanName = 'DEPL_R', geomCmpName(2) = (/'X', 'Y'/)
    character(len=4) :: typcoe
    character(len=8) :: betaf, model
    character(len=8) :: mesh, cmpName(nbCmpMaxi)
    character(len=8) :: beamNodeName
    character(len=8) :: lpain(2), lpaout(2)
    character(len=9) :: nomte
    character(len=16) :: motcle(4), typmcl(4), option
    character(len=19) :: modelLigrel, ligrel
    character(len=24) :: lchin(2), lchout(2), nolili, valk(2)
    character(len=24) :: jvCellName, jvNodeName
    character(len=24) :: vale1, vale2, grnoma
    integer(kind=8) :: ntypel(nbCmpMaxi), dispCmpNume(6), niv, ifm
    integer(kind=8) :: iop, vali(2)
    integer(kind=8) :: cataCmpNameSize, liliMesh, nodeNume
    integer(kind=8) :: idch1, idch2, jno2, iaprno
    integer(kind=8) :: ival
    integer(kind=8) :: iNode, iLili, iCmp
    integer(kind=8) :: nbNodeBeam, beamNodeNume
    integer(kind=8) :: nbNode, nbCell, nbCmp, nbec, nbLili, nbTerm
    real(kind=8) :: coorig(3), valr(9)
    real(kind=8) :: igzz
    real(kind=8) :: beamCoorX, beamCoorY, sectionArea
    real(kind=8) :: ax, ay, axx, ayy, beta, dnorme
    real(kind=8) :: xg, yg
    real(kind=8) :: un, eps
    complex(kind=8) :: betac, ccmp(3)
    complex(kind=8), pointer :: coec(:) => null()
    real(kind=8), pointer :: coer(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    real(kind=8), pointer :: inertie_raccord(:) => null()
    character(len=8), pointer :: lisddl(:) => null()
    character(len=8), pointer :: lisno(:) => null()
    real(kind=8), pointer :: nodeCoor(:) => null()
    character(len=8) :: load
    character(len=14) :: numeDof
    character(len=19) :: listRela
    integer(kind=8), pointer :: cellNume(:) => null()
    character(len=8), pointer :: nodeName(:) => null(), cataCmpName(:) => null()
    character(len=24), parameter :: mapSection = '&&RAPO2D.CAORIGE'
    character(len=24), parameter :: fieldSection = '&&RAPO2D.PSECT'
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
    numeDof = numeDofZ
    load = loadZ
    listRela = listRelaZ
!
    call getvtx(factorKeyword, 'OPTION', iocc=iocc, scal=option, nbret=iop)

! - INITIALISATIONS
    typcoe = 'REEL'
    betaf = '&FOZERO'
    beta = 0.0d0
    betac = (0.0d0, 0.0d0)
    eps = 1.d-2
    un = 1.0d0
    ccmp = (0.0d0, 0.0d0)
    dispCmpNume = 0

!
    ligrel = '&&RAPO2D'
    jvNodeName = '&&RAPO2D.LISTE_NOEUDS'
    jvCellName = '&&RAPO2D.LISTE_MAILLES'
    motcle(1) = 'GROUP_MA_1'
    motcle(2) = 'MAILLE_1'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'

! - Main parameters
    call dismoi('NOM_MODELE', load, 'CHARGE', repk=model)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NOM_MAILLA', modelLigrel, 'LIGREL', repk=mesh)
    grnoma = mesh//'.GROUPENO'
    call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=nodeCoor)

! - Get components
    nomte = 'DDEPL_R_'
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', physQuanName), 'L', vk8=cataCmpName)
    call jelira(jexnom('&CATA.GD.NOMCMP', physQuanName), 'LONMAX', cataCmpNameSize)
    nbCmp = cataCmpNameSize-1
    ASSERT(nbCmp .le. nbCmpMaxi)
    do iCmp = 1, nbCmp
        cmpName(iCmp) = cataCmpName(iCmp)
        call jenonu(jexnom('&CATA.TE.NOMTE', nomte//cmpName(iCmp) (1:7)), ntypel(iCmp))
    end do
    call dismoi('NB_EC', physQuanName, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 11)

! - Access to DOF numbering
    call jelira(numeDof//'.NUME.PRNO', 'NMAXOC', nbLili)
    liliMesh = 0
    do iLili = 1, nbLili
        call jenuno(jexnum(numeDof//'.NUME.LILI', iLili), nolili)
        if (nolili(1:8) .ne. '&MAILLA') goto 30
        liliMesh = iLili
30      continue
    end do
    ASSERT(liliMesh .ne. 0)
    call jeveuo(jexnum(numeDof//'.NUME.PRNO', liliMesh), 'L', iaprno)

! - Get list of nodes
    call malin1(factorKeyword, load, iocc, 1, jvNodeName, nbNode)
    call jeveuo(jvNodeName, 'L', vk8=nodeName)

! - Get list of cells
    call reliem(' ', mesh, 'NU_MAILLE', factorKeyword, iocc, &
                2, motcle(1), typmcl(1), jvCellName, nbCell)
    call jeveuo(jvCellName, 'L', vi=cellNume)

! - Create reduced list of elements
    call exlim1(cellNume, nbCell, model, 'V', ligrel)

! - Get node for beam
    motcle(1) = 'GROUP_NO_2'
    motcle(2) = 'NOEUD_2'
    motcle(3) = 'GROUP_MA_2'
    motcle(4) = 'MAILLE_2'
    typmcl(1) = 'GROUP_NO'
    typmcl(2) = 'NOEUD'
    typmcl(3) = 'GROUP_MA'
    typmcl(4) = 'MAILLE'
    call reliem(' ', mesh, 'NO_NOEUD', factorKeyword, iocc, &
                4, motcle, typmcl, '&&RAPO2D.NO2', nbNodeBeam)
    if (nbNodeBeam .ne. 1) then
        call utmess('F', 'CHARGES10_13')
    end if
    call jeveuo('&&RAPO2D.NO2', 'L', jno2)
    beamNodeName = zk8(jno2)
    beamNodeNume = char8_to_int(beamNodeName)
    call jedetr('&&RAPO2D.NO2')

! - COORDONNEES DU NOEUD POUTRE
    beamCoorX = nodeCoor(3*(beamNodeNume-1)+1)
    beamCoorY = nodeCoor(3*(beamNodeNume-1)+2)

! - CALCUL SUR CHAQUE ELEMENT DE BORD A RELIER A LA POUTRE
! - DES CARACTERISTIQUES GEOMETRIQUES SUIVANTES :
! - SOMME/B_ELEMENT(1,X,Y,X*X,Y*Y,X*Y)DS
    lpain(1) = 'PGEOMER'
    lchin(1) = mesh//'.COORDO'
    lpaout(1) = 'PCASECT'
    lchout(1) = fieldSection
!
    call calcul('S', 'CARA_SECT_POUT3', ligrel, 1, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')

! - VECTEUR DES QUANTITES GEOMETRIQUES PRECITEES SOMMEES
! - SUR LA SURFACE DE RACCORD, CES QUANTITES SERONT NOTEES :
! - A1 = S,AX,AY,AXX,AYY
    AS_ALLOCATE(vr=inertie_raccord, size=6)

! - SOMMATION DES QUANTITES GEOMETRIQUES ELEMENTAIRES
    call mesomm(fieldSection, 6, vr=inertie_raccord)
    sectionArea = inertie_raccord(1)
    ax = inertie_raccord(2)
    ay = inertie_raccord(3)
    axx = inertie_raccord(4)
    ayy = inertie_raccord(5)
    if (abs(sectionArea) .lt. r8prem()) then
        call utmess('F', 'CHARGES10_11', sk=option)
    end if

! - COORDONNEES DU CENTRE GEOMETRIQUE G DE LA SECTION DE RACCORD
    xg = (1.d0/sectionArea)*ax
    yg = (1.d0/sectionArea)*ay

! - VERIFICATION DE L'IDENTITE GEOMETRIQUE DE G AVEC LE NOEUD POUTRE A RACCORDER
    dnorme = sqrt((beamCoorX-xg)*(beamCoorX-xg)+ &
                  (beamCoorY-yg)*(beamCoorY-yg))/(sectionArea/2.)
    if (dnorme .gt. eps) then
        valr(1) = xg
        valr(2) = yg
        valr(3) = 0.d0
        valr(4) = beamCoorX
        valr(5) = beamCoorY
        valr(6) = 0.d0
        valr(7) = eps*100.0d0
        valr(8) = sectionArea
        valr(9) = dnorme
        valk(1) = option
        vali(1) = iocc
        call utmess('A', 'CHARGES10_12', sk=valk(1), si=vali(1), nr=9, valr=valr)
    end if

! - CALCUL DU TENSEUR D'INERTIE EN G
    igzz = axx+ayy-sectionArea*(xg*xg+yg*yg)

! - COORDONNEES DU CENTRE GEOMETRIQUE G DE LA SECTION DE RACCORD
    coorig(1) = xg
    coorig(2) = yg
!
    call mecact('V', mapSection, 'LIGREL', ligrel, 'GEOM_R', &
                ncmp=2, lnomcmp=geomCmpName, vr=coorig)

! - DETERMINATION DE 2 LISTES  DE VECTEURS PAR ELEMENT PRENANT
! - LEURS VALEURS AUX NOEUDS DES ELEMENTS.
! - LA PREMIERE LISTE DE NOM 'VECT_NI' A POUR VALEURS AU NOEUD
! - I D'UN ELEMENT :
! - SOMME/S_ELEMENT(NI,0)DS
! - LA SECONDE LISTE DE NOM 'VECT_XYNI' A POUR VALEURS AU NOEUD
! - I D'UN ELEMENT : SOMME/S_ELEMENT(X*NI,Y*NI)DS
! - AVEC X = XM - XG = NJ*XJ - XG
! -      Y = YM - YG = NJ*YJ - YG
    lpain(1) = 'PGEOMER'
    lchin(1) = mesh//'.COORDO'
    lpain(2) = 'PORIGIN'
    lchin(2) = mapSection
    lpaout(1) = 'PVECTU1'
    lpaout(2) = 'PVECTU2'
    lchout(1) = '&&RAPO2D.VECT_NI'
    lchout(2) = '&&RAPO2D.VECT_XYZNI'
!
    call calcul('S', 'CARA_SECT_POUT4', ligrel, 2, lchin, &
                lpain, 2, lchout, lpaout, 'V', &
                'OUI')

! - CREATION DES .RERR DES VECTEURS EN SORTIE DE CALCUL
    call vemare('V', '&&RAPO2D', model)

! - ASSEMBLAGE DE LCHOUT(1) DANS LE CHAMNO DE NOM 'CH_DEPL_1'
    call jedetr('&&RAPO2D           .RELR')
    call reajre('&&RAPO2D', lchout(1), 'V')
    call assvec('V', 'CH_DEPL_1', 1, '&&RAPO2D           .RELR', [1.d0], numeDof)

! - ASSEMBLAGE DE LCHOUT(2) DANS LE CHAMNO DE NOM 'CH_DEPL_2'
    call jedetr('&&RAPO2D           .RELR')
    call reajre('&&RAPO2D', lchout(2), 'V')
    call assvec('V', 'CH_DEPL_2', 1, '&&RAPO2D           .RELR', [1.d0], numeDof)
    vale1 = 'CH_DEPL_1          .VALE'
    vale2 = 'CH_DEPL_2          .VALE'
    call jeveuo(vale1, 'L', idch1)
    call jeveuo(vale2, 'L', idch2)

! - CREATION DES TABLEAUX NECESSAIRES A L'AFFECTATION DE LISREL
! - MAJORANT DU NOMBRE DE TERMES DANS UNE RELATION
    nbterm = 2*nbNode+2
    AS_ALLOCATE(vk8=lisno, size=nbterm)
    AS_ALLOCATE(vk8=lisddl, size=nbterm)
    AS_ALLOCATE(vr=coer, size=nbterm)
    AS_ALLOCATE(vc=coec, size=nbterm)
    AS_ALLOCATE(vr=direct, size=2*nbterm)
    AS_ALLOCATE(vi=dime, size=nbterm)

! =================================================================================================
! - First group: SOMME/S_RACCORD(U_2D) = S_RACCORD*U_NOEUD_POUTRE
! =================================================================================================

! - First relation: -S.DX(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DX(NOEUD_I) = 0
    nbterm = nbNode+1
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
        ival = zi(iaprno+(nodeNume-1)*(nbec+2))
        lisno(iNode) = nodeName(iNode)
        lisddl(iNode) = 'DX'
        coer(iNode) = zr(idch1+ival-1)
    end do
!
    lisno(1+nbNode+1-1) = beamNodeName
    lisddl(1+nbNode+1-1) = 'DX'
    coer(1+nbNode+1-1) = -sectionArea
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, valeType, 0.d0, listRela)
    call imprel(factorKeyword, nbterm, coer, lisddl, lisno, &
                beta)

! - Second relation: -S.DY(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DY(NOEUD_I) = 0
    nbterm = nbNode+1
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
        ival = zi(iaprno+(nodeNume-1)*(nbec+2))
        lisno(iNode) = nodeName(iNode)
        lisddl(iNode) = 'DY'
        coer(iNode) = zr(idch1+ival-1)
    end do
!
    lisno(1+nbNode+1-1) = beamNodeName
    lisddl(1+nbNode+1-1) = 'DY'
    coer(1+nbNode+1-1) = -sectionArea
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, valeType, 0.d0, listRela)
    call imprel(factorKeyword, nbterm, coer, lisddl, lisno, &
                beta)

! =================================================================================================
! - Second group: SOMME/S_RACCORD(GM X U_2D) = I.OMEGA(NOEUD_POUTRE)
! =================================================================================================
! - Third relation:
!        (SOMME/S_RACCORD(X*NI.DS)).DY(NOEUD_I) -
!        (SOMME/S_RACCORD(Y*NI.DS)).DX(NOEUD_I) -
!        IZZ.DRZ(NOEUD_POUTRE) = 0
    nbterm = 2*nbNode+1
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
        ival = zi(iaprno+(nodeNume-1)*(nbec+2))
        lisno(1+2*(iNode-1)) = nodeName(iNode)
        lisno(1+2*(iNode-1)+1) = nodeName(iNode)
        lisddl(1+2*(iNode-1)) = 'DY'
        lisddl(1+2*(iNode-1)+1) = 'DX'
        coer(1+2*(iNode-1)) = zr(idch2+ival-1)
        coer(1+2*(iNode-1)+1) = -zr(idch2+ival-1+1)
    end do
!
    lisno(1+2*nbNode+1-1) = beamNodeName
    lisddl(1+2*nbNode+1-1) = 'DRZ'
    coer(1+2*nbNode+1-1) = -igzz
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, valeType, 0.d0, listRela)
    call imprel(factorKeyword, nbterm, coer, lisddl, lisno, &
                beta)

! - Clean
    call detrsd('LIGREL', ligrel)
    call jedetr('&&RAPO2D.LISTE_NOEUDS')
    call jedetr('&&RAPO2D.LISTE_MAILLES')
    call detrsd('CHAMP_GD', fieldSection)
    AS_DEALLOCATE(vr=inertie_raccord)
    call detrsd('CARTE', mapSection)
    call detrsd('RESUELEM', '&&RAPO2D.VECT_NI')
    call detrsd('RESUELEM', '&&RAPO2D.VECT_XYZNI')
    call jedetr('&&RAPO2D           .RELR')
    call jedetr('&&RAPO2D           .RERR')
    AS_DEALLOCATE(vk8=lisno)
    AS_DEALLOCATE(vk8=lisddl)
    AS_DEALLOCATE(vr=coer)
    AS_DEALLOCATE(vc=coec)
    AS_DEALLOCATE(vr=direct)
    AS_DEALLOCATE(vi=dime)
    call detrsd('CHAMP_GD', 'CH_DEPL_1')
    call detrsd('CHAMP_GD', 'CH_DEPL_2')
!
    call jedema()
end subroutine
