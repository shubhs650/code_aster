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
subroutine rapoco(numeDofZ, iocc, listRelaZ, loadZ)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterc/indik8.h"
#include "asterc/r8pi.h"
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
#include "asterfort/exisdg.h"
#include "asterfort/exlim1.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprel.h"
#include "asterfort/infniv.h"
#include "asterfort/int_to_char8.h"
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
#include "asterfort/racotu.h"
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
! For Beam/shell and beam/pipe
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
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    integer(kind=8) :: nbFieldIn, nbFieldOut
    character(len=16), parameter :: factorKeyword = "LIAISON_ELEM"
    integer(kind=8), parameter :: nbCmpMaxi = 330
    character(len=4), parameter :: valeType = "REEL"
    character(len=8), parameter :: physQuanName = 'DEPL_R', geomCmpName(3) = (/'X', 'Y', 'Z'/)
    character(len=8), parameter :: dispCmpName(6) = (/'DX ', 'DY ', 'DZ ', &
                                                      'DRX', 'DRY', 'DRZ'/)
    character(len=4) :: typcoe
    character(len=8) :: betaf, model, k8bid, caraElem
    character(len=8) :: mesh, cmpName(nbCmpMaxi)
    character(len=8) :: beamNodeName
    character(len=9) :: nomte
    character(len=16) :: motcle(2), typmcl(2), option
    character(len=19) :: modelLigrel, ligrel
    character(len=24) :: nolili, valk(2)
    character(len=24) :: jvCellName, jvNodeName
    character(len=24) :: vale1, vale2, grnoma
    character(len=24) :: nodeGroupName
    integer(kind=8) :: ntypel(nbCmpMaxi), dispCmpNume(6), niv, ifm
    integer(kind=8) :: iop, vali(2)
    integer(kind=8) :: cataCmpNameSize, liliMesh, nodeNume
    integer(kind=8) :: idch1, idch2, iaprno, nbval
    integer(kind=8) :: naxe, ival, dg
    integer(kind=8) :: geomDime, ncara
    integer(kind=8) :: nbgno, jgro
    integer(kind=8) :: iNode, iLili, iCmp
    integer(kind=8) :: nbNodeBeam, beamNodeNume
    integer(kind=8) :: nbNode, nbCell, nbCmp, nbec, nbLili, nbTerm
    real(kind=8) :: coorig(3), valr(9)
    real(kind=8) :: beamAxis(3)
    real(kind=8) :: ig(6)
    real(kind=8) :: beamCoorX, beamCoorY, beamCoorZ, xnorm, sectionArea
    real(kind=8) :: ax, ay, az, axx, ayy, azz, axy, axz, ayz, beta, dnorme
    real(kind=8) :: xg, yg, zg
    real(kind=8) :: un, pi, eps
    complex(kind=8) :: betac, ccmp(3)
    complex(kind=8), pointer :: coec(:) => null()
    real(kind=8), pointer :: coer(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    real(kind=8), pointer :: inertie_raccord(:) => null()
    character(len=8), pointer :: lisddl(:) => null()
    character(len=8), pointer :: lisno(:) => null()
    real(kind=8), pointer :: nodeCoor(:) => null()
    integer(kind=8), pointer :: prnm(:) => null()
    character(len=8) :: load
    character(len=14) :: numeDof
    character(len=19) :: listRela
    integer(kind=8), pointer :: cellNume(:) => null()
    character(len=8), pointer :: nodeName(:) => null(), cataCmpName(:) => null()
    character(len=24), parameter :: mapBeamAxis = '&&RAPOCO.CAXE_POU'
    character(len=24), parameter :: mapSection = '&&RAPOCO.CAORIGE'
    character(len=24), parameter :: fieldSection = '&&RAPOCO.PSECT'
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
    numeDof = numeDofZ
    load = loadZ
    listRela = listRelaZ
    pi = r8pi()
!
    call getvtx(factorKeyword, 'OPTION', iocc=iocc, scal=option, nbret=iop)

! - Initialiations
    typcoe = 'REEL'
    betaf = '&FOZERO'
    beta = 0.0d0
    betac = (0.0d0, 0.0d0)
    eps = 1.d-2
    un = 1.0d0
    ccmp = (0.0d0, 0.0d0)
    dispCmpNume = 0
    lchin = ' '
    lpain = ' '
    lchout = ' '
    lpaout = ' '

!
    ligrel = '&&RAPOCO'
    jvNodeName = '&&RAPOCO.LISTE_NOEUDS'
    jvCellName = '&&RAPOCO.LISTE_MAILLES'
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

! - No super element !
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)
    if (geomDime .ge. 1000) then
        call utmess("F", "CHARGES10_1", sk=option)
    end if

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
    call jeveuo(modelLigrel//'.PRNM', 'L', vi=prnm)
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
    nbgno = 0
    call getvem(mesh, 'GROUP_NO', factorKeyword, 'GROUP_NO_2', iocc, 0, k8bid, nbgno)
    if (nbgno .eq. 0) then
        valk(1) = factorKeyword
        valk(2) = option
        call utmess('F', 'CHARGES10_2', nk=2, valk=valk)
    end if
    if (nbgno .ne. 0) then
        nbgno = -nbgno
        if (nbgno .ne. 1) then
            call utmess('F', 'CHARGES10_3')
        end if
        call getvem(mesh, 'GROUP_NO', factorKeyword, 'GROUP_NO_2', iocc, &
                    nbgno, nodeGroupName, nbval)
        call jelira(jexnom(grnoma, nodeGroupName), 'LONUTI', nbNodeBeam)
        if (nbNodeBeam .ne. 1) then
            call utmess('F', 'CHARGES10_4', sk=nodeGroupName)
        else
            call jeveuo(jexnom(grnoma, nodeGroupName), 'L', jgro)
            beamNodeNume = zi(jgro+1-1)
            beamNodeName = int_to_char8(beamNodeNume)
        end if
    end if

! - Get vector for beam axis
    call getvr8(factorKeyword, 'AXE_POUTRE', iocc=iocc, nbval=3, vect=beamAxis, nbret=naxe)
    if (naxe .eq. 0) then
        call utmess('F', 'CHARGES10_5', sk="AXE_POUTRE")
    end if
    xnorm = sqrt(beamAxis(1)*beamAxis(1)+beamAxis(2)*beamAxis(2)+beamAxis(3)*beamAxis(3))
    if (xnorm .le. r8prem()) then
        call utmess('F', 'CHARGES10_6', sk="AXE_POUTRE")
    end if
    beamAxis(1) = beamAxis(1)/xnorm
    beamAxis(2) = beamAxis(2)/xnorm
    beamAxis(3) = beamAxis(3)/xnorm

! - Set map for beam axis
    call mecact('V', mapBeamAxis, 'LIGREL', ligrel, 'GEOM_R', &
                ncmp=3, lnomcmp=geomCmpName, vr=beamAxis)

! - RECUPERATION DES CARACTERISTIQUES ELEMENTAIRES
    call getvid(factorKeyword, 'CARA_ELEM', iocc=iocc, scal=caraElem, nbret=ncara)
    if (ncara .eq. 0) then
        call utmess('F', 'CHARGES10_7')
    end if

! - COORDONNEES DU NOEUD POUTRE
    beamCoorX = nodeCoor(3*(beamNodeNume-1)+1)
    beamCoorY = nodeCoor(3*(beamNodeNume-1)+2)
    beamCoorZ = nodeCoor(3*(beamNodeNume-1)+3)

! - Check rotation dof for pipe/shell elements
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
        dg = prnm((nodeNume-1)*nbec+1)
        do iCmp = 4, 6
            dispCmpNume(iCmp) = indik8(cmpName, dispCmpName(iCmp), 1, nbCmp)
            if (.not. exisdg([dg], dispCmpNume(iCmp))) then
                call utmess('F', 'CHARGES10_9', sk=dispCmpName(iCmp))
            end if
        end do
    end do

! - Check all dofs for beam element
    dg = prnm((beamNodeNume-1)*nbec+1)
    do iCmp = 1, 6
        dispCmpNume(iCmp) = indik8(cmpName, dispCmpName(iCmp), 1, nbCmp)
        if (.not. exisdg([dg], dispCmpNume(iCmp))) then
            call utmess('F', 'CHARGES10_10', sk=dispCmpName(iCmp))
        end if
    end do

! - CALCUL SUR CHAQUE ELEMENT DE BORD A RELIER A LA POUTRE
! - DES CARACTERISTIQUES GEOMETRIQUES SUIVANTES :
! - SOMME/S_ELEMENT(1,X,Y,Z,X*X,Y*Y,Z*Z,X*Y,X*Z,Y*Z)DS

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = mesh//'.COORDO'
    nbFieldIn = 1

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem, caorienZ_=mapBeamAxis)

! - Set output field
    lpaout(1) = 'PCASECT'
    lchout(1) = fieldSection(1:19)
    nbFieldOut = 1
!
    call calcul('S', 'CARA_SECT_POUT3', ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')

! - VECTEUR DES QUANTITES GEOMETRIQUES PRECITEES SOMMEES
! - SUR LA SURFACE DE RACCORD, CES QUANTITES SERONT NOTEES :
! - A1 = S,AX,AY,AZ,AXX,AYY,AZZ,AXY,AXZ,AYZ
    AS_ALLOCATE(vr=inertie_raccord, size=10)

! - SOMMATION DES QUANTITES GEOMETRIQUES ELEMENTAIRES
    call mesomm(fieldSection, 10, vr=inertie_raccord)
    sectionArea = inertie_raccord(1)
    ax = inertie_raccord(2)
    ay = inertie_raccord(3)
    az = inertie_raccord(4)
    axx = inertie_raccord(5)
    ayy = inertie_raccord(6)
    azz = inertie_raccord(7)
    axy = inertie_raccord(8)
    axz = inertie_raccord(9)
    ayz = inertie_raccord(10)
    if (abs(sectionArea) .lt. r8prem()) then
        call utmess('F', 'CHARGES10_11', sk=option)
    end if

! - COORDONNEES DU CENTRE GEOMETRIQUE G DE LA SECTION DE RACCORD
    xg = (1.d0/sectionArea)*ax
    yg = (1.d0/sectionArea)*ay
    zg = (1.d0/sectionArea)*az

! - VERIFICATION DE L'IDENTITE GEOMETRIQUE DE G AVEC LE NOEUD POUTRE A RACCORDER
    dnorme = sqrt((beamCoorX-xg)*(beamCoorX-xg)+ &
                  (beamCoorY-yg)*(beamCoorY-yg)+ &
                  (beamCoorZ-zg)*(beamCoorZ-zg))/sqrt(sectionArea/pi)
    if (dnorme .gt. eps) then
        valr(1) = xg
        valr(2) = yg
        valr(3) = zg
        valr(4) = beamCoorX
        valr(5) = beamCoorY
        valr(6) = beamCoorZ
        valr(7) = eps*100.0d0
        valr(8) = sqrt(sectionArea/pi)
        valr(9) = dnorme
        valk(1) = option
        vali(1) = iocc
        call utmess('A', 'CHARGES10_12', sk=valk(1), si=vali(1), nr=9, valr=valr)
    end if

! - CALCUL DU TENSEUR D'INERTIE EN G
    ig(1) = ayy+azz-sectionArea*(yg*yg+zg*zg)
    ig(2) = -axy+sectionArea*xg*yg
    ig(3) = -axz+sectionArea*xg*zg
    ig(4) = azz+axx-sectionArea*(zg*zg+xg*xg)
    ig(5) = -ayz+sectionArea*yg*zg
    ig(6) = axx+ayy-sectionArea*(xg*xg+yg*yg)

! - COORDONNEES DU CENTRE GEOMETRIQUE G DE LA SECTION DE RACCORD
    coorig(1) = xg
    coorig(2) = yg
    coorig(3) = zg
!
    call mecact('V', mapSection, 'LIGREL', ligrel, 'GEOM_R', &
                ncmp=3, lnomcmp=geomCmpName, vr=coorig)

! - DETERMINATION DE 2 LISTES  DE VECTEURS PAR ELEMENT PRENANT
! - LEURS VALEURS AUX NOEUDS DES ELEMENTS.
! - LA PREMIERE LISTE DE NOM 'VECT_EINI' A POUR VALEURS AU NOEUD
! - I D'UN ELEMENT :
! - SOMME/S_ELEMENT(E1(1)*NI,E1(2)*NI,E1(3)*NI,
! -                 E2(1)*NI,E2(2)*NI,E2(3)*NI)DS
! - OU E1 EST UN VECTEUR UNITAIRE PERPENDICULAIRE A L'ELEMENT
! - DE BORD ORIENTE DE LA COQUE VERS LA POUTRE ET
! - E2 EST LE VECTEUR TANGENT A LA FIBRE MOYENNE DE L'ELEMENT DE BORD
! - LA SECONDE LISTE DE NOM 'VECT_XYZNI' A POUR VALEURS AU NOEUD
! - I D'UN ELEMENT :
! - SOMME/S_ELEMENT(X*NI,Y*NI,Z*NI,NI,0,0)DS
! - AVEC X = XM - XG = NJ*XJ - XG
! -      Y = YM - YG = NJ*YJ - YG
! -      Z = ZM - ZG = NJ*ZJ - ZG

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = mesh//'.COORDO'
    lpain(2) = 'PORIGIN'
    lchin(2) = mapSection(1:19)
    nbFieldIn = 2

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem, caorienZ_=mapBeamAxis)

! - Set output fields
    lpaout(1) = 'PVECTU1'
    lpaout(2) = 'PVECTU2'
    lchout(1) = '&&RAPOCO.VECT_XYZNI'
    lchout(2) = '&&RAPOCO.VECT2'
    nbFieldOut = 2
!
    call calcul('S', 'CARA_SECT_POUT4', ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')

! - CREATION DES .RERR DES VECTEURS EN SORTIE DE CALCUL
    call vemare('V', '&&RAPOCO', model)

! - ASSEMBLAGE DE LCHOUT(1) DANS LE CHAMNO DE NOM 'CH_DEPL_1'
    call jedetr('&&RAPOCO           .RELR')
    call reajre('&&RAPOCO', lchout(1), 'V')
    call assvec('V', 'CH_DEPL_1', 1, '&&RAPOCO           .RELR', [1.d0], numeDof)

! - ASSEMBLAGE DE LCHOUT(2) DANS LE CHAMNO DE NOM 'CH_DEPL_2'
    call jedetr('&&RAPOCO           .RELR')
    call reajre('&&RAPOCO', lchout(2), 'V')
    call assvec('V', 'CH_DEPL_2', 1, '&&RAPOCO           .RELR', [1.d0], numeDof)
    vale1 = 'CH_DEPL_1          .VALE'
    vale2 = 'CH_DEPL_2          .VALE'
    call jeveuo(vale1, 'L', idch1)
    call jeveuo(vale2, 'L', idch2)

! - CREATION DES TABLEAUX NECESSAIRES A L'AFFECTATION DE LISREL
! - MAJORANT DU NOMBRE DE TERMES DANS UNE RELATION
    nbterm = 5*nbNode+3
    AS_ALLOCATE(vk8=lisno, size=nbterm)
    AS_ALLOCATE(vk8=lisddl, size=nbterm)
    AS_ALLOCATE(vr=coer, size=nbterm)
    AS_ALLOCATE(vc=coec, size=nbterm)
    AS_ALLOCATE(vr=direct, size=6*nbterm)
    AS_ALLOCATE(vi=dime, size=nbterm)

! =================================================================================================
! - First group: SOMME/S_RACCORD(U_COQUE) = S_RACCORD*U_NOEUD_POUTRE
! =================================================================================================

! - First relation: -S.DX(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DX(NOEUD_I) = 0
    nbterm = nbNode+1
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
        ival = zi(iaprno+(nodeNume-1)*(nbec+2)+1-1)-1
        lisno(iNode) = nodeName(iNode)
        lisddl(iNode) = 'DX'
        coer(iNode) = zr(idch1-1+ival+4)
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
        ival = zi(iaprno+(nodeNume-1)*(nbec+2)+1-1)-1
        lisno(iNode) = nodeName(iNode)
        lisddl(iNode) = 'DY'
        coer(iNode) = zr(idch1-1+ival+4)
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

! - Third relation: -S.DZ(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DZ(NOEUD_I) = 0
    nbterm = nbNode+1
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
        ival = zi(iaprno+(nodeNume-1)*(nbec+2)+1-1)-1
        lisno(iNode) = nodeName(iNode)
        lisddl(iNode) = 'DZ'
        coer(iNode) = zr(idch1-1+ival+4)
    end do
!
    lisno(1+nbNode+1-1) = beamNodeName
    lisddl(1+nbNode+1-1) = 'DZ'
    coer(1+nbNode+1-1) = -sectionArea
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, valeType, 0.d0, listRela)
    call imprel(factorKeyword, nbterm, coer, lisddl, lisno, &
                beta)

! =================================================================================================
! - Second group: SOMME/S_RACCORD(GM X U_COQUE) = I.OMEGA(NOEUD_POUTRE)
! =================================================================================================
    nbterm = 5*nbNode+3
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
        ival = zi(iaprno+(nodeNume-1)*(nbec+2)+1-1)-1
        lisno(1+5*(iNode-1)+1-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+2-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+3-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+4-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+5-1) = nodeName(iNode)

        lisddl(1+5*(iNode-1)+1-1) = 'DZ'
        lisddl(1+5*(iNode-1)+2-1) = 'DY'
        lisddl(1+5*(iNode-1)+3-1) = 'DRX'
        lisddl(1+5*(iNode-1)+4-1) = 'DRY'
        lisddl(1+5*(iNode-1)+5-1) = 'DRZ'
        coer(1+5*(iNode-1)+1-1) = zr(idch1-1+ival+2)
        coer(1+5*(iNode-1)+2-1) = -zr(idch1-1+ival+3)
        coer(1+5*(iNode-1)+3-1) = zr(idch2-1+ival+1)
        coer(1+5*(iNode-1)+4-1) = zr(idch2-1+ival+2)
        coer(1+5*(iNode-1)+5-1) = zr(idch2-1+ival+3)
    end do
!
    lisno(1+5*nbNode+1-1) = beamNodeName
    lisno(1+5*nbNode+2-1) = beamNodeName
    lisno(1+5*nbNode+3-1) = beamNodeName
!
    lisddl(1+5*nbNode+1-1) = 'DRX'
    lisddl(1+5*nbNode+2-1) = 'DRY'
    lisddl(1+5*nbNode+3-1) = 'DRZ'
!
    coer(1+5*nbNode+1-1) = -ig(1)
    coer(1+5*nbNode+2-1) = -ig(2)
    coer(1+5*nbNode+3-1) = -ig(3)
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, valeType, 0.d0, listRela)
    call imprel(factorKeyword, nbterm, coer, lisddl, lisno, &
                beta)

    nbterm = 5*nbNode+3
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
        ival = zi(iaprno+(nodeNume-1)*(nbec+2)+1-1)-1
!
        lisno(1+5*(iNode-1)+1-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+2-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+3-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+4-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+5-1) = nodeName(iNode)
!
        lisddl(1+5*(iNode-1)+1-1) = 'DX'
        lisddl(1+5*(iNode-1)+2-1) = 'DZ'
        lisddl(1+5*(iNode-1)+3-1) = 'DRX'
        lisddl(1+5*(iNode-1)+4-1) = 'DRY'
        lisddl(1+5*(iNode-1)+5-1) = 'DRZ'
!
        coer(1+5*(iNode-1)+1-1) = zr(idch1-1+ival+3)
        coer(1+5*(iNode-1)+2-1) = -zr(idch1-1+ival+1)
        coer(1+5*(iNode-1)+3-1) = zr(idch2-1+ival+2)
        coer(1+5*(iNode-1)+4-1) = zr(idch2-1+ival+4)
        coer(1+5*(iNode-1)+5-1) = zr(idch2-1+ival+5)
    end do
!
    lisno(1+5*nbNode+1-1) = beamNodeName
    lisno(1+5*nbNode+2-1) = beamNodeName
    lisno(1+5*nbNode+3-1) = beamNodeName
!
    lisddl(1+5*nbNode+1-1) = 'DRX'
    lisddl(1+5*nbNode+2-1) = 'DRY'
    lisddl(1+5*nbNode+3-1) = 'DRZ'
!
    coer(1+5*nbNode+1-1) = -ig(2)
    coer(1+5*nbNode+2-1) = -ig(4)
    coer(1+5*nbNode+3-1) = -ig(5)
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, valeType, 0.d0, listRela)
    call imprel(factorKeyword, nbterm, coer, lisddl, lisno, &
                beta)

    nbterm = 5*nbNode+3
    do iNode = 1, nbNode
        nodeNume = char8_to_int(nodeName(iNode))
! ---    ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES CHAMNO
        ival = zi(iaprno+(nodeNume-1)*(nbec+2)+1-1)-1
!
        lisno(1+5*(iNode-1)+1-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+2-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+3-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+4-1) = nodeName(iNode)
        lisno(1+5*(iNode-1)+5-1) = nodeName(iNode)
!
        lisddl(1+5*(iNode-1)+1-1) = 'DY'
        lisddl(1+5*(iNode-1)+2-1) = 'DX'
        lisddl(1+5*(iNode-1)+3-1) = 'DRX'
        lisddl(1+5*(iNode-1)+4-1) = 'DRY'
        lisddl(1+5*(iNode-1)+5-1) = 'DRZ'
!
        coer(1+5*(iNode-1)+1-1) = zr(idch1-1+ival+1)
        coer(1+5*(iNode-1)+2-1) = -zr(idch1-1+ival+2)
        coer(1+5*(iNode-1)+3-1) = zr(idch2-1+ival+3)
        coer(1+5*(iNode-1)+4-1) = zr(idch2-1+ival+5)
        coer(1+5*(iNode-1)+5-1) = zr(idch2-1+ival+6)
    end do
!
    lisno(1+5*nbNode+1-1) = beamNodeName
    lisno(1+5*nbNode+2-1) = beamNodeName
    lisno(1+5*nbNode+3-1) = beamNodeName
!
    lisddl(1+5*nbNode+1-1) = 'DRX'
    lisddl(1+5*nbNode+2-1) = 'DRY'
    lisddl(1+5*nbNode+3-1) = 'DRZ'
!
    coer(1+5*nbNode+1-1) = -ig(3)
    coer(1+5*nbNode+2-1) = -ig(5)
    coer(1+5*nbNode+3-1) = -ig(6)
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, valeType, 0.d0, listRela)
    call imprel(factorKeyword, nbterm, coer, lisddl, lisno, &
                beta)
!
    if ((option .eq. 'COQ_TUYAU')) then
        call racotu(zi(iaprno), nbNode, nodeName, beamNodeName, mesh, &
                    ligrel, model, caraElem, numeDof, listRela, coorig)
    end if

! - Clean
    call detrsd('LIGREL', ligrel)
    call jedetr('&&RAPOCO.LISTE_NOEUDS')
    call jedetr('&&RAPOCO.LISTE_MAILLES')
    call detrsd('CARTE', mapBeamAxis)
    call detrsd('CHAMP_GD', fieldSection)
    AS_DEALLOCATE(vr=inertie_raccord)
    call detrsd('CARTE', mapSection)
    call detrsd('RESUELEM', '&&RAPOCO.VECT2')
    call detrsd('RESUELEM', '&&RAPOCO.VECT_XYZNI')
    call jedetr('&&RAPOCO           .RELR')
    call jedetr('&&RAPOCO           .RERR')
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
