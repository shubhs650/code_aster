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
subroutine eclpgm(ma2, model1, cham1, ligrelIn, shrink, &
                  edgeMin, nbField, listFieldType)
!
    implicit none
!
#include "asterfort/alchml.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celfpg.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/eclapp.h"
#include "asterfort/eclaty.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: ma2, model1
    character(len=19), intent(in) :: cham1, ligrelIn
    real(kind=8), intent(in) :: shrink, edgeMin
    integer(kind=8), intent(in) :: nbField
    character(len=16), intent(in) :: listFieldType(nbField)
!
! --------------------------------------------------------------------------------------------------
!
! TRAITEMENT DU MOT CLE CREA_MAILLAGE/ECLA_PG
!
! BUT : CREER LE MAILLAGE (MA2) CORRESPONDANT AUX POINTS DE GAUSS
!       DES ELEMENTS DU LIGREL (LIGREL).
!
! --------------------------------------------------------------------------------------------------
!
!  IN/JXIN  MO : MODELE DONT ON VEUT "ECLATER" LES POINTS DE GAUSS
!  IN/JXIN  LIGREL : NOM DU LIGREL (D'ELEMENTS DU MODELE) CORRESPONDANT
!           AUX MAILLES QUI INTERESSENT L'UTILISATEUR.
!           SI LIGREL=' ', ON PREND LE LIGREL DU MODELE (TOUT='OUI')
!  IN/JXIN  CHAM1 : NOM DU CHAM_ELEM/ELGA QUE L'ON VEUT ECLATER
!           (OU ' ' SI ON UTILISE NCH ET LISCH)
!  IN       LISCH(*) : LISTE DES NOMS SYMBOLIQUES DES CHAMPS AUXQUELS
!           ON S'INTERESSE. EX: (SIEF_ELGA, VARI_ELGA, ...)
!  IN NCH : DIMENSION DE LISCH  (0 SI CHAM1/=' ')
!  IN/JXOUT MA2 : NOM DU MAILLAGE A CREER (1 MAILLE PAR POINT DE GAUSS)
!           LES MAILLES DE MA2 SONT TOUTES DECONNECTEES LES UNES DES
!           DES AUTRES.
!  IN  SHRINK : COEFFICENT DE "CONTRACTION" DES MAILLES POUR QUE
!       CELLES-CI NE SE CHEVAUCHENT PAS.
!       SHRINK=0.9 DONNE UN BEAU VITRAIL AVEC 10% DE NOIR.
!  IN  LONMIN : LONGUEUR MINIMALE POUR LES ARETES DES MAILLES CREEES.
!      (IGNORE SI LONMIN <=0)
!      CETTE FONCTIONNALITE PERMET D'"EPAISSIR" DES ELEMENTS DE JOINT
!      POUR QU'ON PUISSE LES VOIR.
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: mxnbn2 = 8, mxnbpi = 64, mxnbte = 8, mxnbse = 27
!     MXNBN2 : MAX DU NOMBRE DE NOEUDS D'UN SOUS-ELEMENT (HEXA8)
!     MXNBPI : MAX DU NOMBRE DE POINT_I (HEXA A 27 POINTS DE GAUSS)
!     MXNBPI = 4X4X4
!     MXNBTE : MAX DU NOMBRE DE TERMES DE LA C.L. DEFINISSANT 1 POINT_I
!              AU PLUS LES 8 SOMMETS D'UN HEXA8
!     MXNBSE : MAX DU NOMBRE DE SOUS-ELEMENTS
    integer(kind=8) :: k, typeElemNume, tabno(MT_NNOMAX), iret1, jobj, numa
    integer(kind=8) :: igr, iel, ilmaco, illiel
    integer(kind=8) :: dimgeo, ibid, ino, ino1, ino2
    integer(kind=8) :: nbmail, nbnoeu, nbcoor, iadime, kse
    integer(kind=8) :: nbno2t, ianno2, iatypm, nuno2, nupoi2, cas
    integer(kind=8) :: npg, nbpi, iagese, nno2, nuse, nbSel
    integer(kind=8) :: ima, nbelgr, nupoin, npoini, iterm, ipoini
    integer(kind=8) :: iret, iField, nch2
    character(len=8) :: ma1, elrefa, fapg, nompar
    character(len=8) :: fieldDisc
    character(len=16) :: typeElemName
    character(len=19) :: modelLigrel, cel, fieldLigrel, ligrel
    character(len=24) :: nomobj
    character(len=24) :: valk(2)
    real(kind=8) :: x, xc, xm
    integer(kind=8) :: nbse, corsel(mxnbse)
    integer(kind=8) :: connx(mxnbn2, mxnbse), nsomm1(mxnbpi, mxnbte)
    integer(kind=8) :: nterm1(mxnbpi), nbno2(mxnbse), typeCellNume(mxnbse)
    real(kind=8) :: csomm1(mxnbpi, mxnbte)
    integer(kind=8) :: ico
    integer(kind=8) :: opt, iadesc, iaoppa, nbin
    integer(kind=8), pointer :: connexse(:) => null()
    real(kind=8), pointer :: geopoini(:) => null()
    integer(kind=8), pointer :: liel(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: maco(:) => null()
! --------------------------------------------------------------------------------------------------
!     NUMAIL(IGR,IEL)=NUMERO DE LA MAILLE ASSOCIEE A L'ELEMENT IEL
#define numail(igr,iel) liel(zi(illiel+igr-1)+iel-1)
!     NBNOMA(IMA)=NOMBRE DE NOEUDS DE LA MAILLE IMA
#define nbnoma(ima) zi(ilmaco-1+ima+1)-zi(ilmaco-1+ima)
!     NUMGLM(IMA,INO)=NUMERO GLOBAL DU NOEUD INO DE LA MAILLE IMA
!                     IMA ETANT UNE MAILLE DU MAILLAGE.
#define numglm(ima,ino) maco(zi(ilmaco+ima-1)+ino-1)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    dimgeo = 3
!
    call dismoi('NOM_MAILLA', model1, 'MODELE', repk=ma1)
    call jeveuo(ma1//'.COORDO    .VALE', 'L', vr=vale)
    call jeveuo(ma1//'.CONNEX', 'L', vi=maco)
    call jeveuo(jexatr(ma1//'.CONNEX', 'LONCUM'), 'L', ilmaco)
!
!
!     -- ON ACCEPTE 2 CAS DE FIGURE :
    if (cham1 .ne. ' ') then
!       -- CAS "PROJ_CHAMP" :
        call dismoi('TYPE_CHAMP', cham1, 'CHAMP', repk=fieldDisc)
        ASSERT(fieldDisc .eq. 'ELGA')
        ASSERT(nbField .eq. 0)
        ASSERT(ligrelIn .ne. ' ')
        call dismoi('NOM_LIGREL', cham1, 'CHAM_ELEM', repk=fieldLigrel)
        ASSERT(fieldLigrel .eq. ligrelIn)
        cas = 1
    else
!       -- CAS "CREA_MAILLAGE/ECLA_PG" :
        cas = 2
    end if
!
!
!
!    0.1 : ON CHERCHE LE NOM DES FAMILLES DE POINTS DE GAUSS A ECLATER
!          SORTIES : NCH2 [+ JOBJ SI NCH2>0]
!    -----------------------------------------------------------------
    nomobj = '&&ECLPGM.NOMOBJ'
!
    if (cas .eq. 1) then
        nch2 = 1
        call celfpg(cham1, nomobj, iret)
        ASSERT(iret .eq. 0)
        call jeveuo(nomobj, 'L', jobj)
!
    else if (cas .eq. 2) then
        cel = '&&ECLPGM.CHAM_ELEM'
        ico = 0
        nch2 = nbField
        call dismoi('NOM_LIGREL', model1, 'MODELE', repk=modelLigrel)
        do iField = 1, nch2
            if (listFieldType(iField) (6:9) .ne. 'ELGA') then
                call utmess('F', 'CALCULEL2_41')
            end if
            call jenonu(jexnom('&CATA.OP.NOMOPT', listFieldType(iField)), opt)
            if (opt .eq. 0) goto 30
!
            call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iadesc)
            call jeveuo(jexnum('&CATA.OP.OPTPARA', opt), 'L', iaoppa)
            nbin = zi(iadesc-1+2)
            nompar = zk8(iaoppa-1+nbin+1)
!
            call alchml(modelLigrel, listFieldType(iField), nompar, 'V', cel, &
                        iret1, ' ')
            if (iret1 .ne. 0) goto 30
!
            ico = ico+1
            call celfpg(cel, nomobj, iret)
            call detrsd('CHAMP', cel)
            if (iret .eq. 1) then
                valk(1) = model1
                valk(2) = listFieldType(iField)
                call utmess('F', 'CALCULEL2_33', nk=2, valk=valk)
            end if
            call jeveuo(nomobj, 'L', jobj)
30          continue
        end do
!        -- ON N'A PAS TROUVE DE CHAMP ELGA CORRECT :
        if (ico .eq. 0) then
            valk(1) = listFieldType(1)
            call utmess('F', 'CALCULEL2_24', nk=1, valk=valk)
        end if
    end if
!
!
!
!     0.2 : ON SE RESTREINT AUX MAILLES EVENTUELLEMENT SPECIFIEES PAR
!           L'UTILISATEUR :
!     ----------------------------------------------------------------
    if (ligrelIn .eq. ' ') then
        ligrel = modelLigrel
    else
        ligrel = ligrelIn
    end if
    call jeveuo(ligrel//'.LIEL', 'L', vi=liel)
    call jeveuo(jexatr(ligrel//'.LIEL', 'LONCUM'), 'L', illiel)
!
!
!
!     1. ON COMPTE LES FUTURS SOUS-ELEMENTS
!        ET LES POINT_I ET LES NOEUDS DU FUTUR MAILLAGE
!     ---------------------------------------------------------------
    nbse = 0
    nbpi = 0
    nbno2t = 0
    do igr = 1, nbgrel(ligrel)
        typeElemNume = typele(ligrel, igr)
        nbelgr = nbelem(ligrel, igr)
        if (nbelgr .eq. 0) cycle
        call jenuno(jexnum('&CATA.TE.NOMTE', typeElemNume), typeElemName)
        ASSERT(nch2 .gt. 0)
        numa = numail(igr, 1)
        elrefa = zk16(jobj-1+numa) (1:8)
        fapg = zk16(jobj-1+numa) (9:16)
        if (fapg .eq. ' ') cycle
!
        call eclaty(typeElemName, &
                    elrefa, fapg, &
                    mxnbn2, mxnbpi, mxnbte, mxnbse, &
                    npg, npoini, &
                    nterm1, nsomm1, csomm1, &
                    typeCellNume, nbno2, connx, &
                    nbsel, corsel, iret)
        if (iret .ne. 0) then
            call utmess("F", "PROJECTION4_7")
        end if

        nbse = nbse+nbelgr*nbSel
        nbpi = nbpi+nbelgr*npoini
        do kse = 1, nbSel
            nbno2t = nbno2t+nbelgr*nbno2(kse)
        end do

    end do
!
!
!     2. ON ALLOUE 4 OBJETS DE TRAVAIL (+ MAILLAGE//'.TYPMAIL') :
!        .GEOPOINI : GEOMETRIE DES POINT_I
!        .CONNEXSE : NUMEROS DES POINT_I DES SOUS-ELEMENTS
!        .GEOSE    : GEOMETRIE DES SOUS-ELEMENTS
!        .NBNO2    : NOMBRE DE NOEUDS DES SOUS-ELEMENTS
!        .TYPMAIL  : TYPE_MAILLE DES SOUS-ELEMENTS
!     ---------------------------------------------------------------
    if (nbpi .eq. 0) then
        call utmess('F', 'CALCULEL2_35')
    end if
    AS_ALLOCATE(vr=geopoini, size=nbpi*dimgeo)
    AS_ALLOCATE(vi=connexse, size=nbno2t)
    call wkvect('&&ECLPGM.GEOSE', 'V V R', nbno2t*dimgeo, iagese)
    call wkvect('&&ECLPGM.NBNO2', 'V V I', nbse, ianno2)
!
    if (nbse .eq. 0) then
        call utmess('F', 'CALCULEL2_36')
    end if
    call wkvect(ma2//'.TYPMAIL', 'G V I', nbse, iatypm)
!
!
!
!     3. ON CALCULE DES COORDONNEES DES SOUS-ELEMENTS
!        ET LEUR CONNECTIVITE
!     ---------------------------------------------------------------
!
    nuse = 0
    nupoin = 0
    nuno2 = 0
    do igr = 1, nbgrel(ligrel)
        typeElemNume = typele(ligrel, igr)
        nbelgr = nbelem(ligrel, igr)
        if (nbelgr .eq. 0) cycle
        call jenuno(jexnum('&CATA.TE.NOMTE', typeElemNume), typeElemName)
!
        ASSERT(nch2 .gt. 0)
        numa = numail(igr, 1)
        elrefa = zk16(jobj-1+numa) (1:8)
        fapg = zk16(jobj-1+numa) (9:16)
        if (fapg .eq. ' ') cycle
!
        call eclaty(typeElemName, &
                    elrefa, fapg, &
                    mxnbn2, mxnbpi, mxnbte, mxnbse, &
                    npg, npoini, &
                    nterm1, nsomm1, csomm1, &
                    typeCellNume, nbno2, connx, &
                    nbsel, corsel, iret)
        if (iret .ne. 0) then
            call utmess("F", "PROJECTION4_7")
        end if
        if (nbSel .eq. 0) cycle
!
        do iel = 1, nbelgr
!          ON RECUPERE LE NUMERO DE LA MAILLE ET LE NUMERO
!          DE SES SOMMETS :
            ima = numail(igr, iel)
            do ino1 = 1, nbnoma(ima)
                tabno(ino1) = numglm(ima, ino1)
            end do
!
!          -- CALCUL DES COORDONNEES DES POINT_I :
            do ipoini = 1, npoini
                nupoin = nupoin+1
                do k = 1, dimgeo
                    x = 0.d0
                    do iterm = 1, nterm1(ipoini)
                        x = x+csomm1(ipoini, iterm)*vale(3*( &
                                                         tabno(nsomm1(ipoini, iterm))-1)+k)
                    end do
                    geopoini((nupoin-1)*dimgeo+k) = x
                end do
            end do
!
!          -- STOCKAGE DES NUMEROS DES POINT_I DES SOUS-ELEMENTS
!             ET DE LEURS COORDONNEES :
            do kse = 1, nbSel
                nuse = nuse+1
                zi(iatypm-1+nuse) = typeCellNume(kse)
                zi(ianno2-1+nuse) = nbno2(kse)
                nno2 = nbno2(kse)
                do ino2 = 1, nno2
                    nuno2 = nuno2+1
                    connexse(nuno2) = connx(ino2, kse)+nupoin-npoini
                    nupoi2 = connexse(nuno2)
                    do k = 1, dimgeo
                        zr(iagese-1+(nuno2-1)*dimgeo+k) = geopoini(( &
                                                                   nupoi2-1)*dimgeo+k)
                    end do
                end do
!           DANS LE CAS DU QUADRILATERE ON CONTROLE L'APPLATISSEMENT
                if (nno2 .eq. 4) then
                    call eclapp(dimgeo, nno2, edgeMin, zr(iagese+(nuno2-4)*dimgeo))
                end if
            end do
!
        end do
    end do
!
    call jedetr(nomobj)
!
!     3. CONSTRUCTION DES OBJETS DU MAILLAGE RESULTAT :
!     -------------------------------------------------
    nbmail = nbse
    nbnoeu = nbno2t
    nbcoor = dimgeo
!
!     3.1 CREATION DE L'OBJET .DIME  :
!     ------------------------------------
    call wkvect(ma2//'.DIME', 'G V I', 6, iadime)
    zi(iadime-1+1) = nbnoeu
    zi(iadime-1+3) = nbmail
    zi(iadime-1+6) = nbcoor
!
!
!     3.3 CREATION DES OBJETS  .CONNEX ET .TYPMAIL
!     ---------------------------------------------
    call jecrec(ma2//'.CONNEX', 'G V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbmail)
    call jeecra(ma2//'.CONNEX', 'LONT', nbnoeu, ' ')
    call jeveuo(ma2//'.CONNEX', 'E', ibid)
!
    nuno2 = 0
    do ima = 1, nbmail
        nno2 = zi(ianno2-1+ima)
        call jecroc(jexnum(ma2//'.CONNEX', ima))
        call jeecra(jexnum(ma2//'.CONNEX', ima), 'LONMAX', nno2)
        do ino2 = 1, nno2
            nuno2 = nuno2+1
            zi(ibid-1+nuno2) = nuno2
        end do
    end do
!
!
!
!     3.4 CREATION DU CHAMP DE GEOMETRIE (.COORDO)
!     ---------------------------------------------
    call copisd('CHAMP_GD', 'G', ma1//'.COORDO', ma2//'.COORDO')
    call jedetr(ma2//'.COORDO    .VALE')
    call wkvect(ma2//'.COORDO    .VALE', 'G V R', 3*nbnoeu, ibid)
!
    do k = 1, dimgeo
        nuno2 = 0
        do ima = 1, nbmail
            nno2 = zi(ianno2-1+ima)
!         -- ON FAIT UN PETIT "SHRINK" SUR LES MAILLES :
            xc = 0.d0
            do ino = 1, nno2
                nuno2 = nuno2+1
                xc = xc+zr(iagese-1+(nuno2-1)*dimgeo+k)/dble(nno2)
            end do
            nuno2 = nuno2-nno2
!
            do ino = 1, nno2
                nuno2 = nuno2+1
                xm = zr(iagese-1+(nuno2-1)*dimgeo+k)
                xm = xc+shrink*(xm-xc)
                zr(ibid-1+(nuno2-1)*3+k) = xm
            end do
        end do
    end do
!
!
!
    call jedetr('&&ECLPGM.NOMOBJ')
    AS_DEALLOCATE(vr=geopoini)
    AS_DEALLOCATE(vi=connexse)
    call jedetr('&&ECLPGM.GEOSE')
    call jedetr('&&ECLPGM.NBBO2')
    call jedema()
!
! FIN ------------------------------------------------------------------
end subroutine
