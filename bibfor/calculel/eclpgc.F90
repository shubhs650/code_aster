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
subroutine eclpgc(ch1, ch2, ligrel, mesh2, numeEqua, &
                  nomfpg)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celcel.h"
#include "asterfort/celver.h"
#include "asterfort/chligr.h"
#include "asterfort/cmpcha.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/eclaty.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: ch1, ch2, ligrel
    character(len=8), intent(in) :: mesh2
    character(len=19), intent(in) :: numeEqua
    character(len=24), intent(in) :: nomfpg
!
! --------------------------------------------------------------------------------------------------
!
! BUT : "ECLATER" LE CHAM_ELEM_ELGA CH1 POUR CREER LE CHAM_NO CH2 SUR
!        SUR LE MAILLAGE MA2.
!
! ARGUMENTS :
!  IN/JXN   CH1 : CHAM_ELEM_ELGA A ECLATER
!  IN/JXOUT CH2 : NOM DU CHAM_NO A CREER
!  IN/JXIN  LIGREL : NOM DU LIGREL CORRESPONDANT AUX MAILLES QUI
!           INTERESSENT L'UTILISATEUR.
!           LIGREL EST EVENTUELLEMMENT UN "SOUS" LIGREL DU LIGREL
!           ASSOCIE A CH1.
!  IN/JXIN  MA2 : NOM DU MAILLAGE QUI "PORTERA" LE CHAM_NO CH2
!  IN/JXIN  NUMEQ : NOM DE LA SD_NUME_EQUA QUI SERA ASSOCIEE A CH2
!  IN/JXIN  NOMFPG : NOM D'UN OBJET JEVEUX CONTENANT LE NOM DE LA
!           FAMILLE DE PG A UTILISER POUR CHAQUE MAILLE DE LIGREL.
!
! REMARQUES :
!  * MA2 IL EST TRES IMPORTANT QUE LE MAILLAGE MA2 FOURNI SOIT CELUI
!      QUI A ETE OBTENU PAR LA ROUTINE ECLPGM (AVEC LE MEME
!      LIGREL EN ENTREE).
!      LA "JUSTESSE" DU CHAMP CREE (CH2) PROVIENT DU FAIT QUE DANS
!      LES 2 ROUTINES ECLPGM ET ECLPGC, ON PARCOURT LES MEMES MAILLES
!      DANS LE MEME ORDRE.
!  * NOMFPG PEUT ETRE OBTENU PAR LA ROUTINE CELFPG
!  * LIGREL PEUT ETRE OBTENU PAR LA ROUTINE EXLIMA OU BIEN ON
!           L'EXTRAIT DE CH1 (DISMOI).
!  * NUMEQ PEUT ETRE ' '. DANS CE CAS, ON EN CREERA UN DIFFERENT
!           A CHAQUE APPEL A ECLPGC.
!           CHOISIR NUMEQ /= ' ' PERMET D'ECONOMISER CETTE SD SI
!           PLUSIEURS CHAMPS PEUVENT LA PARTAGER.
!           C'EST LE CAS EN GENERAL POUR LA BOUCLE SUR LES NUME_ORDRE
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
    integer(kind=8) :: corsel(mxnbse), nbSel
    integer(kind=8) :: connx(mxnbn2, mxnbse), nsomm1(mxnbpi, mxnbte)
    integer(kind=8) :: nterm1(mxnbpi), nbno2(mxnbse), typeCellNume(mxnbse)
    real(kind=8) :: csomm1(mxnbpi, mxnbte)
    integer(kind=8), parameter ::  mxcmp = 200
    aster_logical :: lvari
    integer(kind=8) :: numa, jnofpg, kk, iret
    integer(kind=8) :: k, typeElemNume, npg, npoini, ideca2
    integer(kind=8) :: igr, jcmaco, jcliel, jcnsl2
    integer(kind=8) :: ibid, nbpg, ino, nbgr, inogl, kse
    integer(kind=8) :: iamol1, jcnsv2
    integer(kind=8) :: ima, nbelgr, jval2, nbno, nddl, iddl, adiel
    integer(kind=8) :: iipg, jceld1, moloc1, nb_cmp_mx, jcmpgd
    integer(kind=8) :: nuddl(mxcmp), mxvari, iel, nb_cmp
    character(len=8) :: nomg1, nomg2, elrefa, fapg
    character(len=16) :: typeElemName
    character(len=16) :: optio, param
    character(len=19) :: ch2s, ch1b
    character(len=24) :: valk(2)
    character(len=24), pointer :: celk(:) => null()
    real(kind=8), pointer :: celv(:) => null()
    integer(kind=8), pointer :: liel(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: cata_to_field(:) => null()
    integer(kind=8), pointer :: field_to_cata(:) => null()
    character(len=8), pointer :: cmp_name(:) => null()
! --------------------------------------------------------------------------------------------------
!     NBNOMA(IMA)=NOMBRE DE NOEUDS DE LA MAILLE IMA
#define nbnoma(ima) zi(jcmaco-1+ima+1)-zi(jcmaco-1+ima)
!     NUMGLM(IMA,INO)=NUMERO GLOBAL DU NOEUD INO DE LA MAILLE IMA
!                     IMA ETANT UNE MAILLE DU MAILLAGE.
#define numglm(ima,ino) connex(zi(jcmaco+ima-1)+ino-1)
#define numail(igr,iel) liel(zi(jcliel+igr-1)+iel-1)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
!
    ch1b = ch1
    call dismoi('NOM_GD', ch1b, 'CHAMP', repk=nomg1)
    if (nomg1(5:6) .ne. '_R') then
        call utmess('F', 'CALCULEL2_39')
    end if
    nomg2 = nomg1
    lvari = (nomg1 .eq. 'VARI_R')
    if (lvari) nomg2 = 'VAR2_R'
!
!
!
!     -- ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
    call celver(ch1b, 'NBSPT_1', 'COOL', kk)
    if (kk .eq. 1) then
        call utmess('I', 'PREPOST_36', sk=nomg1)
        call celcel('PAS_DE_SP', ch1b, 'V', '&&ECLPGC.CH1B1')
        ch1b = '&&ECLPGC.CH1B1'
    end if
!
!     -- PROJECTION SUR LE LIGREL REDUIT SI NECESSAIRE :
    call jeveuo(ch1b//'.CELK', 'L', vk24=celk)
    if (celk(3) (1:4) .ne. 'ELGA') then
        call utmess('F', 'CALCULEL2_41')
    end if
    if (celk(1) (1:19) .ne. ligrel) then
        optio = celk(2) (1:16)
        param = celk(6) (1:16)
        call chligr(ch1b, ligrel, optio, param, 'V', &
                    '&&ECLPGC.CH1B2')
        ch1b = '&&ECLPGC.CH1B2'
    end if
!
    call jeexin(ch1b//'.CELD', ibid)
    if (ibid .eq. 0) goto 90
!
    call jeveuo(ch1b//'.CELV', 'L', vr=celv)
    call jeveuo(ch1b//'.CELD', 'L', jceld1)
    call jeveuo(ch1b//'.CELD', 'L', jceld1)
!     -- MXVARI : NOMBRE MAXI DE VXX SI VARI_R
    mxvari = max(1, zi(jceld1-1+4))
!
!
!     -- ON CHERCHE LES CMPS PRESENTES DANS LE CHAM_ELEM CH1B :
!         NCMP : NOMBRE DE CMPS PRESENTES
!         '&&ECLPGC.CORR1': CONTIENT LA CORRESPONDANCE ENTRE LE
!                           NUMERO D'1 CMP DU CHAM_ELEM ET LE
!                           NUMERO D'1 CMP DU CHAM_ELEM_S
!         '&&ECLPGC.NOM_CMP': CONTIENT LES NOMS DES CMPS DU CHAM_ELEM_S
!     -----------------------------------------------------------------
    if (.not. lvari) then
!
! ----- Create objects for global components (catalog) <=> local components (field)
!
        call cmpcha(ch1b, cmp_name, cata_to_field, field_to_cata, nb_cmp, &
                    nb_cmp_mx)
    else
!       -- POUR VARI_R :
        ASSERT(nomg2 .eq. 'VAR2_R')
        call dismoi('NB_CMP_MAX', nomg2, 'GRANDEUR', repi=nb_cmp_mx)
        ASSERT(mxvari .le. nb_cmp_mx)
        nb_cmp = mxvari
        call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg2), 'L', jcmpgd)
        AS_ALLOCATE(vk8=cmp_name, size=nb_cmp)
        AS_ALLOCATE(vi=cata_to_field, size=nb_cmp)
        do k = 1, nb_cmp
            cata_to_field(k) = k
            cmp_name(k) = zk8(jcmpgd-1+k)
        end do
    end if
    if (nb_cmp .gt. mxcmp) call utmess('F', 'CALCULEL2_51', ni=2, vali=[mxcmp, nb_cmp])

! - CREATION D'UN CHAM_NO_S : CH2S
    ch2s = '&&ECLPGC.CH2S'
    call cnscre(mesh2, nomg2, nb_cmp, cmp_name, 'V', &
                ch2s)
    call jeveuo(ch2s//'.CNSV', 'E', jcnsv2)
    call jeveuo(ch2s//'.CNSL', 'E', jcnsl2)

! - REMPLISSAGE DU CHAM_NO
    call jeveuo(nomfpg, 'L', jnofpg)
    call jeveuo(ligrel//'.LIEL', 'L', vi=liel)
    call jeveuo(jexatr(ligrel//'.LIEL', 'LONCUM'), 'L', jcliel)
    call jeveuo(mesh2//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(mesh2//'.CONNEX', 'LONCUM'), 'L', jcmaco)

    ima = 0
    nbgr = nbgrel(ligrel)
    do igr = 1, nbgr
        nbelgr = nbelem(ligrel, igr)
        if (nbelgr .eq. 0) cycle
        moloc1 = zi(jceld1-1+zi(jceld1-1+4+igr)+2)
        if (moloc1 .eq. 0) goto 80
!
        if (.not. lvari) then
            ASSERT(mxvari .eq. 1)
        end if
!
        call jeveuo(jexnum('&CATA.TE.MODELOC', moloc1), 'L', iamol1)
        ASSERT(zi(iamol1-1+1) .eq. 3)
        nbpg = zi(iamol1-1+4)
!
        numa = numail(igr, 1)
        elrefa = zk16(jnofpg-1+numa) (1:8)
        fapg = zk16(jnofpg-1+numa) (9:16)
        if (fapg .eq. ' ') goto 80
!
!           -- ON VERIFIE QUE C'EST UN CHAMP "ELGA/IDEN" :
!           ----------------------------------------------
        ASSERT(.not. ((nbpg .lt. 0) .or. (nbpg .gt. 10000)))
!
!           -- ON ECLATE LE TYPE_ELEM :
!           ---------------------------
        typeElemNume = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', typeElemNume), typeElemName)
        call eclaty(typeElemName, &
                    elrefa, fapg, &
                    mxnbn2, mxnbpi, mxnbte, mxnbse, &
                    npg, npoini, &
                    nterm1, nsomm1, csomm1, &
                    typeCellNume, nbno2, connx, &
                    nbsel, corsel, iret)
        if (iret .ne. 0) then
            call utmess("A", "PROJECTION4_6", sk=typeElemName, si=nbelgr)
        end if
        if (npg .ne. 0) then
            if (nbpg .ne. npg) then
                valk(1) = typeElemName
                call utmess('F', 'CALCULEL2_42', sk=valk(1))
            end if
        else
!            -- ON IGNORE LES AUTRES ELEMENTS :
            ASSERT(nbSel .eq. 0)
            goto 80
        end if

!            -- QUELLES SONT LES CMPS PORTEES PAR LES POINTS DE GAUSS ?
!            ----------------------------------------------------------
        if (lvari) then
            nddl = mxvari
            do k = 1, nddl
                nuddl(k) = k
            end do
        else
            nddl = 0
            do k = 1, nb_cmp_mx
                if (exisdg(zi(iamol1-1+4+1), k)) then
                    nddl = nddl+1
                    nuddl(nddl) = cata_to_field(k)
                end if
            end do
        end if
!
!
!          -- BOUCLE SUR TOUS LES POINTS DE GAUSS DU GREL :
!          ------------------------------------------------
        do iel = 1, nbelgr
            if (lvari) nddl = zi(jceld1-1+zi(jceld1-1+4+igr)+4+4*(iel-1)+2)
!
            do kse = 1, nbSel
!            -- AU POINT DE GAUSS IIPG CORRESPOND LA MAILLE NUMERO IMA
!               DANS MA2.
                iipg = corsel(kse)
                ima = ima+1
!
                nbno = nbnoma(ima)
                if (nbno .gt. MT_NNOMAX) then
                    call utmess('F', 'CALCULEL2_43')
                end if
                do ino = 1, nbno
                    inogl = numglm(ima, ino)
                    do iddl = 1, nddl
                        ideca2 = nb_cmp*(inogl-1)+nuddl(iddl)
                        jval2 = jcnsv2-1+ideca2
                        zl(jcnsl2-1+ideca2) = .true.
                        adiel = zi(jceld1-1+zi(jceld1-1+4+igr)+4+4*(iel- &
                                                                    1)+4)
                        zr(jval2) = celv(adiel-1+nddl*(iipg-1)+ &
                                         iddl)
                    end do
                end do
            end do
        end do
80      continue
    end do
!
!     -- ON ESSAYE DE FAIRE UN PEU DE PLACE EN MEMOIRE :
    call detrsd('CHAMP_GD', '&&ECLPGC.CH1B1')
    call detrsd('CHAMP_GD', '&&ECLPGC.CH1B2')
    call jelibe(ch1//'.CELD')
    call jelibe(ch1//'.CELV')
    call jelibe(mesh2//'.CONNEX')
!
!
!         -- ON TRANSFORME CH2S EN VRAI CHAM_NO :
!         ----------------------------------------
!     -- 2 JELIBE POUR ECONOMISER LA MEMOIRE EN DESSOUS (MARQUE):
    call jelibe(ch2s//'.CNSL')
    call jelibe(ch2s//'.CNSV')
    call cnscno(ch2s, numeEqua, 'NON', 'G', ch2, &
                'F', ibid)
    call detrsd('CHAM_NO_S', ch2s)
!
    AS_DEALLOCATE(vi=cata_to_field)
    AS_DEALLOCATE(vi=field_to_cata)
    AS_DEALLOCATE(vk8=cmp_name)
!
!
90  continue
    call detrsd('CHAMP_GD', '&&ECLPGC.CH1B1')
    call detrsd('CHAMP_GD', '&&ECLPGC.CH1B2')
    call jedema()
end subroutine
