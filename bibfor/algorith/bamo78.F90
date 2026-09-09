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
subroutine bamo78(nomres, trange, typres)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesfus.h"
#include "asterfort/compStress.h"
#include "asterfort/copmod.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/dyna_comp_fuse.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mdgeph.h"
#include "asterfort/mechti.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/refdcp.h"
#include "asterfort/rs_getlast.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/rsorac.h"
#include "asterfort/rstran.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/vtcreb.h"
!
    character(len=8) :: nomres
    character(len=16) :: typres
    character(len=19) :: trange
! IN  : NOMRES : NOM UTILISATEUR POUR LA COMMANDE REST_COND_TRAN
! IN  : TYPRES : TYPE DE RESULTAT : 'DYNA_TRANS'
! IN  : TRANGE : NOM UTILISATEUR DU CONCEPT TRAN_GENE AMONT
!
!
!
!
    character(len=8) :: k8bid
    integer(kind=8) :: ibid, iret, iretou
    integer(kind=8) :: icham, iarch
    complex(kind=8) :: c16bid
    integer(kind=8) :: nbcham, nume
    character(len=16) :: champ(3)
    integer(kind=8) :: n0, n1
    character(len=8) :: basemo
    integer(kind=8) :: neq
    integer(kind=8) :: nbinst
    integer(kind=8) :: nbmode
    integer(kind=8) :: jrestr, ldnew, linst
    character(len=14) :: numddl
    character(len=24) :: numedd
    character(len=19), parameter :: ches1 = '&&BAMO78.CHES1'
    character(len=19), parameter :: ches2 = '&&BAMO78.CHES2'
    character(len=19), parameter :: ches3 = '&&BAMO78.CHES3'
    character(len=19), parameter :: chel2 = '&&BAMO78.CHEL2'
    character(len=16), parameter :: opti(2) = (/'SIEF_ELGA', 'VARI_ELGA'/)
    character(len=24), parameter :: chvarc = '&&BAMO78.VARC', chvref = '&&BAMO78.VREF'
    character(len=19) :: chamel, chamgd, chamno, chgene, modelLigrel, chs(2)
    character(len=19) :: chel1
    character(len=16) :: nosy, option
    character(len=24) :: chgeom, chharm, chtime
    character(len=19) :: krefe
    character(len=19), parameter :: knume = '&&BAMO78.NUM_RANG', kinst = '&&BAMO78.INSTANT'
    integer(kind=8) :: jnume, jinst
    character(len=8) :: ctype, sdnoli, k8bla, model, materField, crit, mesh, answer
    character(len=1) :: typcoe
    character(len=2) :: codret
    character(len=24) :: trgene
    integer(kind=8) :: jtrgen, tmod(1)
    character(len=24) :: materCode, compor, caraElem
    real(kind=8) :: lcoer(2)
    complex(kind=8) :: lcoec(2)
    aster_logical :: lcumu(2), lcoc(2)
    integer(kind=8) :: iarc2, ievnew, iopt, jvPara, n, nbins2
    integer(kind=8) :: nbtrou, nc, nncp, num0, nume0
    real(kind=8) :: epsi, rundf, time
    real(kind=8), pointer :: base(:) => null()
    integer(kind=8), pointer :: ordr(:) => null()

    integer(kind=8), parameter :: numeHarm = 0
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - INITIALISATIONS
    basemo = ' '
    ctype = 'K24'
    sdnoli = trange(1:8)
    krefe = nomres
    lcoc = ASTER_FALSE
    lcumu = ASTER_FALSE
    lcoer = 1.d0
    lcoec = dcmplx(1.d0, 0.d0)

! --- RECUPERATION BASE MODALE
    call getvid(' ', 'BASE_MODALE', scal=basemo, nbret=ibid)
    call getvid(' ', 'RESU_FINAL', scal=k8bid, nbret=ievnew)

! - Get material fields
    materField = ' '
    call getvid(' ', 'CHAM_MATER', scal=materField, nbret=n1)
    if (n1 .ne. 0) then
        call rcmfmc(materField, materCode, l_ther_=ASTER_FALSE)
    else
        materCode = ' '
    end if

! - Get elementary characteristics
    caraElem = ' '
    call getvid(' ', 'CARA_ELEM', scal=caraElem, nbret=n1)

! - NOMBRE DE MODES
    call rs_get_liststore(basemo, nbmode)
!
! --- NUME_DDL ATTACHE A LA BASE MODALE
!
    call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=numedd)
!
! --- NOUVELLE NUMEROTATION PAS NECESSAIRE ENCORE DANS REST_COND_TRAN
!
    numddl = numedd(1:14)
!
! --- RECOPIE DES MODES PROPRES DANS UN VECTEUR DE TRAVAIL
!
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
    AS_ALLOCATE(vr=base, size=nbmode*neq)
    call copmod(basemo, bmodr=base, numer=numddl)
!
    call dismoi('NOM_MODELE', numddl, 'NUME_DDL', repk=model)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
!
! - No POUX beams
!
    call dismoi('EXI_POUX', modelLigrel, 'LIGREL', repk=answer)
    if (answer .eq. 'OUI') then
        call utmess('F', 'DYNAPOST_1')
    end if
!
! --- CHAMPS SUR LESQUELS ON RESTITUE
!
    call getvtx(' ', 'TOUT_CHAM', nbval=0, nbret=n0)
    if (n0 .ne. 0) then
        nbcham = 3
        champ(1) = 'DEPL'
        champ(2) = 'VITE'
        champ(3) = 'ACCE'
    else
        call getvtx(' ', 'NOM_CHAM', nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            nbcham = -n1
            if (nbcham .gt. 3) then
                ASSERT(.false.)
            end if
            call getvtx(' ', 'NOM_CHAM', nbval=nbcham, vect=champ, nbret=n1)
        else
            call utmess('A', 'DYNAPOST_2')
            goto 999
        end if
    end if

! --- RECUPERATION DES INSTANTS ET DES NUMEROS DE RANGEMENT
    call rstran('NON', trange, ' ', 1, kinst, &
                knume, nbinst, iretou)
    if (iretou .ne. 0) then
        call utmess('F', 'DYNAPOST_3')
    end if
    call jeexin(kinst, iret)
    if (iret .gt. 0) then
        call jeveuo(kinst, 'L', jinst)
        call jeveuo(knume, 'L', jnume)
    end if
    call jeveuo(trange//'.ORDR', 'L', vi=ordr)
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=n)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n)
!
! --- CREATION DE LA SD RESULTAT EVOL_NOLI
!
    nume0 = 0
    if (ievnew .eq. 0) then
        call rscrsd('G', nomres, typres, nbinst)
    else
        call rs_getlast(nomres, nume0)
        call rsorac(nomres, 'INST', ibid, zr(jinst), k8bid, &
                    c16bid, epsi, crit, tmod, 1, &
                    nbtrou)
        nume = tmod(1)
        if (nbtrou .ne. 0) nume0 = nume
        nbins2 = nbinst+nume0
        call rsagsd(nomres, nbins2)
    end if
!
! --- PROJECTION SUR BASE PHYSIQUE
!
    do icham = 1, nbcham
        do iarch = 1, nbinst
            time = zr(jinst+iarch-1)
            num0 = zi(jnume+iarch-1)
            nume = ordr(num0)
            iarc2 = iarch+nume0-1
!
!         --- RECUP POINTEUR SUR CHAMP GENERALISE
!
!
            call rsadpa(sdnoli, 'L', 1, 'TRAN_GENE_NOLI', nume, &
                        1, sjv=jtrgen, styp=ctype)
            trgene = zk24(jtrgen)
!
            if (champ(icham) .eq. 'DEPL') then
                chgene = trgene(1:18)//'D'
            else if (champ(icham) .eq. 'VITE') then
                chgene = trgene(1:18)//'V'
            else if (champ(icham) .eq. 'ACCE') then
                chgene = trgene(1:18)//'A'
            else
                call utmess('A', 'DYNAPOST_2')
                goto 300
            end if
!
            call jeexin(chgene, iret)
            if (iret .eq. 0) then
                call utmess('F', 'DYNAPOST_4')
            else
                call jeveuo(chgene, 'L', jrestr)
            end if
!
!
!         --- RECUP POINTEUR SUR CHAMP PHYSIQUE DANS SD RESULTAT
!
            call rsexch(' ', nomres, champ(icham) (1:4), iarc2, chamno, iret)
!
!         --- CREATION DU CHAMP
            if (iret .eq. 0) call detrsd('CHAM_NO', chamno)
!
            call vtcreb(chamno, 'G', 'R', &
                        nume_ddlz=numedd, &
                        nb_equa_outz=neq)
            call jeveuo(chamno(1:19)//'.VALE', 'E', ldnew)
!
!         --- TRANSFERT EFFECTIF SUR BASE PHYSIQUE
!
            call mdgeph(neq, nbmode, base, zr(jrestr), zr(ldnew))
!
!         --- STOCKAGE CHAMP PHYSIQUE
!
            call rsnoch(nomres, champ(icham) (1:4), iarc2)
            if (icham .eq. 1) then
                call rsadpa(nomres, 'E', 1, 'INST', iarc2, 0, sjv=linst)
                zr(linst) = zr(jinst+iarch-1)
                call rsadpa(nomres, 'E', 1, 'MODELE', iarc2, 0, sjv=jvPara)
                zk8(jvPara) = model
                call rsadpa(nomres, 'E', 1, 'CHAMPMAT', iarc2, 0, sjv=jvPara)
                zk8(jvPara) = materField
                call rsadpa(nomres, 'E', 1, 'CARAELEM', iarc2, 0, sjv=jvPara)
                zk8(jvPara) = caraElem(1:8)
            end if
!
            call jelibe(chgene)
!
        end do
300     continue
    end do
!
! --- ENRICHISSEMENT SD TRAN_GENE -> EVOL_NOLI SD_VERI = 'NON' !!!
!
    if (typres .ne. 'EVOL_NOLI') then
        call refdcp(basemo, krefe(1:8))
        goto 999
    end if
!
    chtime = ' '
    typcoe = ' '
    k8bla = ' '
    rundf = r8vide()
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    compor = materCode(1:8)//'.COMPOR'

! - Get geometry field
    call megeom(model, chgeom)

! - Create field for Fourier
    call meharm(model, numeHarm, chharm)

    do iarch = 1, nbinst
        num0 = zi(jnume+iarch-1)
        nume = ordr(num0)
        time = zr(jinst+iarch-1)
        call mechti(chgeom(1:8), time, rundf, rundf, chtime)
        call vrcins(model, materCode, caraElem, time, chvarc(1:19), &
                    codret)
        call vrcref(model, materCode(1:8), caraElem(1:8), chvref(1:19))
        iarc2 = iarch+nume0-1
!
!         --- RECUP POINTEUR SUR CHAMP PHYSIQUE DANS SD RESULTAT
        do iopt = 1, 2
!
            option = opti(iopt)
            call rsexch(' ', sdnoli, option, nume, chel1, iret)
            call rsexch(' ', nomres, option, iarc2, chamel, iret)
!
            if (iopt .eq. 1) then
                nosy = 'SIEF_ELGA'
                call rsexch(' ', nomres, 'DEPL', iarc2, chamgd, iret)
                call compStress(model, modelLigrel, &
                                materCode, caraElem, compor, &
                                chamgd, chgeom, &
                                chtime, chharm, &
                                chvarc, chvref, &
                                'V', chel2, iret)
                call celces(chel2, 'V', ches2)
                nc = 2
                chs(1) = ches2
                chs(2) = ches1
            end if
            if (iopt .eq. 2) then
                nosy = ' '
                nc = 1
                chs(1) = ches1
            end if
!         --- CREATION DU CHAMP
!
            call celces(chel1, 'V', ches1)
            call cesfus(nc, chs, lcumu, lcoer, lcoec, &
                        lcoc(1), 'V', ches3)
            call cescel(ches3, modelLigrel, nosy, ' ', 'OUI', &
                        nncp, 'G', chamel, 'F', ibid)
!
!         --- STOCKAGE CHAMP PHYSIQUE
!
            call rsnoch(nomres, option, iarc2)
        end do
!
        call rsexch('F', sdnoli, 'COMPORTEMENT', nume, chel1, iret)
        call rsexch(' ', nomres, 'COMPORTEMENT', iarc2, chamel, iret)
        if (iret .eq. 0) call detrsd('CHAMP_GD', chamel)
!
        call dyna_comp_fuse(mesh, chel1, chamel)
!
        call rsnoch(nomres, 'COMPORTEMENT', iarc2)
    end do
!
999 continue

! - MENAGE
    AS_DEALLOCATE(vr=base)
    call jedetr(knume)
    call jedetr(kinst)
!
    call jedema()
end subroutine
