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
! aslint: disable=W0413
!
subroutine peweib(resu, modele, mate, mateco, cara, chmat, &
                  nh, nbocc, iresu, nomcmd)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/calcul.h"
#include "asterfort/chmrck.h"
#include "asterfort/chpve2.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/exlim3.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mecact.h"
#include "asterfort/mecham.h"
#include "asterfort/mesomm.h"
#include "asterfort/rccome.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: iresu, nh, nbocc
    character(len=*) :: resu, modele, mate, mateco, cara, nomcmd
!     OPERATEUR   POST_ELEM
!     ( TRAITEMENT DU MOT CLE-FACTEUR "WEIBULL" )
!     OPERATEUR   RECA_WEIBULL
!     ( RECALAGE DE LA METHODE DE WEIBULL )
!     ------------------------------------------------------------------
!
    real(kind=8) :: valr(3)
    integer(kind=8) :: nbparr, nbpard, nbmtcm, inum, nbout
    integer(kind=8) :: vali
    integer(kind=8) :: ibid, ibik, mxvale, ifm, niv
    integer(kind=8) :: nd, ng, ni, nm, np, nq, nr, nt, i, n1, n2, n3
    integer(kind=8) :: no, nlo, nli
    integer(kind=8) :: iret, nbordr, jord, jins, nc, nbgrma, jgr, ig, nbma, jad
    integer(kind=8) :: nbmail, jma, im, nume, imc, ier
    integer(kind=8) :: numord, iainst, iord, nbmtrc, nbin, iocc, nbMaiT
    parameter(mxvale=3, nbparr=7, nbpard=5)
    real(kind=8) :: rtval(mxvale), prec, inst, valer(4), vref, coesym, mref, sref, probaw
    real(kind=8) :: sigmaw
    character(len=1) :: base
    character(len=2) :: codret
    character(len=11) :: k11
    character(len=8) :: k8b, noma, resul, crit, chmat, nommai, typarr(nbparr), typard(nbpard)
    character(len=8) :: lpain(9), lpaout(2), valek(2)
    character(len=16) :: typres, option, optio2, optcal(2), toptca(2), nomrc, noparr(nbparr)
    character(len=16) :: nopard(nbpard), motcl1, motcl2, motcl3
    character(len=19) :: chelem, knum, kins, tabtyp(3), chvarc
    character(len=24) :: chgeom, chharm
    character(len=24) :: valk(2), nomgrm
    character(len=24) :: mlggma, ligrel, lchin(9), compor
    character(len=24) :: lchout(2), contg, defog, varig, ssoup
    character(len=24) :: kvalrc, kvalrk, vale2(2)
    aster_logical :: opti
    complex(kind=8) :: c16b
    character(len=8), pointer :: l_nom_mat(:) => null()
    real(kind=8), pointer :: trav1(:) => null()
!
    data noparr/'NUME_ORDRE', 'INST', 'LIEU', 'ENTITE', 'SIGMA_WEIBULL',&
     &     'PROBA_WEIBULL', 'SIGMA_WEIBULL**M'/
    data typarr/'I', 'R', 'K24', 'K8', 'R', 'R', 'R'/
    data nopard/'LIEU', 'ENTITE', 'SIGMA_WEIBULL', 'PROBA_WEIBULL',&
     &     'SIGMA_WEIBULL**M'/
    data typard/'K8', 'K8', 'R', 'R', 'R'/
    data tabtyp/'NOEU#DEPL_R', 'NOEU#TEMP_R', 'ELEM#ENER_R'/
    data chvarc/'&&PEWEIB.CHVARC'/
!     ------------------------------------------------------------------
    call jemarq()
    c16b = (0.d0, 0.d0)
!
! --- RECUPERATION DU NIVEAU D'IMPRESSION
    call infniv(ifm, niv)
!
    if (nomcmd(1:12) .eq. 'RECA_WEIBULL') then
        opti = .true.
        motcl1 = ' '
        motcl2 = 'RESU'
        motcl3 = 'RESU'
        base = 'V'
    else if (nomcmd(1:9) .eq. 'POST_ELEM') then
        opti = .false.
        motcl1 = 'WEIBULL'
        motcl2 = ' '
        motcl3 = 'WEIBULL'
        base = 'G'
    end if
!
    inst = 0.d0
!
    nd = 0
    if (.not. opti) then
        call getvid(' ', 'CHAM_GD', scal=contg, nbret=nd)
    end if
    if (nd .ne. 0) then
        call chpve2(contg, 3, tabtyp, ier)
    end if
    ni = 0
    no = 0
    nlo = 0
    nli = 0
    if (.not. opti) then
        call getvr8(' ', 'NUME_ORDRE', scal=inst, nbret=no)
        call getvr8(' ', 'LIST_ORDRE', scal=inst, nbret=nlo)
        call getvr8(' ', 'LIST_INST', scal=inst, nbret=nli)
        inst = 0.d0
        call getvr8(' ', 'INST', scal=inst, nbret=ni)
        if ((ni .ne. 0) .or. (no .ne. 0) .or. (nli .ne. 0) .or. (nlo .ne. 0)) then
            call utmess('I', 'UTILITAI6_63')
        end if
    end if
!
    if (.not. opti) then
        call getvid(motcl2, 'RESULTAT', iocc=1, scal=resul, nbret=nr)
    else
        call getvid(motcl2, 'EVOL_NOLI', iocc=iresu, scal=resul, nbret=nr)
    end if
    call getvtx(motcl1, 'OPTION', iocc=1, scal=optcal(1), nbret=np)
    call getvtx(motcl1, 'CORR_PLAST', iocc=1, scal=optcal(2), nbret=nq)
    if (nbocc .gt. 1) then
        do i = 2, nbocc
            call getvtx(motcl1, 'OPTION', iocc=i, scal=toptca(1), nbret=n1)
            call getvtx(motcl1, 'CORR_PLAST', iocc=i, scal=toptca(2), nbret=n2)
            if ((toptca(1) .ne. optcal(1)) .or. (toptca(2) .ne. optcal(2))) then
                call utmess('F', 'UTILITAI3_83')
            end if
        end do
    end if
!
    option = 'WEIBULL'
    call mecham(option, modele, nh, &
                chgeom, chharm, iret)
    if (iret .ne. 0) goto 100
    noma = chgeom(1:8)
    mlggma = noma//'.GROUPEMA'
!
    call exlim3(motcl3, 'V', modele, ligrel)
!
!     CREATION CARTE CONSTANTE ET NULLE SUR TOUT LE MAILLAGE
    call mecact('V', '&&PEWEIB.SIGIE', 'MAILLA', noma, 'DOMA_R', &
                ncmp=1, nomcmp='DOMA', sr=0.d0)
!
    knum = '&&PEWEIB.NUME_ORDRE'
    kins = '&&PEWEIB.INSTANT'
    if (nd .ne. 0) then
        nbordr = 1
        call wkvect(knum, 'V V I', nbordr, jord)
        zi(jord) = 1
        call wkvect(kins, 'V V R', nbordr, jins)
        zr(jins) = inst
        call tbcrsd(resu, base)
        call tbajpa(resu, nbpard, nopard, typard)
    else
        call gettco(resul, typres)
        if (typres(1:9) .ne. 'EVOL_NOLI') then
            call utmess('F', 'UTILITAI3_84')
        end if
!
        np = 0
        if (.not. opti) then
            call getvr8(motcl2, 'PRECISION', scal=prec, nbret=np)
        else
            call getvr8(motcl2, 'PRECISION', iocc=iresu, scal=prec, nbret=np)
        end if
        nc = 0
        if (.not. opti) then
            call getvtx(motcl2, 'CRITERE', scal=crit, nbret=nc)
        else
            call getvtx(motcl2, 'CRITERE', iocc=iresu, scal=crit, nbret=nc)
        end if
!
        if (.not. opti) then
            call rsutnu(resul, motcl2, 1, knum, nbordr, &
                        prec, crit, iret)
        else
            call rsutnu(resul, motcl2, iresu, knum, nbordr, &
                        prec, crit, iret)
        end if
!
        if (iret .ne. 0) goto 90
        call jeveuo(knum, 'L', jord)
!        --- ON RECUPERE LES INSTANTS ---
        call wkvect(kins, 'V V R', nbordr, jins)
        call jenonu(jexnom(resul//'           .NOVA', 'INST'), iret)
        if (iret .ne. 0) then
            do iord = 1, nbordr
                numord = zi(jord+iord-1)
                call rsadpa(resul, 'L', 1, 'INST', numord, &
                            0, sjv=iainst, styp=k8b)
                zr(jins+iord-1) = zr(iainst)
            end do
        end if
        call tbcrsd(resu, base)
        call tbajpa(resu, nbparr, noparr, typarr)
    end if
!
!     --- VERIF D'HOMOGENEITE WEIBULL ---
!
    if (.not. opti) then
        call getvid(motcl2, 'CHAM_MATER', iocc=iresu, scal=chmat, nbret=n3)
        if (n3 .eq. 0) chmat = mate(1:8)
    end if
    call jelira(chmat//'.CHAMP_MAT .VALE', 'LONMAX', nbmtcm)
    AS_ALLOCATE(vk8=l_nom_mat, size=nbmtcm)
    nomrc = 'WEIBULL         '
    call chmrck(chmat, nomrc, l_nom_mat, nbmtrc)
    if (nbmtrc .gt. 1) then
        vali = nbmtrc
        call utmess('A', 'UTILITAI6_60', sk=l_nom_mat(1), si=vali)
    else if (nbmtrc .eq. 0) then
        valk(1) = nomrc
        valk(2) = chmat
        call utmess('F', 'UTILITAI6_59', nk=2, valk=valk)
    end if
!
!     --- RECUPERATION DES PARAMETRES DE LA RC WEIBULL ---
    kvalrc(1:8) = l_nom_mat(1)
    call rccome(kvalrc(1:8), 'WEIBULL', iret, k11_ind_nomrc=k11)
    kvalrc = l_nom_mat(1)//k11//'.VALR'
    kvalrk = l_nom_mat(1)//k11//'.VALK'
    call jeveuo(kvalrc, 'L', ibid)
    call jeveuo(kvalrk, 'L', ibik)
    call jelira(kvalrk, 'LONMAX', imc)
    sref = 0.d0
    do i = 1, imc
        if (zk16(ibik+i-1) .eq. 'SIGM_CNV') sref = zr(ibid+i-1)
        if (zk16(ibik+i-1) .eq. 'M       ') mref = zr(ibid+i-1)
        if (zk16(ibik+i-1) .eq. 'VOLU_REFE') vref = zr(ibid+i-1)
    end do
! CAS WEIBULL_FO
    if (sref .eq. 0.d0) then
        do i = 1, imc
            if (zk16(ibik+i-1) .eq. 'SIGM_REFE') sref = zr(ibid+i-1)
        end do
        valr(1) = mref
        valr(2) = vref
        valr(3) = sref
        call utmess('I', 'UTILITAI6_61', nr=3, valr=valr)
! CAS WEIBULL
    else
        valr(1) = mref
        valr(2) = vref
        valr(3) = sref
        call utmess('I', 'UTILITAI6_62', nr=3, valr=valr)
    end if
!
    AS_ALLOCATE(vr=trav1, size=mxvale)
    do iord = 1, nbordr
        call jemarq()
        call jerecu('V')
        numord = zi(jord+iord-1)
        inst = zr(jins+iord-1)
        valer(1) = inst
!
        call rsexch(' ', resul, 'COMPORTEMENT', numord, compor, &
                    iret)
        if (nr .ne. 0) then
            call rsexch('F', resul, 'SIEF_ELGA', numord, contg, &
                        iret)
            call rsexch('F', resul, 'VARI_ELGA', numord, varig, &
                        iret)
        end if
!
!        --- DANS LE CAS D'UNE OPTION AVEC CORRECTION DE DEFORMATION
!            RECUPERATION DES DEFORMATIONS DE GREEN LAGRANGE ---
!
        if (optcal(2) .eq. 'OUI') then
            call rsexch('F', resul, 'EPSG_ELGA', numord, defog, &
                        iret)
        else
            defog = '&&PEWEIB.EPSG'
        end if
!
!        --- RECUPERATION DU CHAMP DE TEMPERATURE
!            UTILE POUR LE CAS OU SIGU(T)
        call vrcins(modele, mate, cara, inst, chvarc, &
                    codret)
!
!        --- AFFECTATION D'UNE CARTE CONSTANTE SUR LE MAILLAGE :
!            OPTION DE CALCUL WEIBULL ---
!
        ssoup = optcal(1)//optcal(2) (1:8)
        call mecact('V', '&&PEWEIB.CH.SOUSOP', 'MAILLA', noma, 'NEUT_K24', &
                    ncmp=1, nomcmp='Z1', sk=ssoup)
!
        optio2 = 'WEIBULL'
        chelem = '&&PEWEIB.WEIBULL'
        nbin = 9
        lchin(1) = chgeom
        lpain(1) = 'PGEOMER'
        lchin(2) = contg
        lpain(2) = 'PCONTRG'
        lchin(3) = varig
        lpain(3) = 'PVARIPG'
        lchin(4) = defog
        lpain(4) = 'PDEFORR'
        lchin(5) = mateco
        lpain(5) = 'PMATERC'
        lchin(6) = '&&PEWEIB.CH.SOUSOP'
        lpain(6) = 'PSOUSOP'
!  EN ENTREE : CHELEM DES SIGI MAX ATTEINTE AU COURS DU TEMPS
        lchin(7) = '&&PEWEIB.SIGIE'
        lpain(7) = 'PDOMMAG'
!  EN ENTREE : CHELEM DES VARIABLES DE COMMANDES
        lchin(8) = chvarc
        lpain(8) = 'PVARCPR'
        lchin(9) = compor
        lpain(9) = 'PCOMPOR'
        nbout = 2
        lchout(1) = chelem
        lpaout(1) = 'PWEIBUL'
        lchout(2) = '&&PEWEIB.SIGIS'
        lpaout(2) = 'PSIGISG'
        call calcul('S', optio2, ligrel, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, 'V', &
                    'OUI')
!
!        RECOPIE DE SIGIS DANS SIGIE
        call copisd('CHAMP_GD', 'V', '&&PEWEIB.SIGIS', '&&PEWEIB.SIGIE')
!
        do iocc = 1, nbocc
            if (.not. opti) then
                inum = iocc
            else
                inum = iresu
            end if
            call getvtx(motcl3, 'TOUT', iocc=inum, nbval=0, nbret=nt)
            call getvem(noma, 'MAILLE', motcl3, 'MAILLE', inum, &
                        0, k8b, nm)
            call getvem(noma, 'GROUP_MA', motcl3, 'GROUP_MA', inum, &
                        0, k8b, ng)
            call getvr8(motcl3, 'COEF_MULT', iocc=inum, scal=coesym, nbret=n1)
!
            if (nt .ne. 0) then
                call mesomm(chelem, mxvale, vr=trav1)
                probaw = coesym*trav1(1)
                sigmaw = probaw*(sref**mref)
                probaw = 1.0d0-exp(-probaw)
                rtval(3) = sigmaw
                rtval(2) = probaw
                rtval(1) = sigmaw**(1.0d0/mref)
                valek(1) = noma
                valek(2) = 'TOUT'
                if (nr .ne. 0) then
                    valer(2) = rtval(1)
                    valer(3) = rtval(2)
                    valer(4) = rtval(3)
                    call tbajli(resu, nbparr, noparr, [numord], valer, &
                                [c16b], valek, 0)
                else
                    call tbajli(resu, nbpard, nopard, [numord], rtval, &
                                [c16b], valek, 0)
                end if
            end if
!
            if (ng .ne. 0) then
                nbgrma = -ng
                call wkvect('&&PEWEIB_GROUPM', 'V V K24', nbgrma, jgr)
                call getvem(noma, 'GROUP_MA', motcl3, 'GROUP_MA', inum, &
                            nbgrma, zk24(jgr), ng)
                vale2(2) = 'GROUP_MA'
                do ig = 1, nbgrma
                    nomgrm = zk24(jgr+ig-1)
                    call jeexin(jexnom(mlggma, nomgrm), iret)
                    if (iret .eq. 0) then
                        call utmess('A', 'UTILITAI3_46', sk=nomgrm)
                        goto 50
                    end if
                    call jelira(jexnom(mlggma, nomgrm), 'LONUTI', nbma)
                    if (nbma .eq. 0) then
                        call utmess('A', 'UTILITAI3_47', sk=nomgrm)
                        goto 50
                    end if
                    call jeveuo(jexnom(mlggma, nomgrm), 'L', jad)
                    call mesomm(chelem, mxvale, vr=trav1, nbma=nbma, linuma=zi(jad))
                    sigmaw = coesym*trav1(1)*(sref**mref)
                    probaw = sigmaw/(sref**mref)
                    probaw = 1.0d0-exp(-probaw)
                    rtval(3) = sigmaw
                    rtval(2) = probaw
                    rtval(1) = sigmaw**(1.0d0/mref)
                    vale2(1) = nomgrm
                    if (nr .ne. 0) then
                        valer(2) = rtval(1)
                        valer(3) = rtval(2)
                        valer(4) = rtval(3)
                        call tbajli(resu, nbparr, noparr, [numord], valer, &
                                    [c16b], vale2, 0)
                    else
                        call tbajli(resu, nbpard, nopard, [numord], rtval, &
                                    [c16b], vale2, 0)
                    end if
50                  continue
                end do
                call jedetr('&&PEWEIB_GROUPM')
            end if
!
            if (nm .ne. 0) then
                nbmail = -nm
                call wkvect('&&PEWEIB_MAILLE', 'V V K8', nbmail, jma)
                call getvem(noma, 'MAILLE', motcl3, 'MAILLE', inum, &
                            nbmail, zk8(jma), nm)
                valek(2) = 'MAILLE'
                call jelira(noma//'.TYPMAIL', 'LONMAX', nbMaiT)
                do im = 1, nbmail
                    nommai = zk8(jma+im-1)
                    nume = char8_to_int(nommai)
                    if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                        call utmess('A', 'UTILITAI3_49', sk=zk8(jma+im-1))
                        goto 60
                    end if
                    call mesomm(chelem, mxvale, vr=trav1, nbma=1, linuma=[nume])
                    probaw = coesym*trav1(1)
                    sigmaw = probaw*(sref**mref)
                    probaw = 1.0d0-exp(-probaw)
                    rtval(3) = sigmaw
                    rtval(2) = probaw
                    rtval(1) = sigmaw**(1.0d0/mref)
                    valek(1) = nommai
                    if (nr .ne. 0) then
                        valer(2) = rtval(1)
                        valer(3) = rtval(2)
                        valer(4) = rtval(3)
                        call tbajli(resu, nbparr, noparr, [numord], valer, &
                                    [c16b], valek, 0)
                    else
                        call tbajli(resu, nbpard, nopard, [numord], rtval, &
                                    [c16b], valek, 0)
                    end if
60                  continue
                end do
                call jedetr('&&PEWEIB_MAILLE')
            end if
        end do
!
        call detrsd('CHAMP_GD', '&&PEWEIB.EPSG')
        call detrsd('CARTE', '&&PEWEIB.CH.SOUSOP')
        call jedema()
    end do
! FIN BOUCLE SUR LES NUMEROS D ORDRE
90  continue
!
! -- MENAGE
    call jedetr('&&PEWEIB.NUME_ORDRE')
    call jedetr('&&PEWEIB.INSTANT')
    call detrsd('CHAMP_GD', '&&PEWEIB.SIGIE')
    call detrsd('CHAMP_GD', '&&PEWEIB.SIGIS')
    AS_DEALLOCATE(vr=trav1)
    AS_DEALLOCATE(vk8=l_nom_mat)
!
100 continue
!
    call jedema()
end subroutine
