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
subroutine peritr(tablOutZ, &
                  modelZ, numeHarm, &
                  nbFactorKeyword)
!
    implicit none
!
#include "asterfort/calcul.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/chpve2.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gettco.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerazo.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mecact.h"
#include "asterfort/mecham.h"
#include "asterfort/memaxm.h"
#include "asterfort/memoy.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: tablOutZ
    character(len=*), intent(in) :: modelZ
    integer(kind=8), intent(in) :: numeHarm, nbFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "RICE_TRACEY"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = "RICE_TRACEY"
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    integer(kind=8) :: cellNume, long
    integer(kind=8) :: ifm, niv, iFactorKeyword, np, nq, n1, n2, iret, nbret
    integer(kind=8) :: iStore, jvPara, lvale, nt, nm, nc
    integer(kind=8) :: ng, kk, nbgrma, jgr, ig, nbma, jad, nbmail, jma, im, nume, ier
    integer(kind=8) :: numeStore, nbStore, nbMaiT
    integer(kind=8), parameter :: mxvale = 5
    real(kind=8) :: vr(mxvale)
    real(kind=8) :: prec, inst, rsr0, volu, numema, triax, lnrsr0
    real(kind=8) :: rtval(2), valer(3)
    character(len=8) :: k8b, mesh, result, crit, cellName, nommai
    character(len=8) :: valek(2)
    character(len=16) :: resultType, option, optcal(2), toptca(2)
    character(len=19), parameter :: chelem = '&&PERITR.RITR'
    character(len=19), parameter :: varnul = '&&PERITR.VARNUL'
    character(len=19), parameter :: varipr = '&&PERITR.SDRMR'
    character(len=19), parameter :: chRiceOpt = '&&PERITR.CH.SOUSOP'
    character(len=19), parameter :: listStoreJv = '&&PEECIN.NUME_ORDRE'
    integer(kind=8), pointer :: listStore(:) => null()
    character(len=19), parameter :: listTimeJv = '&&PEECIN.INSTANT'
    real(kind=8), pointer :: listTime(:) => null()
    character(len=24) :: chgeom, chharm, ligrel
    character(len=24) :: compor, nomma2
    character(len=24) :: sigm, varipg, varimg, depla, riceOpt
    complex(kind=8) :: c16b
    aster_logical :: lFieldUser, lResultUser
    integer(kind=8), parameter :: nbParaResu = 6, nbParaField = 4
    character(len=16), parameter :: tablRParaName(nbParaResu) = &
                                    (/'NUME_ORDRE      ', 'INST            ', &
                                      'LIEU            ', 'ENTITE          ', &
                                      'TX_CROIS_CAVITES', 'VOLUME_CONCERNE '/)
    character(len=8), parameter :: tablRParaType(nbParaResu) = &
                                   (/'I  ', 'R  ', &
                                     'K24', 'K8 ', &
                                     'R  ', 'R  '/)
    character(len=16), parameter :: tablFParaName(nbParaField) = &
                                    (/'LIEU            ', 'ENTITE          ', &
                                      'TX_CROIS_CAVITES', 'VOLUME_CONCERNE '/)
    character(len=8), parameter :: tablFParaType(nbParaField) = &
                                   (/'K24', 'K8 ', &
                                     'R  ', 'R  '/)
    character(len=8), parameter :: tabcmp(5) = (/'TRIAX ', 'RSR0  ', &
                                                 'VOLU  ', 'NUMEMA', &
                                                 'DEPSEQ'/)
    character(len=16), parameter :: fieldTypePara(3) = &
                                    (/'NOEU#DEPL_R', 'NOEU#TEMP_R', 'ELEM#ENER_R'/)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)

! - Initializations
    c16b = (0.d0, 0.d0)
    inst = 0.d0
    option = factorKeyword
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=ligrel)

! - Get stress from user
    call getvid(' ', 'CHAM_GD', scal=sigm, nbret=nbret)
    lFieldUser = nbret .ne. 0
    if (lFieldUser) then
        call chpve2(sigm, 3, fieldTypePara, ier)
    end if

! - Get parameters
    call getvr8(' ', 'INST', scal=inst, nbret=nbRet)
    call getvid(' ', 'RESULTAT', scal=result, nbret=nbRet)
    lResultUser = nbret .ne. 0

    call getvtx(factorKeyword, 'OPTION', iocc=1, scal=optcal(1), nbret=np)
    call getvtx(factorKeyword, 'LOCAL', iocc=1, scal=optcal(2), nbret=nq)
    if (nbFactorKeyword .gt. 1) then
        do iFactorKeyword = 2, nbFactorKeyword
            call getvtx(factorKeyword, 'OPTION', iocc=iFactorKeyword, scal=toptca(1), nbret=n1)
            call getvtx(factorKeyword, 'LOCAL', iocc=iFactorKeyword, scal=toptca(2), nbret=n2)
            if ((toptca(1) .ne. optcal(1)) .or. (toptca(2) .ne. optcal(2))) then
                call utmess('F', 'UTILITAI3_83')
            end if
        end do
    end if

! - Check and prepare input fields
    call mecham(option, modelZ, numeHarm, &
                chgeom, chharm, iret)
    if (iret .ne. 0) goto 110
    mesh = chgeom(1:8)

! - Create list of time steps and storing
    if (lFieldUser) then
        nbStore = 1
        call wkvect(listStoreJv, 'V V I', nbStore, vi=listStore)
        liststore(1) = 1
        call wkvect(listTimeJv, 'V V R', nbStore, vr=listTime)
        listTime(1) = inst
        call tbcrsd(tablOutZ, 'G')
        call tbajpa(tablOutZ, nbParaField, tablFParaName, tablFParaType)
    else
        call gettco(result, resultType)
        if (resultType(1:9) .ne. 'EVOL_NOLI') then
            call utmess('F', 'UTILITAI3_84')
        end if
        call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
        call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
        call rsutnu(result, ' ', 0, listStoreJv, nbStore, prec, crit, iret)
        if (iret .ne. 0) goto 80
        call jeveuo(listStoreJv, 'L', vi=listStore)

! ----- Create list of time steps
        call wkvect(listTimeJv, 'V V R', nbStore, vr=listTime)

! ----- Get time steps
        call jenonu(jexnom(result//'           .NOVA', 'INST'), iret)
        if (iret .ne. 0) then
            do iStore = 1, nbStore
                numeStore = listStore(iStore)
                call rsadpa(result, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
                listTime(iStore) = zr(jvPara)
            end do
        end if
        call tbcrsd(tablOutZ, 'G')
        call tbajpa(tablOutZ, nbParaResu, tablRParaName, tablRParaType)
    end if

! - Create field
    lnrsr0 = 0.d0
    call mecact('V', varipr, 'MAILLA', mesh, 'NEUT_R', &
                ncmp=1, nomcmp='X1', sr=0.d0)
    call wkvect('&&PERITR.TRAV1', 'V V R', mxvale, lvale)
    do iStore = 1, nbStore
        call jemarq()
        call jerecu('V')

! ----- Current storing index
        numeStore = listStore(iStore)
        inst = listTime(iStore)
        valer(1) = inst
!
        call rsexch(' ', result, 'COMPORTEMENT', numeStore, compor, iret)
        if (lResultUser) then
            call rsexch('F', result, 'SIEF_ELGA', numeStore, sigm, iret)
            call rsexch('F', result, 'VARI_ELGA', numeStore, varipg, iret)
            if (numeStore .ge. 1) then
                call rsexch('F', result, 'VARI_ELGA', numeStore-1, varimg, iret)
            else
                call copisd('CHAMP_GD', 'V', varipg, varnul)
                call jelira(varnul//'.CELV', 'LONUTI', long)
                call jerazo(varnul//'.CELV', long, 1)
            end if
            call rsexch('F', result, 'DEPL', numeStore, depla, iret)
        end if

! ----- Create field for option of Rice Tracey
        riceOpt = optcal(1)//optcal(2) (1:8)
        call mecact('V', chRiceOpt, 'MAILLA', mesh, 'NEUT_K24', &
                    ncmp=1, nomcmp='Z1', sk=riceOpt)

! ----- Add input fields
        lchin(1) = chgeom(1:19)
        lpain(1) = 'PGEOMER'
        lchin(2) = sigm(1:19)
        lpain(2) = 'PCONTPR'
        if (numeStore .ge. 1) then
            lchin(3) = varimg(1:19)
        else
            lchin(3) = varnul
        end if
        lpain(3) = 'PVARIMR'
        lchin(4) = varipg(1:19)
        lpain(4) = 'PVARIPR'
        lchin(5) = varipr
        lpain(5) = 'PSDRMR'
        lchin(6) = chRiceOpt
        lpain(6) = 'PSOUSOP'
        lchin(7) = compor(1:19)
        lpain(7) = 'PCOMPOR'
        nbFieldIn = 7

! ----- Add output fields
        lchout(1) = chelem
        lpaout(1) = 'PRICTRA'
        lchout(2) = '&&PERITR.SDRPR'
        lpaout(2) = 'PSDRPR'
        call calcul('S', option, ligrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    'V', 'OUI')
!
        do iFactorKeyword = 1, nbFactorKeyword
            call getvtx(option(1:11), 'TOUT', iocc=iFactorKeyword, nbval=0, nbret=nt)
            call getvem(mesh, 'MAILLE', option(1:11), 'MAILLE', iFactorKeyword, &
                        0, k8b, nm)
            call getvem(mesh, 'GROUP_MA', option(1:11), 'GROUP_MA', iFactorKeyword, &
                        0, k8b, ng)
            if (nt .ne. 0) then
                if (optcal(2) .eq. 'OUI') then
                    call memaxm('MAX', chelem, 'RSR0', mxvale, tabcmp, vr, 0, [0])
                    do kk = 1, mxvale
                        zr(lvale+kk-1) = vr(kk)
                    end do
                else if (optcal(2) .eq. 'NON') then
                    call memoy(chelem, 1, chelem, 3, vr, &
                               0, [0])
                    zr(lvale) = vr(1)
                    zr(lvale+2) = vr(2)
                    triax = zr(lvale)
                    call memoy(chelem, 5, chelem, 3, vr, &
                               0, [0])
                    zr(lvale+4) = vr(1)
                    lnrsr0 = lnrsr0+0.283d0*sign(1.d0, triax)*exp(1.5d0*abs(triax))*zr(lvale+4)
                    zr(lvale+1) = exp(lnrsr0)
                    zr(lvale+3) = 0.d0
                end if
                rsr0 = zr(lvale+1)
                volu = zr(lvale+2)
                numema = zr(lvale+3)
                if (optcal(2) .eq. 'OUI') then
                    cellNume = nint(numema)
                    cellName = int_to_char8(cellNume)
                    valek(1) = cellName
                    valek(2) = 'MAILLE'
                else
                    valek(1) = mesh
                    valek(2) = 'TOUT'
                end if
                rtval(1) = rsr0
                rtval(2) = volu
                if (lResultUser) then
                    valer(2) = rtval(1)
                    valer(3) = rtval(2)
                    call tbajli(tablOutZ, nbParaResu, tablRParaName, [numeStore], valer, &
                                [c16b], valek, 0)
                else
                    call tbajli(tablOutZ, nbParaField, tablFParaName, [numeStore], rtval, &
                                [c16b], valek, 0)
                end if
            end if
            if (ng .ne. 0) then
                nbgrma = -ng
                call wkvect('&&PERITR_GROUPM', 'V V K24', nbgrma, jgr)
                call getvem(mesh, 'GROUP_MA', option(1:11), 'GROUP_MA', iFactorKeyword, &
                            nbgrma, zk24(jgr), ng)
                do ig = 1, nbgrma
                    nomma2 = zk24(jgr+ig-1)
                    call jeexin(jexnom(mesh//'.GROUPEMA', nomma2), iret)
                    if (iret .eq. 0) then
                        call utmess('A', 'UTILITAI3_46', sk=nomma2)
                        goto 50
                    end if
                    call jelira(jexnom(mesh//'.GROUPEMA', nomma2), 'LONUTI', nbma)
                    if (nbma .eq. 0) then
                        call utmess('A', 'UTILITAI3_47', sk=nomma2)
                        goto 50
                    end if
                    call jeveuo(jexnom(mesh//'.GROUPEMA', nomma2), 'L', jad)
                    if (optcal(2) .eq. 'OUI') then
                        call memaxm('MAX', chelem, 'RSR0', mxvale, tabcmp, &
                                    vr, nbma, zi(jad))
                        do kk = 1, mxvale
                            zr(lvale+kk-1) = vr(kk)
                        end do
                    else if (optcal(2) .eq. 'NON') then
                        call memoy(chelem, 1, chelem, 3, vr, &
                                   nbma, zi(jad))
                        zr(lvale) = vr(1)
                        zr(lvale+2) = vr(2)
                        triax = zr(lvale)
                        call memoy(chelem, 5, chelem, 3, vr, &
                                   nbma, zi(jad))
                        zr(lvale+4) = vr(1)
                        lnrsr0 = lnrsr0+0.283d0*sign(1.d0, triax)*exp(1.5d0*abs(triax))*zr(lval&
                                 &e+4)
                        zr(lvale+1) = exp(lnrsr0)
                        zr(lvale+3) = 0.d0
                    end if
                    rsr0 = zr(lvale+1)
                    volu = zr(lvale+2)
                    numema = zr(lvale+3)
                    if (optcal(2) .eq. 'OUI') then
                        cellNume = nint(numema)
                        cellName = int_to_char8(cellNume)
                        valek(1) = cellName
                        valek(2) = 'MAILLE'
                    else
                        valek(1) = mesh
                        valek(2) = 'TOUT'
                    end if
                    rtval(1) = rsr0
                    rtval(2) = volu
                    if (lResultUser) then
                        valer(2) = rtval(1)
                        valer(3) = rtval(2)
                        call tbajli(tablOutZ, nbParaResu, tablRParaName, [numeStore], valer, &
                                    [c16b], valek, 0)
                    else
                        call tbajli(tablOutZ, nbParaField, tablFParaName, [numeStore], rtval, &
                                    [c16b], valek, 0)
                    end if
50                  continue
                end do
                call jedetr('&&PERITR_GROUPM')
            end if
            if (nm .ne. 0) then
                nbmail = -nm
                call wkvect('&&PERITR_MAILLE', 'V V K8', nbmail, jma)
                call getvem(mesh, 'MAILLE', option(1:11), 'MAILLE', iFactorKeyword, &
                            nbmail, zk8(jma), nm)
                call jelira(mesh//'.TYPMAIL', 'LONMAX', nbMaiT)
                do im = 1, nbmail
                    nommai = zk8(jma+im-1)
                    nume = char8_to_int(nommai)
                    if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                        call utmess('A', 'UTILITAI3_49', sk=zk8(jma+im-1))
                        goto 70
                    end if
                    if (optcal(2) .eq. 'OUI') then
                        call memaxm('MAX', chelem, 'RSR0', mxvale, tabcmp, &
                                    vr, 1, [nume])
                        do kk = 1, mxvale
                            zr(lvale+kk-1) = vr(kk)
                        end do
                    else if (optcal(2) .eq. 'NON') then
                        call memoy(chelem, 1, chelem, 3, vr, &
                                   1, [nume])
                        zr(lvale) = vr(1)
                        zr(lvale+2) = vr(2)
                        triax = zr(lvale)
                        call memoy(chelem, 5, chelem, 3, vr, &
                                   1, [nume])
                        zr(lvale+4) = vr(1)
                        lnrsr0 = lnrsr0+0.283d0*sign(1.d0, triax)*exp(1.5d0*abs(triax))*zr(lval&
                                 &e+4)
                        zr(lvale+1) = exp(lnrsr0)
                        zr(lvale+3) = 0.d0
                    end if
                    rsr0 = zr(lvale+1)
                    volu = zr(lvale+2)
                    numema = zr(lvale+3)
                    if (optcal(2) .eq. 'OUI') then
                        cellNume = nint(numema)
                        cellName = int_to_char8(cellNume)
                        valek(1) = cellName
                        valek(2) = 'MAILLE'
                    else
                        valek(1) = mesh
                        valek(2) = 'TOUT'
                    end if
                    rtval(1) = rsr0
                    rtval(2) = volu
                    if (lResultUser) then
                        valer(2) = rtval(1)
                        valer(3) = rtval(2)
                        call tbajli(tablOutZ, nbParaResu, tablRParaName, [numeStore], valer, &
                                    [c16b], valek, 0)
                    else
                        call tbajli(tablOutZ, nbParaField, tablFParaName, [numeStore], rtval, &
                                    [c16b], valek, 0)
                    end if
70                  continue
                end do
                call jedetr('&&PERITR_MAILLE')
            end if
        end do
        call copisd('CHAMP_GD', 'V', '&&PERITR.SDRPR', varipr)
        call detrsd('CARTE', chRiceOpt)
        call detrsd('CHAM_ELEM', chelem)
        call jedema()
    end do
!
80  continue
!
! - MENAGE
    call jedetr(listStoreJv)
    call jedetr(listTimeJv)
    call jedetr('&&PERITR.TRAV1')
    call detrsd('CHAMP_GD', varnul)
    call detrsd('CHAMP_GD', '&&PERITR.SDRPR')
    call detrsd('CHAMP_GD', varipr)
!
110 continue
    call jedema()
end subroutine
