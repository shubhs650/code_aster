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
subroutine peecin(tablOutZ, &
                  modelZ, materFieldZ, materCodeZ, caraElemZ, &
                  numeHarm, nbFactorKeyword)
!
    implicit none
!
#include "asterc/r8depi.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/chpve2.h"
#include "asterfort/compEnergyKinetic.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim3.h"
#include "asterfort/gettco.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
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
#include "asterfort/mechti.h"
#include "asterfort/meharm.h"
#include "asterfort/peenca.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: tablOutZ
    character(len=*), intent(in) :: modelZ, materFieldZ, materCodeZ, caraElemZ
    integer(kind=8), intent(in) :: numeHarm, nbFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "ENERCIN"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'ENER_CIN'
    integer(kind=8) :: iret, np, nc, jad, nbStore, iStore, numeStore, jvPara, jnmo, ibid
    integer(kind=8) :: ie, nt, nm, ng, nbgrma, ig, jgr, nbma, nume, im, lfreq
    integer(kind=8) :: iFactorKeyword, jma, massDiagIndx, ier, nbMaiT, nbret
    integer(kind=8), parameter :: nbpaep = 2
    real(kind=8) :: varpep(nbpaep)
    real(kind=8) :: prec, xfreq, valer(3), inst
    real(kind=8) :: rundf
    character(len=1), parameter :: jvBase = "V"
    character(len=2) :: codret
    character(len=8) :: k8b, mesh, result, crit, nommai, nommas
    character(len=8) :: valk(2), physQuanName
    character(len=16) :: resultType, optmas
    character(len=19) :: field, ligrel
    character(len=19) :: ecinElemUser, ecinElem
    character(len=19) :: chdisp, chvite
    character(len=24) :: fieldType, chtime, chgeom
    character(len=24), parameter :: chmasd = '&&PEECIN.MASD', chfreq = '&&PEECIN.OMEGA2'
    character(len=24) :: chtemp, opt, chharm, nomgrm, valk2(2)
    aster_logical :: l_modal
    complex(kind=8) :: c16b
    character(len=19), parameter :: chvarc = '&&PEECIN.VARC'
    character(len=19), parameter :: listStoreJv = '&&PEECIN.NUME_ORDRE'
    integer(kind=8), pointer :: listStore(:) => null()
    character(len=19), parameter :: listTimeJv = '&&PEECIN.INSTANT'
    real(kind=8), pointer :: listTime(:) => null()
    aster_logical :: lFieldUser, lResultUser, lHasFreq, lHasTime
    integer(kind=8), parameter :: nbParaResu = 6, nbParaField = 4
    character(len=16) :: tablRParaName(nbParaResu)
    character(len=8), parameter :: tablRParaType(nbParaResu) = &
                                   (/'I  ', 'R  ', &
                                     'K24', 'K8 ', 'R  ', 'R  '/)
    character(len=16), parameter :: tablFParaName(nbParaField) = &
                                    (/'LIEU      ', 'ENTITE    ', 'TOTALE    ', 'POUR_CENT '/)
    character(len=8), parameter :: tablFParaType(nbParaField) = &
                                   (/'K24', 'K8 ', 'R  ', 'R  '/)
    character(len=16), parameter :: fieldTypePara(3) = &
                                    (/'NOEU#DEPL_R', 'NOEU#TEMP_R', 'ELEM#ENER_R'/)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    c16b = (0.d0, 0.d0)
    rundf = r8vide()
    lHasTime = ASTER_FALSE
    inst = 0.d0
    chdisp = ' '
    chvite = ' '
    chtemp = ' '
    resultType = ' '

! - Get field from user
    call getvid(' ', 'CHAM_GD', scal=field, nbret=nbret)
    lFieldUser = nbret .ne. 0
    if (lFieldUser) then
        call chpve2(field, 3, fieldTypePara, ier)
        call dismoi('TYPE_SUPERVIS', field, 'CHAMP', repk=fieldType)
        call dismoi('NOM_GD', field, 'CHAMP', repk=physQuanName)
    end if

! - Get parameters
    call getvr8(' ', 'FREQ', scal=xfreq, nbret=nbRet)
    lHasFreq = nbRet .ne. 0
    call getvr8(' ', 'INST', scal=inst, nbret=nbRet)
    lHasTime = nbRet .ne. 0
    call getvid(' ', 'RESULTAT', scal=result, nbret=nbRet)
    lResultUser = nbret .ne. 0

! - Parameters in output table
    tablRParaName(1) = 'NUME_ORDRE'
    tablRParaName(3) = 'LIEU      '
    tablRParaName(4) = 'ENTITE    '
    tablRParaName(5) = 'TOTALE    '
    tablRParaName(6) = 'POUR_CENT '
    if (lResultUser) then
        call gettco(result, resultType)
        if (resultType(1:9) .eq. 'MODE_MECA') then
            tablRParaName(2) = 'FREQ'
        else if (resultType(1:9) .eq. 'EVOL_THER' .or. resultType(1:9) .eq. 'EVOL_ELAS' .or. &
                 resultType(1:9) .eq. 'EVOL_NOLI' .or. resultType(1:10) .eq. 'DYNA_TRANS') then
            tablRParaName(2) = 'INST'
        else
            ASSERT(ASTER_FALSE)
        end if
    end if

! - Check and prepare input fields
    call mecham(option, modelZ, numeHarm, &
                chgeom, chharm, iret)
    if (iret .ne. 0) goto 90
    mesh = chgeom(1:8)
!
    call exlim3(option, 'V', modelZ, ligrel)

! - Create list of time steps and storing
    if (lFieldUser) then
        if (lHasFreq) then
            call utmess('I', 'UTILITAI3_70')
            xfreq = (r8depi()*xfreq)**2
        else
            xfreq = 1.d0
            call utmess('I', 'UTILITAI3_69')
        end if
        nbStore = 1
        call wkvect(listStoreJv, 'V V I', nbStore, vi=listStore)
        liststore(1) = 1
        call wkvect(listTimeJv, 'V V R', nbStore, vr=listTime)
        listTime(1) = inst
        call tbcrsd(tablOutZ, 'G')
        call tbajpa(tablOutZ, nbParaField, tablFParaName, tablFParaType)
    else
        ASSERT(lResultUser)
        call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
        call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
        call rsutnu(result, ' ', 0, listStoreJv, nbStore, prec, crit, iret)
        if (iret .ne. 0) goto 80
        call jeveuo(listStoreJv, 'L', vi=listStore)

! ----- Create list of time steps
        call wkvect(listTimeJv, 'V V R', nbStore, vr=listTime)

! ----- Get frequencies
        call jenonu(jexnom(result//'           .NOVA', 'FREQ'), iret)
        if (iret .ne. 0) then
            do iStore = 1, nbStore
                numeStore = listStore(iStore)
                call rsadpa(result, 'L', 1, 'FREQ', numeStore, 0, sjv=jvPara, istop=0)
                listTime(iStore) = zr(jvPara)
            end do
        end if

! ----- Get time steps
        call jenonu(jexnom(result//'           .NOVA', 'INST'), iret)
        if (iret .ne. 0) then
            lHasTime = ASTER_TRUE
            do iStore = 1, nbStore
                numeStore = listStore(iStore)
                call rsadpa(result, 'L', 1, 'INST', numeStore, 0, sjv=jvPara, istop=0)
                listTime(iStore) = zr(jvPara)
            end do
        end if
        call tbcrsd(tablOutZ, 'G')
        call tbajpa(tablOutZ, nbParaResu, tablRParaName, tablRParaType)
    end if

! - Get option of mass
    massDiagIndx = 1
    if (lResultUser) then
        if (resultType(1:9) .ne. 'EVOL_NOLI') then
            call dismoi('REF_MASS_PREM', result, 'RESU_DYNA', repk=nommas, arret='C')
            if (nommas .ne. ' ') then
                call dismoi('SUR_OPTION', nommas, 'MATR_ASSE', repk=opt, arret='C', ier=ie)
                if (ie .ne. 0) then
                    call utmess('A', 'UTILITAI3_71')
                else
                    if (opt(1:14) .eq. 'MASS_MECA_DIAG') then
                        massDiagIndx = 0
                    end if
                end if
            end if
        end if
        call getvtx(option(1:9), 'OPTION', iocc=1, scal=optmas, nbret=nt)
        if (optmas(1:14) .eq. 'MASS_MECA_DIAG') then
            massDiagIndx = 0
            call utmess('I', 'UTILITAI3_72')
        end if
    end if

! - Create input field for lumped mass
    call mecact('V', chmasd, 'MAILLA', mesh, 'POSI', &
                ncmp=1, nomcmp='POS', si=massDiagIndx)
!
    do iStore = 1, nbStore
        call jemarq()
        call jerecu('V')
        l_modal = ASTER_FALSE

! ----- Current storing index
        numeStore = listStore(iStore)
        inst = listTime(iStore)
        ASSERT(inst .ne. rundf)
        valer(1) = inst
        if (resultType .eq. 'FOURIER_ELAS') then
            call rsadpa(result, 'L', 1, 'NUME_MODE', numeStore, 0, sjv=jnmo)
            call meharm(modelZ, zi(jnmo), chharm)
        end if
        chtime = ' '
        if (lHasTime) then
            call mechti(mesh, inst, rundf, rundf, chtime)
        end if
!
        if (lResultUser) then
            call rsexch(' ', result, 'ECIN_ELEM', numeStore, ecinElemUser, iret)
            if (iret .gt. 0) then
                if (lHasTime) then
                    call rsexch(' ', result, 'VITE', numeStore, field, iret)
                    if (iret .gt. 0) goto 72
                    call dismoi('NOM_GD', field, 'CHAMP', repk=physQuanName)
                    call dismoi('TYPE_SUPERVIS', field, 'CHAMP', repk=fieldType)
                else
                    l_modal = ASTER_TRUE
                    call rsexch(' ', result, 'DEPL', numeStore, field, iret)
                    if (iret .gt. 0) goto 72
                    call dismoi('NOM_GD', field, 'CHAMP', repk=physQuanName)
                    call dismoi('TYPE_SUPERVIS', field, 'CHAMP', repk=fieldType)
                end if
            else
                call dismoi('NOM_GD', ecinElemUser, 'CHAMP', repk=physQuanName)
                call dismoi('TYPE_SUPERVIS', ecinElemUser, 'CHAMP', repk=fieldType)
            end if
            if (lHasTime) then
                xfreq = 1.d0
            else
                call rsadpa(result, 'L', 1, 'OMEGA2', numeStore, 0, sjv=lfreq)
                xfreq = zr(lfreq)
            end if
        end if

! ----- Create field for frequency
        call mecact('V', chfreq, 'MAILLA', mesh, 'OME2_R', &
                    ncmp=1, nomcmp='OMEG2', sr=xfreq)
!
        if (fieldType(1:7) .eq. 'CHAM_NO') then
            if (physQuanName(1:4) .eq. 'DEPL') then
                call vrcins(modelZ, materFieldZ, caraElemZ, inst, chvarc, codret)
            else
                call utmess('F', 'UTILITAI3_73')
            end if

        else if (fieldType(1:9) .eq. 'CHAM_ELEM') then
            if (physQuanName(1:4) .eq. 'ENER') then
                ecinElem = ecinElemUser
                goto 30
            else
                call utmess('F', 'UTILITAI3_73')
            end if
        else
            call utmess('F', 'UTILITAI3_73')
        end if
        ecinElem = '&&PEECIN.CHAM_ELEM'
        ibid = 0
        if (l_modal) then
            chvite = ' '
            chdisp = field
        else
            chvite = field
            chdisp = ' '
        end if

! ----- Compute kinetic energy
        call compEnergyKinetic(l_modal, modelZ, materCodeZ, caraElemZ, &
                               chdisp, chvite, chFreq, chgeom, &
                               chmasd, chvarc, &
                               ligrel, jvBase, ecinElem, iret)
30      continue
!
!        --- ON CALCULE L'ENERGIE TOTALE ---
        call peenca(ecinElem, nbpaep, varpep, 0, [ibid])
!
        do iFactorKeyword = 1, nbFactorKeyword
            call getvtx(option(1:9), 'TOUT', iocc=iFactorKeyword, nbval=0, nbret=nt)
            call getvem(mesh, 'MAILLE', option(1:9), 'MAILLE', iFactorKeyword, &
                        0, k8b, nm)
            call getvem(mesh, 'GROUP_MA', option(1:9), 'GROUP_MA', iFactorKeyword, &
                        0, k8b, ng)
            if (nt .ne. 0) then
                call peenca(ecinElem, nbpaep, varpep, 0, [ibid])
                valk(1) = mesh
                valk(2) = 'TOUT'
                if (lResultUser) then
                    valer(2) = varpep(1)
                    valer(3) = varpep(2)
                    call tbajli(tablOutZ, nbParaResu, tablRParaName, [numeStore], valer, &
                                [c16b], valk, 0)
                else
                    call tbajli(tablOutZ, nbParaField, tablFParaName, [numeStore], varpep, &
                                [c16b], valk, 0)
                end if
            end if
            if (ng .ne. 0) then
                nbgrma = -ng
                call wkvect('&&PEECIN_GROUPM', 'V V K24', nbgrma, jgr)
                call getvem(mesh, 'GROUP_MA', option(1:9), 'GROUP_MA', iFactorKeyword, &
                            nbgrma, zk24(jgr), ng)
                valk2(2) = 'GROUP_MA'
                do ig = 1, nbgrma
                    nomgrm = zk24(jgr+ig-1)
                    call jeexin(jexnom(mesh//'.GROUPEMA', nomgrm), iret)
                    if (iret .eq. 0) then
                        call utmess('A', 'UTILITAI3_46', sk=nomgrm)
                        goto 40
                    end if
                    call jelira(jexnom(mesh//'.GROUPEMA', nomgrm), 'LONUTI', nbma)
                    if (nbma .eq. 0) then
                        call utmess('A', 'UTILITAI3_47', sk=nomgrm)
                        goto 40
                    end if
                    call jeveuo(jexnom(mesh//'.GROUPEMA', nomgrm), 'L', jad)
                    call peenca(ecinElem, nbpaep, varpep, nbma, zi(jad))
                    valk2(1) = nomgrm
                    if (lResultUser) then
                        valer(2) = varpep(1)
                        valer(3) = varpep(2)
                        call tbajli(tablOutZ, nbParaResu, tablRParaName, [numeStore], valer, &
                                    [c16b], valk2, 0)
                    else
                        call tbajli(tablOutZ, nbParaField, tablFParaName, [numeStore], varpep, &
                                    [c16b], valk2, 0)
                    end if
40                  continue
                end do
                call jedetr('&&PEECIN_GROUPM')
            end if
            if (nm .ne. 0) then
                nbma = -nm
                call wkvect('&&PEECIN_MAILLE', 'V V K8', nbma, jma)
                call getvem(mesh, 'MAILLE', option(1:9), 'MAILLE', iFactorKeyword, &
                            nbma, zk8(jma), nm)
                valk(2) = 'MAILLE'
                call jelira(mesh//'.TYPMAIL', 'LONMAX', nbMaiT)
                do im = 1, nbma
                    nommai = zk8(jma+im-1)
                    nume = char8_to_int(nommai)
                    if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                        call utmess('A', 'UTILITAI3_49', sk=nommai)
                        goto 50
                    end if
                    call peenca(ecinElem, nbpaep, varpep, 1, [nume])
                    valk(1) = nommai
                    if (lResultUser) then
                        valer(2) = varpep(1)
                        valer(3) = varpep(2)
                        call tbajli(tablOutZ, nbParaResu, tablRParaName, [numeStore], valer, &
                                    [c16b], valk, 0)
                    else
                        call tbajli(tablOutZ, nbParaField, tablFParaName, [numeStore], varpep, &
                                    [c16b], valk, 0)
                    end if
50                  continue
                end do
                call jedetr('&&PEECIN_MAILLE')
            end if
        end do
        call jedetr('&&PEECIN.PAR')
72      continue
        call jedema()
    end do
!
80  continue
    call jedetr(listStoreJv)
    call jedetr(listTimeJv)
!
90  continue
    call jedema()
end subroutine
