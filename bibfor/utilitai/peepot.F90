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
subroutine peepot(tablOutZ, &
                  modelZ, materFieldZ, materCodeZ, caraElemZ, &
                  numeHarm, nbFactorKeyword)
!
    implicit none
!
#include "asterc/asmpi_comm.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/celver.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/chpve2.h"
#include "asterfort/compEnergyPotential.h"
#include "asterfort/digdel.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim3.h"
#include "asterfort/gettco.h"
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
#include "asterfort/mecham.h"
#include "asterfort/mechti.h"
#include "asterfort/meharm.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/peenca2.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/scalai.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/vecint.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: tablOutZ, modelZ, materFieldZ, materCodeZ, caraElemZ
    integer(kind=8), intent(in) :: numeHarm, nbFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "ENER_POT"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'ENER_POT'
    integer(kind=8) :: iret, np, nc, jad, nbStore, iStore, numeStore, jvPara
    integer(kind=8) :: ire1, ire2, nt, nm, ng, nbgrma, ig, jgr, nbma, nume, im, nbret
    integer(kind=8) :: iFactorKeyword, jma, icheml, ier, nbMaiT, jnmo, ibid
    integer(kind=8), parameter :: nbpaep = 2
    real(kind=8) :: varpep(nbpaep)
    real(kind=8) :: prec, inst, valer(3), rundf
    character(len=1), parameter :: jvBase = "V"
    character(len=2) :: codret
    character(len=8) :: k8b, mesh, result, crit, nommai, valk(2)
    character(len=8) :: physQuanName
    character(len=16) :: resultType, optio2
    character(len=19) :: epotElem, ligrel, ligrel2
    character(len=19) :: field, epotElemUser
    character(len=24) :: chtime, fieldType, chgeom, chtemp, chharm, chdisp
    character(len=24) :: compor, nomgrm, valk2(2)
    aster_logical :: l_temp
    complex(kind=8) :: c16b
    character(len=19), parameter :: chvarc = '&&PEECIN.VARC', chvref = '&&PEEPOT.VARC_REF'
    character(len=19), parameter :: listStoreJv = '&&PEECIN.NUME_ORDRE'
    integer(kind=8), pointer :: listStore(:) => null()
    character(len=19), parameter :: listTimeJv = '&&PEECIN.INSTANT'
    real(kind=8), pointer :: listTime(:) => null()
    aster_logical :: lFieldUser, lResultUser, lHasTime
    mpi_int :: mpicow, mrang, mnbproc, mpicou
    integer(kind=8) :: rang, nbproc, k, ntsum, nmsum, nmmax, ngsum, ngmax
    integer(kind=8) :: decalig, decalim, jmntmg, jmigk, jmigi, jmim, niv, ifm, nbgr
    integer(kind=8) :: numpas, numloc
    integer(kind=8) :: longt, icoef, mode, nel, idecgr, j, nbmasum, jnp, ind
    real(kind=8) :: ztot
    character(len=4) :: docu
    character(len=8) :: k8X, scal
    character(len=24) :: k24X
    character(len=24), pointer :: celk(:) => null()
    integer(kind=8), pointer :: celd(:) => null()
    real(kind=8), pointer :: celv(:) => null()
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
    call infniv(ifm, niv)

! - Initializations
    c16b = (0.d0, 0.d0)
    rundf = r8vide()
    lHasTime = ASTER_FALSE
    inst = 0.d0
    ntsum = 0
    nmsum = 0
    ngsum = 0
    nmmax = 0
    ngmax = 0
    k8X = 'XXXXXXXX'
    k24X = 'XXXXXXXXXXXXXXXXXXXXXXXX'
    chtemp = ' '
    chdisp = ' '
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
                 resultType(1:9) .eq. 'MULT_ELAS' .or. resultType(1:9) .eq. 'EVOL_NOLI' .or. &
                 resultType(1:10) .eq. 'DYNA_TRANS') then
            tablRParaName(2) = 'INST'
        else
            ASSERT(ASTER_FALSE)
        end if
    end if

! - Prepare input fields
    call mecham(option, modelZ, numeHarm, &
                chgeom, chharm, iret)
    if (iret .ne. 0) goto 90
    mesh = chgeom(1:8)
!
    call exlim3(option, 'V', modelZ, ligrel)
!
    if (lFieldUser) then
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

!        --- ON RECUPERE LES INSTANTS ---
        call wkvect(listTimeJv, 'V V R', nbStore, vr=listTime)
        call jenonu(jexnom(result//'           .NOVA', 'INST'), iret)
        if (iret .ne. 0) then
            lHasTime = ASTER_TRUE
            do iStore = 1, nbStore
                numeStore = listStore(iStore)
                call rsadpa(result, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
                listTime(iStore) = zr(jvPara)
            end do
        else
            call jenonu(jexnom(result//'           .NOVA', 'FREQ'), iret)
            if (iret .ne. 0) then
                do iStore = 1, nbStore
                    numeStore = listStore(iStore)
                    call rsadpa(result, 'L', 1, 'FREQ', numeStore, 0, sjv=jvPara)
                    listTime(iStore) = zr(jvPara)
                end do
            end if
        end if
        call tbcrsd(tablOutZ, 'G')
        call tbajpa(tablOutZ, nbParaResu, tablRParaName, tablRParaType)
    end if
!-----------------------------------------------------------------------------
! MUTUALISATION POUR APPELS GETVTX
! AFIN DE NE PAS LE REFAIRE POUR CHAQUE PAS DE TEMPS
!-----------------------------------------------------------------------------
    call wkvect('&&PEEPOT_jmntmg', 'V V I', 3*nbFactorKeyword, jmntmg)

    do iFactorKeyword = 1, nbFactorKeyword
        call getvtx(option(1:9), 'TOUT', iocc=iFactorKeyword, nbval=0, nbret=nt)
        call getvem(mesh, 'MAILLE', option(1:9), 'MAILLE', iFactorKeyword, 0, k8b, nm)
        call getvem(mesh, 'GROUP_MA', option(1:9), 'GROUP_MA', iFactorKeyword, 0, k8b, ng)
        zi(jmntmg+3*(iFactorKeyword-1)) = nt
        ntsum = ntsum+abs(nt)
        zi(jmntmg+3*(iFactorKeyword-1)+1) = nm
        nmmax = max(nmmax, abs(nm))
        nmsum = nmsum+abs(nm)
        zi(jmntmg+3*(iFactorKeyword-1)+2) = ng
        ngmax = max(ngmax, abs(ng))
        ngsum = ngsum+abs(ng)
    end do
!
    ASSERT((ntsum .ge. 0) .and. (ngsum .ge. 0) .and. (nmsum .ge. 0))
    ASSERT((ngmax .ge. 0) .and. (ngmax .le. ngsum))
    ASSERT((nmmax .ge. 0) .and. (nmmax .le. nmsum))
    if (ngsum .gt. 0) then
        call wkvect('&&PEEPOT_jmigk', 'V V K24', ngsum, jmigk)
        call wkvect('&&PEEPOT_jmigi', 'V V I', ngsum, jmigi)
    end if
    if (nmsum .gt. 0) then
        call wkvect('&&PEEPOT_jmim', 'V V K8', nmsum, jmim)
    end if
    decalig = 0
    decalim = 0
    nbmasum = 0
    do iFactorKeyword = 1, nbFactorKeyword
        ng = zi(jmntmg+3*(iFactorKeyword-1)+2)
        if (ng .ne. 0) then
            nbgrma = -ng
            call wkvect('&&PEEPOT_GROUPM', 'V V K24', nbgrma, jgr)
            call getvem(mesh, 'GROUP_MA', option(1:9), 'GROUP_MA', iFactorKeyword, &
                        nbgrma, zk24(jgr), ng)
            do ig = 1, nbgrma
                nomgrm = zk24(jgr+ig-1)
                zk24(jmigk-1+ig+decalig) = nomgrm
                call jeexin(jexnom(mesh//'.GROUPEMA', nomgrm), iret)
                if (iret .eq. 0) then
                    call utmess('A', 'UTILITAI3_46', sk=nomgrm)
                    zk24(jmigk-1+ig+decalig) = k24X
                    zi(jmigi-1+ig+decalig) = -999
                    goto 140
                end if
                call jelira(jexnom(mesh//'.GROUPEMA', nomgrm), 'LONUTI', nbma)
                if (nbma .eq. 0) then
                    call utmess('A', 'UTILITAI3_47', sk=nomgrm)
                    zk24(jmigk-1+ig+decalig) = k24X
                    zi(jmigi-1+ig+decalig) = -999
                    goto 140
                else
                    zi(jmigi-1+ig+decalig) = nbma
                    nbmasum = nbmasum+nbma
                end if
140             continue
            end do
            call jedetr('&&PEEPOT_GROUPM')
            decalig = decalig+nbgrma
! fin if sur nm (groupe de mailles)
        end if
        nm = zi(jmntmg+3*(iFactorKeyword-1)+1)
        if (nm .ne. 0) then
            nbma = -nm
            call wkvect('&&PEEPOT_MAILLE', 'V V K8', nbma, jma)
            call getvem(mesh, 'MAILLE', option(1:9), 'MAILLE', iFactorKeyword, &
                        nbma, zk8(jma), nm)
            nbmasum = nbmasum+nbma
            call jelira(mesh//'.TYPMAIL', 'LONMAX', nbMaiT)
            do im = 1, nbma
                nommai = zk8(jma+im-1)
                nume = char8_to_int(zk8(jma+im-1))
                if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                    call utmess('A', 'UTILITAI3_49', sk=nommai)
                    zk8(jmim-1+im+decalim) = k8X
                    goto 150
                else
                    zk8(jmim-1+im+decalim) = nommai
                end if
150             continue
            end do
            call jedetr('&&PEEPOT_MAILLE')
            decalim = decalim+nbma
        end if
    end do
!
!-----------------------------------------------------------------------------
! PREPARATION DE LA DISTRIBUTION DE TACHES MPI VIA
! FILTRE &PEECA2_vldist (POUR CELUI EN ESPACE-CONNECTIVITE INVERSE DE PEENCA2)
!-----------------------------------------------------------------------------
! Recuperation des donnees MPI pour le //isme en espace de peenca2 (actif par defaut)
    call asmpi_comm('GET_WORLD', mpicow)
    call asmpi_comm('GET', mpicou)
    ASSERT(mpicow .eq. mpicou)
    call asmpi_info(mpicow, mrang, mnbproc)
    rang = to_aster_int(mrang)
    ASSERT(rang .ge. 0)
    nbproc = to_aster_int(mnbproc)
    ASSERT(nbproc .ge. 1)
!
!-----------------------------------------------------------------------------
! MUTUALISATION POUR APPEL PEENCA (step 1): OBJET POUR STOCKER LA CONNECTIVITE INVERSE
! POUR LA KIEME MAILLE DU GROUP_MA: (cf. PEENCA)
! &&PEEPOT_peenca(2*(k-1)+1)=numero du GREL
! &&PEEPOT_peenca(2*(k-1)+2)=indice de l element dans ce GREL
! AFIN DE NE PAS LE REFAIRE POUR CHAQUE MAILLE DE CALCUL DES GROUP_MA OU DES LISTE
! DE MAILLES ET POUR CHAQUE PAS DE TEMPS
!-----------------------------------------------------------------------------
    if (nbmasum .gt. 0) then
        call wkvect('&&PEEPOT_peenca', 'V V I', 2*nbmasum, jnp)
        call vecint(2*nbmasum, 0, zi(jnp))
    end if
    numpas = 0
!
!-----------------------------------------------------------------------------
! BOUCLE PRINCIPALE: PAS DE TEMPS OU MODES OU...
!-----------------------------------------------------------------------------
!
    do iStore = 1, nbStore
        numpas = numpas+1
        numloc = iStore-(numpas-1)*nbproc
        call jemarq()
        call jerecu('V')
        icheml = 0
        numeStore = listStore(iStore)
        inst = listTime(iStore)
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
            call rsexch(' ', result, 'EPOT_ELEM', numeStore, epotElemUser, iret)
            if (iret .gt. 0) then
                call rsexch(' ', result, 'DEPL', numeStore, field, ire1)
                if (ire1 .gt. 0) then
                    call rsexch(' ', result, 'TEMP', numeStore, field, ire2)
                    if (ire2 .gt. 0) goto 72
                    call dismoi('NOM_GD', field, 'CHAMP', repk=physQuanName)
                    call dismoi('TYPE_SUPERVIS', field, 'CHAMP', repk=fieldType)
                else
                    call dismoi('NOM_GD', field, 'CHAMP', repk=physQuanName)
                    call dismoi('TYPE_SUPERVIS', field, 'CHAMP', repk=fieldType)
                end if
            else
                call dismoi('TYPE_SUPERVIS', epotElemUser, 'CHAMP', repk=fieldType)
                call dismoi('NOM_GD', epotElemUser, 'CHAMP', repk=physQuanName)
            end if
        end if
!
        if (fieldType(1:7) .eq. 'CHAM_NO') then
            call vrcins(modelZ, materFieldZ, caraElemZ, inst, chvarc, codret)
            call vrcref(modelZ(1:8), materFieldZ(1:8), caraElemZ(1:8), chvref(1:19))
            if (physQuanName(1:4) .eq. 'DEPL') then
                optio2 = 'EPOT_ELEM'
                l_temp = ASTER_FALSE
            else if (physQuanName(1:4) .eq. 'TEMP') then
                optio2 = 'ETHE_ELEM'
                l_temp = ASTER_TRUE
            else
                call utmess('F', 'UTILITAI3_73')
            end if

        else if (fieldType(1:9) .eq. 'CHAM_ELEM') then
            if (physQuanName(1:4) .eq. 'ENER') then
                epotElem = epotElemUser
                goto 30
            else
                call utmess('F', 'UTILITAI3_73')
            end if
        else
            call utmess('F', 'UTILITAI3_73')
        end if
        icheml = 1
        epotElem = '&&PEEPOT.CHAM_ELEM'
        compor = materFieldZ(1:8)//'.COMPOR'
        ibid = 0
        if (l_temp) then
            chtemp = field
            chdisp = ' '
        else
            chdisp = field
            chtemp = ' '
        end if

        call compEnergyPotential(optio2, &
                                 modelZ, materCodeZ, caraElemZ, compor, &
                                 chdisp, chharm, chgeom, &
                                 chtime, chvarc, chvref, &
                                 l_temp, chtemp, &
                                 ligrel, jvBase, epotElem, &
                                 iret)
30      continue
!
!-----------------------------------------------------------------------------
! MUTUALISATION POUR APPEL PEENCA (step 2) AFIN DE NE PAS LE REFAIRE POUR
! CHAQUE MAILLE DE CALCUL DES GROUP_MA OU DES LISTE DE MAILLES.
!
! VERIFICATIONS AU PREMIER PAS DE TEMPS, CALCUL #GREL (NBGR), NOM DU LIGREL (LIGREL2),
! TYPE DE CHAMPS (SCAL), ENERGIE TOTALE (ZTOT)
!-----------------------------------------------------------------------------
        if (iStore .eq. 1) then
! on fait ces verifications qu'au premier pas de temps, cela suffit ici
            call celver(epotElem, 'NBVARI_CST', 'STOP', ibid)
            call celver(epotElem, 'NBSPT_1', 'STOP', ibid)
            call jelira(epotElem//'.CELD', 'DOCU', cval=docu)
            if (docu .ne. 'CHML') then
                call utmess('F', 'CALCULEL3_52')
            end if
        end if
        call jeveuo(epotElem//'.CELK', 'L', vk24=celk)
        call jeveuo(epotElem//'.CELD', 'L', vi=celd)
        call jeveuo(epotElem//'.CELV', 'L', vr=celv)
        ligrel2 = celk(1) (1:19)
        nbgr = nbgrel(ligrel2)
        if (iStore .eq. 1) then
            scal = scalai(celd(1))
            if (scal(1:1) .ne. 'R') then
                call utmess('F', 'CALCULEL3_74', sk=scal)
            end if
        end if
        ztot = 0.d0
        do j = 1, nbgr
            mode = celd(celd(4+j)+2)
            if (mode .eq. 0) goto 34
            longt = digdel(mode)
            icoef = max(1, celd(4))
            longt = longt*icoef
            nel = nbelem(ligrel2, j)
            idecgr = celd(celd(4+j)+8)
            do k = 1, nel
                ztot = ztot+celv(idecgr+(k-1)*longt)
            end do
34          continue
        end do

!
! CALCUL ENERGIE TOTALE DEJA DISPONIBLE
        varpep(1) = ztot
        varpep(2) = 100.d0
        decalig = 0
        decalim = 0
        if (numpas .eq. 1) then
            ind = -1
        else
            ind = 1
        end if
!
!-----------------------------------------------------------------------------
! BOUCLE SECONDAIRE: LISTE DE GROUP_MA OU DE MAILLES
!-----------------------------------------------------------------------------
!
        do iFactorKeyword = 1, nbFactorKeyword
! Resultats getvtx deja lus une fois pour toute
            nt = zi(jmntmg+3*(iFactorKeyword-1))
            nm = zi(jmntmg+3*(iFactorKeyword-1)+1)
            ng = zi(jmntmg+3*(iFactorKeyword-1)+2)
!
! Calcul sur 'TOUT'
            if (nt .ne. 0) then
                varpep(1) = ztot
                varpep(2) = 100.d0
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
!
! Calcul sur GROUP_MA
            if (ng .ne. 0) then
                nbgrma = -ng
                valk2(2) = 'GROUP_MA'
                do ig = 1, nbgrma
                    nomgrm = zk24(jmigk-1+ig+decalig)
                    if (nomgrm(1:24) .ne. k24X) then
                        nbma = zi(jmigi-1+ig+decalig)
                        ASSERT(nbma .ne. -999)
                        call jeveuo(jexnom(mesh//'.GROUPEMA', nomgrm), 'L', jad)
                        call peenca2(epotElem, nbpaep, varpep, nbma, zi(jad), &
                                     ligrel2, nbgr, ztot, ind, nbproc, &
                                     rang)
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
                    end if
                end do
                decalig = decalig+nbgrma
            end if
!
! Calcul sur liste de MA
            if (nm .ne. 0) then
                nbma = -nm
                valk(2) = 'MAILLE'
                do im = 1, nbma
                    nommai = zk8(jmim-1+im+decalim)
                    if (nommai .ne. k8X) then
                        nume = char8_to_int(nommai)
                        call peenca2(epotElem, nbpaep, varpep, 1, [nume], &
                                     ligrel2, nbgr, ztot, ind, nbproc, &
                                     rang)
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
                    end if
                end do
                decalim = decalim+nbma
            end if
        end do
        call jedetr('&&PEEPOT.PAR')
        if (icheml .ne. 0) call jedetr(epotElem)
72      continue
        call jedema()
    end do
!
! Nettoyage des objets de mutualisations et des buffers de com mpi
    call jedetr('&&PEEPOT_jmntmg')
    if (ngsum .gt. 0) then
        call jedetr('&&PEEPOT_jmigk')
        call jedetr('&&PEEPOT_jmigi')
    end if
    if (nmsum .gt. 0) then
        call jedetr('&&PEEPOT_jmim')
    end if
    if (nbmasum .gt. 0) then
        call jedetr('&&PEEPOT_peenca')
    end if
!
80  continue
    call jedetr(listStoreJv)
    call jedetr(listTimeJv)
!
90  continue
    call jedema()
end subroutine
