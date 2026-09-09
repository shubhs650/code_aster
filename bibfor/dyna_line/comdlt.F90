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
subroutine comdlt()
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cochre.h"
#include "asterfort/compStrx.h"
#include "asterfort/dismoi.h"
#include "asterfort/dladap.h"
#include "asterfort/dldiff.h"
#include "asterfort/dlnewi.h"
#include "asterfort/dltali.h"
#include "asterfort/dltlec.h"
#include "asterfort/dvcrob.h"
#include "asterfort/dyGetKineLoad.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lobs.h"
#include "asterfort/mecham.h"
#include "asterfort/mechti.h"
#include "asterfort/nmch1p.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmobse.h"
#include "asterfort/nmobsw.h"
#include "asterfort/nonlinDSEnergyClean.h"
#include "asterfort/nonlinDSInOutClean.h"
#include "asterfort/nonlinDSInOutInit.h"
#include "asterfort/nonlinDSInOutRead.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! DYNA_VIBRA
!
! TYPE_CALCUL = 'TRAN' + BASE_CALCUL = 'PHYS'
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: numeHarm = 0
    integer(kind=8) :: nbVectAsse, nbLoad
    integer(kind=8) :: jvMatr(3), nume, niv, ifm, jvLoadWave, jvPara, numrep
    integer(kind=8) :: jvVectFunc, jvVectAsse, nbWave, ifexte, ifamor, ifliai
    integer(kind=8) :: neq, idepl0, ivite0, iacce0, iwk, iStore
    integer(kind=8) :: iinteg, iret, nbpas, nbpas_min, nbpas_max
    integer(kind=8) :: nbStore, jchar, jinst, pasar, nbar
    integer(kind=8) :: lresu, lcrre, iresu, nbexre, l, ncomu
    integer(kind=8) :: nbchre, iocc, nfon, nbexcl, i, counter, lsize
    real(kind=8) :: t0, time, rundf, coefMultR, tinit
    real(kind=8) :: tfin, dt, dtmin, dtmax, cdivi
    real(kind=8) :: nbpas_max_r, epsi
    character(len=1) :: jvBase
    character(len=2) :: codret
    character(len=8) :: existPoux, masse, rigid, amort, result
    character(len=8) :: existStr2, nomfon, loadPres
    character(len=9) :: nomsym(6)
    character(len=19) :: solveu, listLoad, ligrel, linst
    character(len=12) :: allschemes(4), schema, schtyp
    character(len=24) :: model, caraElem, loadNameJv, loadFuncJv, materCode, materField
    character(len=24) :: numedd, disp
    character(len=24) :: loadInfoJv, criter
    character(len=24) :: chgeom, chharm, chtime
    character(len=24) :: chvarc, chvref, strx, compor, kineLoad
    character(len=19) :: force0, force1
    character(len=19) :: sd_obsv
    type(NL_DS_InOut) :: ds_inout
    character(len=46) :: champs
    type(NL_DS_Energy) :: ds_energy
    aster_logical :: lamort, lcrea, lprem, lPoux
    integer(kind=8), pointer :: listStore(:) => null()
    character(len=8), pointer :: chexc(:) => null()
    data model/'                        '/
    data allschemes/'NEWMARK', 'WILSON', 'DIFF_CENTRE', 'ADAPT_ORDRE2'/
    character(len=8) :: mesh
    aster_logical :: l_obsv
    integer(kind=8), parameter :: zvalin = 28
    character(len=19) :: valinc(zvalin)
    character(len=19) :: depmoi, vitmoi, accmoi
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call titre()
    call infmaj()
    call infniv(ifm, niv)
!
! - Initializations
!
    rundf = r8vide()
    epsi = 1000.d0*r8prem()
    lprem = .true.
    lamort = .true.
    amort = ' '
    criter = '&&RESGRA_GCPC'
    coefMultR = 0.d0
    nfon = 0
    ncomu = 0
    loadPres = ' '
    chtime = ' '
    jvBase = 'G'
    nbpas = 0
    nbpas_min = 0
    nbpas_max = 0
    dt = 0.d0
    dtmin = 0.d0
    dtmax = 0.d0
    chvarc = '&&COMDLT.VARC'
    chvref = '&&COMDLT.VREF'
!
! - Get parameters
    call dltlec(result, model, numedd, materField, materCode, &
                caraElem, jvMatr, masse, rigid, amort, &
                lamort, nbLoad, nbVectAsse, listLoad, loadNameJv, &
                loadInfoJv, loadFuncJv, jvVectAsse, jvVectFunc, nbWave, &
                jvLoadWave, solveu, iinteg, t0, nume, &
                numrep, ds_inout)
    if (nbWave .gt. 0 .and. materField .eq. ' ') then
        call utmess('F', 'DYNALINE1_4')
    end if
!
! - Get kinematic loads
    call dyGetKineLoad(rigid, masse, amort, lamort, listLoad, &
                       kineLoad, iinteg)
!
    neq = zi(jvMatr(1)+2)

    ! WRITE (6, *) "CHARGEMENTS"
    ! WRITE (6, *) "==========="
    ! WRITE (6, *) "nbLoad: ", nbLoad, listLoad
    ! WRITE (6, *) "nbVectAsse: ", nbVectAsse, jvVectAsse, jvVectFunc
    ! WRITE (6, *) "nbWave: ", nbWave, jvLoadWave

!
!====
! 3. CREATION DES VECTEURS DE TRAVAIL SUR BASE VOLATILE
!====
!   NOUVEAUX OBJETS DEPL VITE ACCE DE DYNA_LINE_TRAN
! ----------------------------------------
!     ancienne routine : nmcrch
! --- CREATION DES CHAMPS DE BASE - ETAT EN T-
!
    call nmch1p(valinc)
    call nmchex(valinc, 'VALINC', 'DEPMOI', depmoi)
    call nmchex(valinc, 'VALINC', 'VITMOI', vitmoi)
    call nmchex(valinc, 'VALINC', 'ACCMOI', accmoi)
!
    call vtcreb(depmoi, 'V', 'R', nume_ddlz=numedd)
    call vtcreb(vitmoi, 'V', 'R', nume_ddlz=numedd)
    call vtcreb(accmoi, 'V', 'R', nume_ddlz=numedd)
!
!---lecture d adresse depl/vite/acce
!
    call jeveuo(depmoi//'.VALE', 'L', idepl0)
    call jeveuo(vitmoi//'.VALE', 'L', ivite0)
    call jeveuo(accmoi//'.VALE', 'L', iacce0)
!
    call wkvect('&&COMDLT.FEXTE', 'V V R', 2*neq, ifexte)
    call wkvect('&&COMDLT.FAMOR', 'V V R', 2*neq, ifamor)
    call wkvect('&&COMDLT.FLIAI', 'V V R', 2*neq, ifliai)
    call wkvect('&&COMDLT.TRAV', 'V V R', neq, iwk)
!
    call getfac('EXCIT_RESU', nbexre)
    if (nbexre .ne. 0) then
        call wkvect('&&COMDLT.COEF_RRE', 'V V R  ', nbexre, lcrre)
        call wkvect('&&COMDLT.LISTRESU', 'V V K8 ', nbexre, lresu)
        do iresu = 1, nbexre
            call getvid('EXCIT_RESU', 'RESULTAT', iocc=iresu, scal=zk8(lresu+iresu-1), nbret=l)
            call getvr8('EXCIT_RESU', 'COEF_MULT', iocc=iresu, scal=zr(lcrre+iresu-1), nbret=l)
        end do
    end if
!
!===
! 4. IMPRESSIONS RECAPITULATIVES POUR L'UTILISATEUR
!===
!
    call utmess('I', 'DYNAMIQUE_55', nk=3, &
                valk=['D Y N A _ V I B R A', 'TRANsitoire        ', 'PHYSique           '])
    call utmess('I', 'DYNAMIQUE_82', sk=model, si=neq)
    call utmess('I', 'DYNAMIQUE_60')
    call utmess('I', 'DYNAMIQUE_61', nk=2, valk=[masse, rigid])
    if (lamort) then
        call utmess('I', 'DYNAMIQUE_62', sk=amort)
    else
        call utmess('I', 'DYNAMIQUE_64')
    end if
!
    schema = allschemes(iinteg)
    schtyp = 'explicite'
    if ((iinteg .eq. 1) .or. (iinteg .eq. 2)) schtyp = 'implicite'
    call getvid('INCREMENT', 'LIST_INST', iocc=1, scal=linst, nbret=iret)
    if (iret .eq. 0) then
        call getvr8('INCREMENT', 'INST_INIT', iocc=1, scal=tinit, nbret=iret)
        if (iret .eq. 0) tinit = 0.d0
        call getvr8('INCREMENT', 'INST_FIN', iocc=1, scal=tfin)
        call getvr8('INCREMENT', 'PAS', iocc=1, scal=dt, nbret=iret)
    else
        call jeveuo(linst//'.VALE', 'L', jinst)
        call jelira(linst//'.VALE', 'LONMAX', nbpas)
        tinit = zr(jinst)
        tfin = zr(jinst+nbpas-1)
        ASSERT(nbpas .gt. 1)
        dt = (tfin-tinit)/real(nbpas-1)

        call getvis('INCREMENT', 'NUME_FIN', iocc=1, scal=iStore, nbret=iret)
        if (iret .ne. 0) then
            if (iStore .ge. nbpas) goto 99
            tfin = zr(jinst+iStore)
        else
            call getvr8('INCREMENT', 'INST_FIN', iocc=1, scal=tfin, nbret=iret)
        end if
    end if
99  continue
!

! -- Check
!
    if (dt .le. 0.0) then
        call utmess('F', 'DYNAMIQUE_42', sr=dt)
    end if
    if (tfin .le. tinit) then
        call utmess('F', 'DYNAMIQUE_43', nr=2, valr=[tfin, tinit])
    end if
!
    if (iinteg .ne. 4) then
        nbpas = nint((tfin-tinit)/dt)
        call utmess('I', 'DYNAMIQUE_70', nk=2, valk=[schema, schtyp], nr=1, &
                    valr=[dt], ni=1, vali=[nbpas])
    else
        call getvr8('SCHEMA_TEMPS', 'PAS_MINI', iocc=1, scal=dtmin, nbret=iret)
        if (iret .eq. 0) dtmin = dt*1.d-6
        call getvr8('SCHEMA_TEMPS', 'PAS_MAXI', iocc=1, scal=dt, nbret=iret)
        if (iret .eq. 0) dtmax = dt*1.d6
        call utmess('I', 'DYNAMIQUE_66', nk=1, valk=[schema], nr=3, &
                    valr=[dt, dtmin, dtmax])
        call getvr8('SCHEMA_TEMPS', 'COEF_DIVI_PAS', iocc=1, scal=cdivi)
        call utmess('I', 'DYNAMIQUE_68', nr=1, valr=[cdivi])
!
        nbpas_min = nint((tfin-tinit)/dtmax)
        nbpas_max = 1000000000
        if (dtmin .gt. epsi) then
            nbpas_max_r = min(1.d0*nbpas_max, (tfin-tinit-epsi)/dtmin)
            nbpas_max = int(nbpas_max_r)+1
        end if
        call utmess('I', 'DYNAMIQUE_69', ni=2, vali=[nbpas_min, nbpas_max])
    end if
!
!
!===
! 5. INITIALISATION DE L'ALGORITHME
!===
!
    force0 = '&&COMDLT.FORCE0'
    force1 = '&&COMDLT.FORCE1'
    call dltali(neq, result, jvMatr, masse, rigid, &
                zi(jvVectAsse), zk24(jvVectFunc), nbLoad, nbVectAsse, lcrea, &
                lprem, lamort, t0, materField, materCode, &
                caraElem, loadNameJv, loadInfoJv, loadFuncJv, model, &
                numedd, nume, solveu, criter, zr(idepl0), &
                zr(ivite0), zr(iacce0), zr(ifexte+neq), zr(ifamor+neq), zr(ifliai+neq), &
                zr(iwk), force0, force1, ds_energy, kineLoad)
!
    call utmess('I', 'DYNAMIQUE_80', nr=2, valr=[tinit, tfin])
!
    call getvis('ARCHIVAGE', 'PAS_ARCH', iocc=1, scal=pasar, nbret=iret)
    if (iret .ne. 0) then
        call utmess('I', 'DYNAMIQUE_85', si=pasar)
    else
        call getvid('ARCHIVAGE', 'LIST_INST', iocc=1, scal=linst, nbret=iret)
        if (iret .ne. 0) then
            call jelira(linst//'.VALE', 'LONMAX', nbar)
        else
            call getvr8('ARCHIVAGE', 'INST', iocc=1, nbval=0, nbret=iret)
            nbar = -iret
        end if
        call utmess('I', 'DYNAMIQUE_86', si=nbar)
    end if
!
    call getvtx('ARCHIVAGE', 'CHAM_EXCLU', iocc=1, nbval=0, nbret=iret)
!
    nomsym = ' '
    nomsym(1) = 'DEPL'
    nomsym(2) = 'VITE'
    nomsym(3) = 'ACCE'
    if (ds_energy%l_comp) then
        nomsym(4) = 'FORC_EXTE'
        nomsym(5) = 'FORC_AMOR'
        nomsym(6) = 'FORC_LIAI'
    end if
!
    if (iret .ne. 0) then
        nbexcl = -iret
        AS_ALLOCATE(vk8=chexc, size=nbexcl)
        call getvtx('ARCHIVAGE', 'CHAM_EXCLU', iocc=1, nbval=nbexcl, vect=chexc)
        do i = 1, nbexcl
            if (chexc(i) (1:4) .eq. 'DEPL') nomsym(1) = ' '
            if (chexc(i) (1:4) .eq. 'VITE') nomsym(2) = ' '
            if (chexc(i) (1:4) .eq. 'ACCE') nomsym(3) = ' '
        end do
        AS_DEALLOCATE(vk8=chexc)
    end if
!
    champs = ' '
    counter = 0
    do i = 1, 6
        if (nomsym(i) (1:4) .ne. '    ') then
            lsize = 9
            if (nomsym(i) (5:5) .eq. ' ') lsize = 4
            champs(counter+1:counter+lsize+1) = nomsym(i) (1:lsize)//' '
            counter = counter+lsize+1
        end if
    end do
    call utmess('I', 'DYNAMIQUE_96', sk=champs)
!
! -  LECTURE DES DONNEES ET INITIALISATION POUR ds_inout
!
    call nonlinDSInOutRead('VIBR', result, ds_inout)
    call nonlinDSInOutInit('VIBR', ds_inout)
!
!--- Create observation datastructure
! -  routine in stat_non_line : nmcrob => dvcrob
!
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dvcrob(mesh, model, ds_inout, materField, sd_obsv)
!
! - Make initial observation
!
    l_obsv = ASTER_FALSE
    call lobs(sd_obsv, nume, t0, l_obsv)
    if (l_obsv) then
        call nmobse(mesh, sd_obsv, t0)
        if (nume .eq. 0) then
            call nmobsw(sd_obsv, ds_inout)
        end if
    end if
!
!====
! 6. INTEGRATION SELON LE TYPE SPECIFIE
!====
!
    if (iinteg .eq. 1) then
!
        call dlnewi(result, force0, force1, lcrea, lamort, &
                    iinteg, neq, jvMatr, masse, rigid, &
                    amort, zr(idepl0), zr(ivite0), zr(iacce0), zr(ifexte), &
                    zr(ifamor), zr(ifliai), t0, nbLoad, nbVectAsse, &
                    zi(jvVectAsse), zk24(jvVectFunc), model, materField, materCode, &
                    caraElem, loadNameJv, loadInfoJv, loadFuncJv, numedd, &
                    nume, solveu, criter, zk8(jvLoadWave), nbWave, &
                    numrep, ds_energy, sd_obsv, mesh, kineLoad)
!
    else if (iinteg .eq. 2) then
!
        call dlnewi(result, force0, force1, lcrea, lamort, &
                    iinteg, neq, jvMatr, masse, rigid, &
                    amort, zr(idepl0), zr(ivite0), zr(iacce0), zr(ifexte), &
                    zr(ifamor), zr(ifliai), t0, nbLoad, nbVectAsse, &
                    zi(jvVectAsse), zk24(jvVectFunc), model, materField, materCode, &
                    caraElem, loadNameJv, loadInfoJv, loadFuncJv, numedd, &
                    nume, solveu, criter, zk8(jvLoadWave), nbWave, &
                    numrep, ds_energy, sd_obsv, mesh, kineLoad)
!
    else if (iinteg .eq. 3) then
!
        call dldiff(result, force1, lcrea, lamort, neq, &
                    jvMatr, masse, rigid, amort, zr(idepl0), &
                    zr(ivite0), zr(iacce0), zr(ifexte), zr(ifamor), zr(ifliai), &
                    t0, nbLoad, nbVectAsse, zi(jvVectAsse), zk24(jvVectFunc), &
                    model, materField, materCode, caraElem, loadNameJv, &
                    loadInfoJv, loadFuncJv, numedd, nume, numrep, &
                    ds_energy, sd_obsv, mesh)
!
    else if (iinteg .eq. 4) then
!
        call dladap(result, t0, lcrea, lamort, neq, &
                    jvMatr, masse, rigid, amort, zr(idepl0), &
                    zr(ivite0), zr(iacce0), zr(ifexte), zr(ifamor), zr(ifliai), &
                    nbLoad, nbVectAsse, zi(jvVectAsse), zk24(jvVectFunc), model, &
                    materField, materCode, caraElem, loadNameJv, loadInfoJv, &
                    loadFuncJv, numedd, nume, numrep, ds_energy, &
                    sd_obsv, mesh)
!
    end if
!
!====
! 7. RESULTATS
!====
!
!
    call jeveuo(result//'           .ORDR', 'L', vi=listStore)
    call jelira(result//'           .ORDR', 'LONUTI', nbStore)
    do iStore = 1, nbStore
        call rsadpa(result, 'E', 1, 'MODELE', listStore(iStore), 0, sjv=jvPara)
        zk8(jvPara) = model(1:8)
        if (materField .ne. ' ') then
            call rsadpa(result, 'E', 1, 'CHAMPMAT', listStore(iStore), 0, sjv=jvPara)
            zk8(jvPara) = materField(1:8)
        else
            call utmess('I', 'DYNALINE1_3')
        end if
        call rsadpa(result, 'E', 1, 'CARAELEM', listStore(iStore), 0, sjv=jvPara)
        zk8(jvPara) = caraElem(1:8)
    end do

! - Compute STRX_ELGA field
    call dismoi('EXI_STR2', model, 'MODELE', repk=existStr2)
    if (existStr2(1:3) .eq. 'OUI') then
        compor = materField(1:8)//'.COMPOR'
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
        lPoux = .false.
        call dismoi('EXI_POUX', model, 'MODELE', repk=existPoux)
        if (existPoux(1:3) .eq. 'OUI') then
            lPoux = .true.
            if (nbLoad .ne. 0) then
                call jeveuo(loadNameJv, 'L', jchar)
                call cochre(zk24(jchar), nbLoad, nbchre, iocc)
                if (nbchre .gt. 1) then
                    call utmess('F', 'DYNAMIQUE_19')
                end if
                if (iocc .gt. 0) then
                    call getvid('EXCIT', 'CHARGE', iocc=iocc, scal=loadPres, nbret=iret)
                    call getvid('EXCIT', 'FONC_MULT', iocc=iocc, scal=nomfon, nbret=nfon)
                    if (nfon .ne. 0) then
                        call getvr8('EXCIT', 'COEF_MULT', iocc=iocc, scal=coefMultR, nbret=ncomu)
                    end if
                end if
            end if
            if (ncomu .eq. 0) coefMultR = 1.d0
        end if
        do iStore = 0, nbStore
! --------- Get displacements
            call rsexch(' ', result, 'DEPL', iStore, disp, iret)
            if (iret .gt. 0) cycle

! --------- Prepare input fields
            call mecham('STRX_ELGA', model, numeHarm, &
                        chgeom, chharm, iret)
            if (iret .ne. 0) cycle

! --------- Prepare field for time
            call rsadpa(result, 'L', 1, 'INST', iStore, 0, sjv=jvPara)
            time = zr(jvPara)
            call mechti(chgeom(1:8), time, rundf, rundf, chtime)

! --------- Get external state varaibles
            call vrcins(model, materField, caraElem(1:8), time, chvarc(1:19), codret)
            call vrcref(model(1:8), materField(1:8), caraElem(1:8), chvref(1:19))
            if (lPoux .and. nfon .ne. 0) then
                call fointe('F ', nomfon, 1, ['INST'], [time], &
                            coefMultR, iret)
            end if

! --------- Get field to create STRX_ELGA
            call rsexch(' ', result, 'STRX_ELGA', iStore, strx, iret)
            if (iret .ne. 0) cycle

! --------- Compute STRX_ELGA
            call compStrx(model, materCode, caraElem, compor, &
                          disp, chgeom, &
                          chvarc, chvref, &
                          lPoux, loadPres, coefMultR, &
                          ligrel, jvBase, strx, iret)

! --------- Save STRX_ELGA
            call rsnoch(result, 'STRX_ELGA', iStore)
        end do
    end if
!
    call nonlinDSEnergyClean(ds_energy)
    call nonlinDSInOutClean(ds_inout)
!
    call jedema()
end subroutine
