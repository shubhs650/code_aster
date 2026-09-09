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
subroutine crtype()
!
    use listLoad_module
    use listLoad_type
!
    implicit none

!
!     COMMANDE:  CREA_RESU /AFFE
!     CREE UNE STRUCTURE DE DONNEE DE TYPE
!           "EVOL_THER"    "EVOL_VARC"       "EVOL_ELAS"
!           "MULT_ELAS"    "FOURIER_ELAS"    "FOURIER_THER"
!           "DYNA_TRANS"   "DYNA_HARMO"      "EVOL_CHAR"
!           "MODE_MECA"    "MODE_MECA_C"
!
! --- ------------------------------------------------------------------
#include "asterc/getexm.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/fonbpa.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/gnomsd.h"
#include "asterfort/idensd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/lrcomm.h"
#include "asterfort/refdaj.h"
#include "asterfort/resuGetLoads.h"
#include "asterfort/resuSaveParameters.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsmxno.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/rssepa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8) :: mxpara, ibid, ier, lg, icompt, iret, nbfac, numini, numfin
    integer(kind=8) :: n0, n1, n2, n3, nis, nbinst, ip, nbval, nume, igd, l, i, j, jc
    integer(kind=8) :: iad, jinst, nbpf
    integer(kind=8) :: ino, nbv(1), jrefe, icmpd, icmpi, nocc
    integer(kind=8) :: nbtrou, jcpt, nbr, ivmx, k, iocc, nbecd, nbeci, nboini
    integer(kind=8) :: valii(2), nfr, n4, jnmo, nmode, nbcmpd, nbcmpi, tnum(1)
    integer(kind=8) :: nbordr1, nbordr2, ier1, nbModel, nbMaterField, nbCaraElem
!
    parameter(mxpara=10)
!
    aster_logical :: lncas, lfonc, lcopy, lReuse, lLireResu
!
    real(kind=8) :: valpu(mxpara), rbid, tps, prec, valrr(3), freq, amor_red, coef(3)
    complex(kind=8) :: cbid
    character(len=4) :: typabs
    character(len=6) :: typegd, typegd2
    character(len=8) :: k8b, resultName, nomf, mesh, typmod, criter, matr, nogdsi, axe
    character(len=8) :: resultNameReuse
    character(len=8) :: model, materField, caraElem, mesh2
    character(len=8) :: modelPrev, materFieldPrev, caraElemPrev
    character(len=14) :: numedd
    character(len=16) :: nomp(mxpara), type, oper, acces, k16b
    character(len=19) :: nomch, champ, listr8, pchn1, resu19, profch
    character(len=19) :: nume_equa_tmp
    character(len=24) :: k24, linst, fieldType, resultType, lcpt, o1, o2, noojb
    character(len=24) :: valkk(4), matric(3), listLoad
    character(len=32) :: kjexn
    character(len=8), pointer :: champs(:) => null()
    real(kind=8), pointer :: coor(:) => null()
    character(len=8), pointer :: vnomf(:) => null()
    real(kind=8), pointer :: val(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    character(len=24), pointer :: prol(:) => null()
    character(len=8) :: answer
    aster_logical :: lAlarm
    aster_logical, parameter :: staticOperator = ASTER_TRUE
    integer(kind=8) :: nbLoad, iLoad, indxLoadInList
    character(len=4), parameter :: phenom = "MECA"
    character(len=16), parameter :: loadApply = "FIXE_CSTE"
    character(len=8), pointer :: loadName(:) => null()
    character(len=8) :: loadFunc
    character(len=8), parameter :: funcCste = '&&NMDOME'
    character(len=16) :: loadCommand
    character(len=13) :: loadPreObject
    aster_logical :: loadIsFunc

!
    data linst, listr8, lcpt/'&&CRTYPE_LINST', '&&CRTYPE_LISR8',&
     &     '&&CPT_CRTYPE'/
! --- ------------------------------------------------------------------
    call jemarq()
!
    listLoad = ' '
    nboini = 10
    nbModel = 0
    nbMaterField = 0
    nbCaraElem = 0
    modelPrev = ' '
    materFieldPrev = ' '
    caraElemPrev = ' '
    answer = ' '
    lLireResu = ASTER_FALSE
!
    call getres(resultName, type, oper)
    resu19 = resultName
    call getfac('AFFE', nbfac)
    call getvtx(' ', 'TYPE_RESU', scal=resultType, nbret=n1)
!
! - Reuse mode
!
    lReuse = ASTER_FALSE
    if (getexm(' ', 'RESULTAT') .eq. 1) then
        call getvid(' ', 'RESULTAT', scal=resultNameReuse, nbret=nocc)
        if (nocc .ne. 0) then
            lReuse = ASTER_TRUE
            if (resultName .ne. resultNameReuse) then
                call utmess('F', 'SUPERVIS2_79', sk='RESULTAT')
            end if
        end if
    end if
!
    call jeexin(resultName//'           .DESC', iret)
    if (iret .eq. 0) call rscrsd('G', resultName, resultType, nboini)
!
    call jelira(resultName//'           .ORDR', 'LONUTI', nbordr1)
!
    lncas = .false.
    if (resultType .eq. 'MULT_ELAS' .or. resultType .eq. 'FOURIER_ELAS' .or. resultType .eq. &
        'FOURIER_THER' .or. resultType .eq. 'MODE_MECA' .or. resultType .eq. 'MODE_MECA_C') then
        lncas = .true.
    end if
!
    numini = -1
    icompt = -1
    profch = ' '
    AS_ALLOCATE(vk8=champs, size=nbfac)
!
    do iocc = 1, nbfac
        model = ' '

        call getvtx('AFFE', 'NOM_CHAM', iocc=iocc, scal=fieldType, nbret=n1)
!
!   on compte les modeles, materiaux et les cara_ele différents d'un pas à l'autre
!   (y compris la chaine ' ' )
!   si on en trouve au moins 2 différents, l'appel final à lrcomm se fera avec ' '
!
! ----- Count number of different model/caraElem/materField
        model = " "
        call getvid('AFFE', 'MODELE', iocc=iocc, scal=model, nbret=n1)
        if (model .ne. ' ' .and. model .ne. modelPrev) then
            nbModel = nbModel+1
        end if
        modelPrev = model

        materField = " "
        call getvid('AFFE', 'CHAM_MATER', iocc=iocc, scal=materField, nbret=n1)
        if (materField .ne. ' ' .and. materField .ne. materFieldPrev) then
            nbMaterField = nbMaterField+1
        end if
        materFieldPrev = materField

        caraElem = " "
        call getvid('AFFE', 'CARA_ELEM', iocc=iocc, scal=caraElem, nbret=n1)
        if (caraElem .ne. ' ' .and. caraElem .ne. caraElemPrev) then
            nbCaraElem = nbCaraElem+1
        end if
        caraElemPrev = caraElem

! ----- Get current loads
        call getvid('AFFE', 'CHARGE', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .lt. 0) then
! --------- Generate name of datastructure to save in result datastructure
            call nameListLoad(listLoad)

! --------- Generate constant function
            call createUnitFunc(funcCste, "G", loadFunc)

! --------- Create list of loads datastructure
            nbLoad = -n1
            call creaListLoad("MECA", "G", nbLoad, listLoad)
            if (nbLoad .ne. 0) then
                AS_ALLOCATE(vk8=loadName, size=nbLoad)
                call getvid('AFFE', 'CHARGE', iocc=iocc, nbval=nbLoad, vect=loadName)
                indxLoadInList = 0
                do iLoad = 1, nbLoad
                    call getLoadParameters(phenom, model, loadName(iLoad), &
                                           loadPreObject, loadCommand, loadIsFunc)
                    call addLoadMeca(staticOperator, listLoad, &
                                     loadName(iLoad), loadFunc, &
                                     loadApply, loadCommand, loadPreObject, &
                                     loadIsFunc, &
                                     indxLoadInList)
                end do
                AS_DEALLOCATE(vk8=loadName)
            end if
        end if

! ----- Get current field
        call getvid('AFFE', 'CHAM_GD', iocc=iocc, scal=champ, nbret=n1)
        champs(iocc) = champ(1:8)
        call dismoi('NOM_MAILLA', champ, 'CHAMP', repk=mesh)
        if (model .ne. ' ') then
            call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh2)
            if (mesh .ne. mesh2) then
                call utmess('F', 'ALGORITH2_1')
            end if
        end if
        call dismoi('NOM_GD', champ, 'CHAMP', repk=nogdsi)
        if (resultType .eq. 'EVOL_CHAR' .and. nogdsi .eq. 'NEUT_R') then
            valkk(1) = 'NEUT_R'
            valkk(2) = 'EVOL_CHAR'
            call utmess('F', 'ALGORITH2_80', nk=2, valk=valkk)
        end if
!
        call dismoi('TYPE_SUPERVIS', champ, 'CHAMP', repk=k24)
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=coor)
!
!        CALCUL DE LFONC ET TYPEGD
        lfonc = .false.
        do i = 24, 1, -1
            if (k24(i:i) .eq. ' ') goto 10
            if (k24(i-1:i) .eq. '_F') then
                if (k24(1:7) .ne. 'CHAM_NO') then
                    call utmess('F', 'ALGORITH2_45', sk=k24)
                end if
                lfonc = .true.
                typegd = k24(i-5:i-2)//'_R'
            else if (k24(i-1:i) .eq. '_R') then
                typegd = k24(i-5:i)
            else if (k24(i-1:i) .eq. '_C') then
                typegd = k24(i-5:i)
            else
                call utmess('F', 'ALGORITH2_46', sk=k24)
            end if
            goto 20
10          continue
        end do
20      continue
!
        if (k24(1:7) .eq. 'CHAM_NO') then
!           -- on cherche a economiser les nume_equa (partage si possible)
            if (profch .eq. ' ') then
                call dismoi('NUME_EQUA', champ, 'CHAM_NO', repk=pchn1)
                noojb = '12345678.NUMEQ00000.PRNO'
                call gnomsd(' ', noojb, 15, 19)
                profch = noojb(1:19)
                lcopy = .true.
!
                if (lcopy) then
                    call copisd('NUME_EQUA', 'G', pchn1, profch)
                end if
            else
                call dismoi('NUME_EQUA', champ, 'CHAM_NO', repk=pchn1)
                if (.not. idensd('NUME_EQUA', profch, pchn1)) then
                    noojb = '12345678.NUMEQ00000.PRNO'
                    call gnomsd(' ', noojb, 15, 19)
                    profch = noojb(1:19)
                    call copisd('NUME_EQUA', 'G', pchn1, profch)
                end if
            end if

            if (lfonc) then
                call dismoi('NOM_GD', profch, 'NUME_EQUA', repk=typegd2)
                if (typegd .ne. typegd2) then
                    noojb = '12345678.NUMEC00000.PRNO'
                    call gnomsd(' ', noojb, 15, 19)
                    nume_equa_tmp = noojb(1:19)
                    call copisd('NUME_EQUA', 'G', profch, nume_equa_tmp)
                    profch = nume_equa_tmp
                    call jeveuo(profch//'.REFN', 'E', jrefe)
                    zk24(jrefe-1+2) = typegd
                end if
            end if
        end if
!
!        MOT CLE "NOM_CAS", "NUME_MODE", "FREQ"  PRESENT :
        if (lncas) then
            call rsorac(resultName, 'LONUTI', 0, rbid, k8b, &
                        cbid, rbid, k8b, tnum, 1, &
                        nbtrou)
            numini = tnum(1)
            j = 0
            if (resultType .eq. 'MODE_MECA') then
                call getvis('AFFE', 'NUME_MODE', iocc=iocc, scal=nume, nbret=n0)
                if (n0 .ne. 0) then
                    do i = 1, numini
                        call rsadpa(resultName, 'L', 1, 'NUME_MODE', i, &
                                    0, sjv=jnmo, styp=k8b)
                        nmode = zi(jnmo)
                        if (nmode .eq. nume) then
                            numini = nume
                            j = j+1
                        end if
                    end do
                end if
            else if (resultType .eq. 'MULT_ELAS') then
                call getvtx('AFFE', 'NOM_CAS', iocc=iocc, scal=acces, nbret=n0)
                if (n0 .gt. 0) then
                    call rsorac(resultName, 'NOM_CAS', ibid, rbid, acces, &
                                cbid, 1.d0, 'ABSOLU', tnum, 1, &
                                nbr)
                    if (nbr .ne. 0) then
                        numini = tnum(1)
                        j = j+1
                    end if
                end if
            end if
!           not found
            if (j .eq. 0) then
                numini = numini+1
            end if
            !
            call rsexch(' ', resultName, fieldType, numini, nomch, &
                        iret)
            if (iret .eq. 0) then
                valkk(1) = champ(1:8)
                valii(1) = numini
                call utmess('A', 'ALGORITH12_74', sk=valkk(1), si=valii(1))
            else if (iret .eq. 110) then
                call rsagsd(resultName, 0)
                call rsexch(' ', resultName, fieldType, numini, nomch, &
                            iret)
            else if (iret .eq. 100) then
!              ON NE FAIT RIEN
            else
                call utmess('F', 'ALGORITH2_47', sk=fieldType)
            end if
!
            call copisd('CHAMP_GD', 'G', champ, nomch)
            if (k24(1:7) .eq. 'CHAM_NO') then
                call dismoi('NUME_EQUA', nomch, 'CHAM_NO', repk=pchn1)
                if (pchn1 .ne. profch) then
                    call detrsd('NUME_EQUA', pchn1)
                    call jeveuo(nomch//'.REFE', 'E', jrefe)
                    zk24(jrefe+1) = profch
                end if
            end if
!
            call rsnoch(resultName, fieldType, numini)
            call rssepa(resultName, numini, model, materField, caraElem, listLoad)
!
            call getvtx('AFFE', 'NOM_CAS', iocc=iocc, scal=acces, nbret=n0)
            if (n0 .ne. 0) then
                call rsadpa(resultName, 'E', 1, 'NOM_CAS', numini, &
                            0, sjv=iad, styp=k8b)
                zk16(iad) = acces
            end if
!
            call getvis('AFFE', 'NUME_MODE', iocc=iocc, scal=nume, nbret=n0)
            if (n0 .ne. 0) then
                call rsadpa(resultName, 'E', 1, 'NUME_MODE', numini, &
                            0, sjv=iad, styp=k8b)
                zi(iad) = nume
            end if
!
            call getvtx('AFFE', 'TYPE_MODE', iocc=iocc, scal=typmod, nbret=n0)
            if (n0 .ne. 0) then
                call rsadpa(resultName, 'E', 1, 'TYPE_MODE', numini, &
                            0, sjv=iad, styp=k8b)
                zk8(iad) = typmod
            end if
!
            call getvr8('AFFE', 'FREQ', iocc=iocc, scal=freq, nbret=n0)
            if (n0 .ne. 0) then
                call rsadpa(resultName, 'E', 1, 'FREQ', numini, &
                            0, sjv=iad, styp=k8b)
                zr(iad) = freq
!               HERE ONE IS IN THE CASE 'MODE_MECA' or 'MODE_MECA_C'
!               SO IF A FREQUENCY IS GIVEN, ONE CONSIDER THAT THE GIVEN CHAM_GD
!               IS A MODAL SHAPE
!               (IN OPPOSITION WITH A STATIC DEFORMED SHAPE)
                call rsadpa(resultName, 'E', 1, 'TYPE_DEFO', numini, &
                            0, sjv=iad, styp=k8b)
                zk16(iad) = 'PROPRE'
            end if
!           pour COMB_SISM_MODAL/MODE_CORR
            call getvtx('AFFE', 'AXE', iocc=iocc, scal=axe, nbret=n0)
            if (n0 .ne. 0) then
                call rsadpa(resultName, 'E', 1, 'NOEUD_CMP', numini, &
                            0, sjv=iad, styp=k8b)
                if (axe(1:1) .eq. 'X') then
                    zk16(iad) = 'ACCE    X       '
                    coef(1) = 1.d0
                    coef(2) = 0.d0
                    coef(3) = 0.d0
                else if (axe(1:1) .eq. 'Y') then
                    zk16(iad) = 'ACCE    Y       '
                    coef(1) = 0.d0
                    coef(2) = 1.d0
                    coef(3) = 0.d0
                else if (axe(1:1) .eq. 'Z') then
                    zk16(iad) = 'ACCE    Z       '
                    coef(1) = 0.d0
                    coef(2) = 0.d0
                    coef(3) = 1.d0
                end if
                call rsadpa(resultName, 'E', 1, 'COEF_X', numini, &
                            0, sjv=iad, styp=k8b)
                zr(iad) = coef(1)
                call rsadpa(resultName, 'E', 1, 'COEF_Y', numini, &
                            0, sjv=iad, styp=k8b)
                zr(iad) = coef(2)
                call rsadpa(resultName, 'E', 1, 'COEF_Z', numini, &
                            0, sjv=iad, styp=k8b)
                zr(iad) = coef(3)
                call rsadpa(resultName, 'E', 1, 'TYPE_DEFO', numini, &
                            0, sjv=iad, styp=k8b)
                zk16(iad) = 'ACCE_IMPO'
            end if
!
            call getvr8('AFFE', 'AMOR_REDUIT', iocc=iocc, scal=amor_red, nbret=n0)
            if (n0 .ne. 0) then
                call rsadpa(resultName, 'E', 1, 'AMOR_REDUIT', numini, &
                            0, sjv=iad, styp=k8b)
                zr(iad) = amor_red
            end if
            goto 80
        end if
!
!        MOT CLE INST/FREQ PRESENT :
        nis = 0
        nfr = 0
        nbinst = 0
        call getvr8('AFFE', 'INST', iocc=iocc, nbval=0, nbret=nis)
        call getvr8('AFFE', 'FREQ', iocc=iocc, nbval=0, nbret=nfr)
        if (nis .ne. 0) then
            typabs = 'INST'
            nbinst = -nis
        end if
        if (nfr .ne. 0) then
            typabs = 'FREQ'
            nbinst = -nfr
        end if
!
        if ((nis .ne. 0) .or. (nfr .ne. 0)) then
            call wkvect(lcpt, 'V V I', nbinst, jcpt)
            call wkvect(linst, 'V V R', nbinst, jinst)
            call getvr8('AFFE', typabs, iocc=iocc, nbval=nbinst, vect=zr(jinst), &
                        nbret=n1)
            call getvr8('AFFE', 'PRECISION', iocc=iocc, scal=prec, nbret=ibid)
            call getvtx('AFFE', 'CRITERE', iocc=iocc, scal=criter, nbret=ibid)
            call rsorac(resultName, 'LONUTI', 0, rbid, k8b, &
                        cbid, rbid, k8b, nbv, 1, &
                        ibid)
!
            ivmx = rsmxno(resultName)
            do k = 1, nbinst
                if (nbv(1) .gt. 0) then
                    call rsorac(resultName, typabs, ibid, zr(jinst+k-1), k8b, &
                                cbid, prec, criter, tnum, 1, &
                                nbr)
                    nume = tnum(1)
                else
                    nbr = 0
                end if
                if (nbr .lt. 0) then
                    call utmess('F', 'ALGORITH2_48')
                else if (nbr .eq. 0) then
                    zi(jcpt+k-1) = ivmx+1
                    ivmx = ivmx+1
                else
                    zi(jcpt+k-1) = nume
                end if
            end do
        else
!           MOT CLE LIST_INST/LIST_FREQ PRESENT :
            n1 = 0
            n4 = 0
            call getvid('AFFE', 'LIST_INST', iocc=iocc, scal=listr8, nbret=n1)
            call getvid('AFFE', 'LIST_FREQ', iocc=iocc, scal=listr8, nbret=n4)
            if (n1 .ne. 0) then
                typabs = 'INST'
            end if
            if (n4 .ne. 0) then
                typabs = 'FREQ'
            end if
!
            call getvr8('AFFE', 'PRECISION', iocc=iocc, scal=prec, nbret=ibid)
            call getvtx('AFFE', 'CRITERE', iocc=iocc, scal=criter, nbret=ibid)
            call jelira(listr8//'.VALE', 'LONMAX', nbval)
!
            nbinst = nbval
            numini = 1
            numfin = nbinst
            call getvis('AFFE', 'NUME_INIT', iocc=iocc, scal=numini, nbret=n2)
            call getvis('AFFE', 'NUME_FIN', iocc=iocc, scal=numfin, nbret=n3)
            if (numfin .gt. nbval) numfin = nbval
            if (n2 .ne. 0 .and. n3 .ne. 0) then
                if (numfin .lt. numini) then
                    call utmess('F', 'ALGORITH2_49')
                end if
                nbinst = numfin-numini+1
!
            else if (n2 .ne. 0) then
                nbinst = nbval-numini+1
            else if (n3 .ne. 0) then
                nbinst = numfin
            else
                nbinst = nbval
            end if
            nbinst = min(nbinst, nbval)
!
            call wkvect(linst, 'V V R', nbinst, jinst)
            call jeveuo(listr8//'.VALE', 'L', vr=val)
            call rsorac(resultName, 'LONUTI', 0, rbid, k8b, &
                        cbid, rbid, k8b, nbv, 1, &
                        ibid)
            call wkvect(lcpt, 'V V I', nbinst, jcpt)
            ivmx = rsmxno(resultName)
            j = 0
            do k = 1, nbval
                if (k .lt. numini) goto 40
                if (k .gt. numfin) goto 40
                j = j+1
                zr(jinst-1+j) = val(k)
                if (nbv(1) .gt. 0) then
                    call rsorac(resultName, typabs, ibid, val(k), k8b, &
                                cbid, prec, criter, tnum, 1, &
                                nbr)
                    nume = tnum(1)
                else
                    nbr = 0
                end if
                if (nbr .lt. 0) then
                    call utmess('F', 'ALGORITH2_48')
                else if (nbr .eq. 0) then
                    zi(jcpt+j-1) = ivmx+1
                    ivmx = ivmx+1
                else
                    zi(jcpt+j-1) = nume
                end if
40              continue
            end do
        end if
!
!        DANS LE CAS DES FONCTIONS, LA PROGRAMMATION N'EST VALABLE QUE
!        SI POUR LES GRANDEURS XXXX_F ET YYYY_R :
!           * XXXX = YYYY
!           * ONT LE MEME NOMBRE D'ENTIERS CODES
!           * QUE LE RANG (DANS LE CATALOGUE) DE CHAQUE CMP
!             DE XXXX_F SOIT LE MEME QUE DANS YYYY_R
        if (lfonc) then
!           POUR EVOL_VARC : MEME GRANDEUR XXXX ET SOUS NOM_CHAM
            if (resultType .eq. 'EVOL_VARC') then
                if (fieldType(1:4) .ne. nogdsi(1:4)) then
                    valkk(1) = fieldType(1:4)
                    valkk(2) = nogdsi(1:4)
                    call utmess('F', 'CALCULEL2_79', nk=2, valk=valkk)
                end if
            end if
!           DANS TOUS LES AUTRES CAS, MEME GRANDEUR XXXX = YYYY
            if (typegd(1:4) .ne. nogdsi(1:4)) then
                valkk(1) = typegd(1:4)
                valkk(2) = nogdsi(1:4)
                call utmess('F', 'CALCULEL2_90', nk=2, valk=valkk)
            end if
!           NOMBRE D'ENTIER CODE
            call dismoi('NB_EC', typegd, 'GRANDEUR', repi=nbecd)
            call dismoi('NB_EC', nogdsi, 'GRANDEUR', repi=nbeci)
            if (nbecd .ne. nbeci) then
                valkk(1) = typegd
                valkk(2) = nogdsi
                valii(1) = nbecd
                valii(2) = nbeci
                call utmess('F', 'CALCULEL2_80', nk=2, valk=valkk, ni=2, &
                            vali=valii)
            end if
!           NOM DES COMPOSANTES DU MEME RANG IDENTIQUE
            kjexn = jexnom('&CATA.GD.NOMCMP', typegd)
            call jeveuo(kjexn, 'L', icmpd)
            call jelira(kjexn, 'LONMAX', nbcmpd)
!
            kjexn = jexnom('&CATA.GD.NOMCMP', nogdsi)
            call jeveuo(kjexn, 'L', icmpi)
            call jelira(kjexn, 'LONMAX', nbcmpi)
            do j = 1, nbcmpi
                if (zk8(icmpi+j-1) .ne. zk8(icmpd+j-1)) then
                    valkk(1) = typegd
                    valkk(2) = nogdsi
                    valkk(3) = zk8(icmpd+j-1)
                    valkk(4) = zk8(icmpi+j-1)
                    call utmess('F', 'CALCULEL2_5', nk=4, valk=valkk)
                end if
            end do
        end if
!
        do j = 1, nbinst
            if (j .ge. 2) call jemarq()
            call jerecu('V')
            icompt = zi(jcpt+j-1)
            tps = zr(jinst+j-1)
            call rsexch(' ', resultName, fieldType, icompt, nomch, &
                        iret)
            if (iret .eq. 0) then
                call rsadpa(resultName, 'L', 1, typabs, icompt, &
                            0, sjv=iad, styp=k8b)
                valkk(1) = fieldType
                valkk(2) = champ(1:8)
                valrr(1) = zr(iad)
                valrr(2) = tps
                valrr(3) = prec
                call utmess('A', 'ALGORITH11_87', nk=2, valk=valkk, nr=3, &
                            valr=valrr)
            else if (iret .eq. 110) then
                call rsagsd(resultName, 0)
                call rsexch(' ', resultName, fieldType, icompt, nomch, &
                            iret)
            end if
!
            if (k24(1:7) .eq. 'CHAM_NO') then
!
                o1 = champ//'.REFE'
                o2 = nomch//'.REFE'
                call jedupo(o1, 'G', o2, .false._1)
!
                o1 = champ//'.VALE'
                o2 = nomch//'.VALE'
                call jedupo(o1, 'G', o2, .false._1)
!
                call jeveuo(nomch//'.REFE', 'E', jrefe)
                zk24(jrefe+1) = profch
            else
                call copisd('CHAMP_GD', 'G', champ, nomch)
            end if
!
            if (lfonc) then
                call jelira(champ//'.VALE', 'LONMAX', lg)
                call jeveuo(champ//'.VALE', 'L', vk8=vnomf)
                call jeveuo(champ//'.REFE', 'L', jrefe)
                call jeveuo(zk24(jrefe+1) (1:19)//'.DEEQ', 'L', vi=deeq)
!
                if (k24(1:7) .ne. 'CHAM_NO') then
                    call jeveuo(nomch//'.DESC', 'E', vi=desc)
                    call jenonu(jexnom('&CATA.GD.NOMGD', typegd), igd)
                    desc(1) = igd
                end if
                call jedetr(nomch//'.VALE')
                call wkvect(nomch//'.VALE', 'G V R', lg, jc)
!              CHAM_NO DE FONCTIONS A EVALUER
                call jeveuo(nomch//'.VALE', 'E', jc)
                do l = 1, lg
                    nomf = vnomf(l)
                    if (nomf .eq. ' ') goto 60
                    call jeveuo(nomf//'           .PROL', 'L', vk24=prol)
                    call fonbpa(nomf, prol, k16b, mxpara, nbpf, &
                                nomp)
                    ino = deeq(1+2*(l-1))
                    if (ino .eq. 0) goto 60
                    do ip = 1, nbpf
                        if (nomp(ip) .eq. 'INST') then
                            valpu(ip) = tps
                        else if (nomp(ip) .eq. 'X') then
                            valpu(ip) = coor(3*(ino-1)+1)
                        else if (nomp(ip) .eq. 'Y') then
                            valpu(ip) = coor(3*(ino-1)+2)
                        else if (nomp(ip) .eq. 'Z') then
                            valpu(ip) = coor(3*(ino-1)+3)
                        else
                            call utmess('F', 'ALGORITH2_50')
                        end if
                    end do
                    call fointe('F', nomf, nbpf, nomp, valpu, &
                                zr(jc+l-1), ier)
60                  continue
                end do
            end if
!
            call rsnoch(resultName, fieldType, icompt)
            call rsadpa(resultName, 'E', 1, typabs, icompt, 0, sjv=iad)
            zr(iad) = tps
            call rssepa(resultName, icompt, model, materField, caraElem, listLoad)
            if (j .ge. 2) call jedema()
!
        end do
        call jedetr(linst)
        call jedetr(lcpt)
80      continue
    end do

!
!     REMPLISSAGE DE .REFD POUR LES MODE_MECA  ET DYNA_*:
    call jelira(resultName//'           .ORDR', 'LONUTI', nbordr2)
    if (nbordr2 .gt. nbordr1) then

        if (resultType(1:9) .eq. 'MODE_MECA' &
            .or. resultType(1:10) .eq. 'DYNA_HARMO' &
            .or. resultType(1:10) .eq. 'DYNA_TRANS') then

            matric(1) = ' '
            matric(2) = ' '
            matric(3) = ' '
            numedd = ' '
            call getvid(' ', 'MATR_RIGI', scal=matr, nbret=n1)
            if (n1 .eq. 1) then
                call dismoi('NOM_NUME_DDL', matr, 'MATR_ASSE', repk=numedd)
                matric(1) = matr
            else
                call getvid(' ', 'MATR_MASS', scal=matr, nbret=n1)
                if (n1 .eq. 1) then
                    call dismoi('NOM_NUME_DDL', matr, 'MATR_ASSE', repk=numedd)
                end if
            end if
            call getvid(' ', 'MATR_MASS', scal=matr, nbret=n1)
            if (n1 .eq. 1) then
                matric(2) = matr
            end if
!           If no numbering information could be found, try to retrieve the
!           information from the fields composing the sd_resultat
            if (numedd .eq. ' ') then
                call getvid('AFFE', 'CHAM_GD', iocc=1, scal=champ, nbret=ier)
                call dismoi('NUME_EQUA', champ, 'CHAMP', repk=profch, arret='C', &
                            ier=ier)
                if (ier .eq. 0) then
                    call refdaj('F', resu19, (nbordr2-nbordr1), profch, 'DYNAMIQUE', matric, ier)
                end if
            else
                call refdaj('F', resu19, (nbordr2-nbordr1), numedd, 'DYNAMIQUE', &
                            matric, ier)
!
!               compare numedd and numeequa of all the new fields added (only DEPL)
!
                do j = nbordr1+1, nbordr2-nbordr1
                    call rsexch(' ', resu19, 'DEPL', j, nomch, ier1)
                    if (ier1 .eq. 0) then
                        call dismoi('NUME_EQUA', nomch, 'CHAMP', repk=profch, &
                                    arret='C', ier=ier)
                        if (ier .eq. 0) then
                            if (.not. idensd('NUME_EQUA', numedd(1:14)//'.NUME', profch)) then
                                valkk(1) = numedd
                                valkk(2) = profch
                                call utmess('A', 'ALGORITH2_51', nk=2, valk=valkk)
                            end if
                        end if
                    end if
                end do
            end if
        end if
    end if
!
    if (resultType .eq. 'EVOL_NOLI' .or. resultType .eq. 'EVOL_ELAS' .or. &
        resultType .eq. 'EVOL_THER') then
        if (nbModel .gt. 1) model = ' '
        if (nbMaterField .gt. 1) materField = ' '
        if (nbCaraElem .gt. 1) caraElem = ' '
! ----- Get loads/BC and create list of loads
        call resuGetLoads(model, resultType, listLoad)

! ----- Save standard parameters in results datastructure
        call resuSaveParameters(resultName, resultType, &
                                model, caraElem, materField, listLoad)
    end if

! - Non-linear behaviour management
    if (resultType .eq. 'EVOL_NOLI') then
        call getvtx(' ', 'VERI_VARI', scal=answer, nbret=n1)
        lAlarm = answer .eq. 'OUI'
        call lrcomm(lReuse, resultName, model, caraElem, materField, lLireResu, lAlarm)
    end if
!
    AS_DEALLOCATE(vk8=champs)
    call jedema()
end subroutine
