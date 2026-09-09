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
subroutine modirepresu(resultOutZ, resultIn)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesfus.h"
#include "asterfort/chrpel.h"
#include "asterfort/chrpno.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/gettco.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/medome_once.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsinfo.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsorac.h"
#include "asterfort/rsutnu.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: resultOutZ, resultIn
!
! --------------------------------------------------------------------------------------------------
!
!     COMMANDE : MODI_REPERE / RESULTAT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n0, n1, nbordr, iret, noccField, i, j, np, iordr, ndim
    integer(kind=8) :: iord, ioccField, ibid, nc
    integer(kind=8) :: jordr, nbnosy, jpa, iadin, iadou
    integer(kind=8) :: nbpara, nbac, nbpa, ifm, niv, nncp
    real(kind=8) :: prec
    real(kind=8) :: lcoer(2)
    real(kind=8) :: r8b
    complex(kind=8) :: c16b
    complex(kind=8) :: lcoec(2)
    character(len=8) :: k8b
    character(len=8) :: crit, tych, nomma, model, modelRefe
    character(len=8) :: caraElem, exipla, exicoq
    character(len=16) :: fieldName, tysd, type, fieldDime, repere, option1, optinit
    character(len=16) :: cham_resu
    character(len=19) :: knum, resultOut
    character(len=19) :: chams1, chams0, chafus, chs(2), ligrelField, ligrelCalc
    character(len=24) :: nompar, fieldIn, champ01, fieldOut, champ2
    character(len=24) :: valk(2)
    integer(kind=8), pointer :: nume_ordre(:) => null()
!
    aster_logical :: lreuse, lcumu(2), lcoc(2), lModelVariable, check
!
    data lcumu/.false., .false./
    data lcoc/.false., .false./
    data lcoer/1.d0, 1.d0/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
!   LE CONCEPT EST REENTRANT SI REPERE = 'COQUE_INTR_UTIL' OU 'COQUE_UTIL_INTR' ou 'COQUE_UTIL_CYL'
!   DANS CE CAS ON CREE UNE SD RESULTAT TEMPORAIRE POUR LES CALCULS ET ENSUITE ON SURCHARGE
!   RESUIN PAR LES CHAMPS MODIFIES STOCKES DANS RESUOU
    lreuse = .false.
    resultOut = resultOutZ
    if (resultIn .eq. resultOut) then
        lreuse = .true.
        resultOut = 'MODIREPE'
    end if
!
    call jelira(resultIn//'.DESC', 'NOMMAX', nbnosy)
    if (nbnosy .eq. 0) goto 999
!
    call gettco(resultIn, tysd)
!   RECUPERATION DU NOMBRE DE CHAMPS SPECIFIE
    call getfac('MODI_CHAM', noccField)
!
!   DEFINITION DU REPERE UTILISE
    call getvtx(' ', 'REPERE', scal=repere, nbret=i)

    if (lreuse) then
        if (i .eq. 0) then
            call utmess('F', 'MODELISA3_14')
        end if
        if ((repere .ne. 'COQUE_INTR_UTIL') .and. &
            (repere .ne. 'COQUE_UTIL_INTR') .and. &
            (repere .ne. 'UTILISATEUR')) then
            call utmess('F', 'MODELISA3_15', nk=1, valk=repere)
        end if
    end if
!
!   RECUPERATION DES NUMEROS D'ORDRE DE LA STRUCTURE DE DONNEES DE TYPE RESULTAT RESU A PARTIR
!   DES VARIABLES D'ACCES UTILISATEUR 'NUME_ORDRE','FREQ','INST','NOEUD_CMP'
!   (VARIABLE D'ACCES 'TOUT_ORDRE' PAR DEFAUT)
    knum = '&&OP0191.NUME_ORDRE'
    call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
    call rsutnu(resultIn, ' ', 1, knum, nbordr, prec, crit, iret)
    if (iret .eq. 10) then
        call utmess('F', 'CALCULEL4_8', sk=resultIn)
    end if
    if (iret .ne. 0) then
        call utmess('F', 'ALGORITH3_41')
    end if
    call jeveuo(knum, 'L', jordr)
    call rscrsd('G', resultOut, tysd, nbordr)
!
    ! QUELQUES INITIALISATIONS
    lModelVariable = ASTER_FALSE
    modelRefe = " "
    ligrelCalc = " "
    option1 = 'EGRU_ELNO       '
    optinit = ' '
    n1 = 0

    ! POUR REPERE = COQUE_* VERIFIER QU'ON N'EST PAS EN MULTI MODELES
    if (repere(1:5) .eq. 'COQUE') then
        AS_ALLOCATE(vi=nume_ordre, size=nbordr)
        call rsorac(resultIn, 'TOUT_ORDRE', ibid, r8b, k8b, &
                    c16b, r8b, k8b, nume_ordre, nbordr, &
                    ndim)
        call medome_once(resultIn, nume_ordre, nbordr)
        AS_DEALLOCATE(vi=nume_ordre)
    end if

    do ioccField = 1, noccField
        call getvtx('MODI_CHAM', 'NOM_CHAM', iocc=ioccField, scal=fieldName, nbret=n0)
        call getvtx('MODI_CHAM', 'TYPE_CHAM', iocc=ioccField, scal=fieldDime, nbret=n0)
        call getvtx('MODI_CHAM', 'NOM_CHAM_RESU', iocc=ioccField, scal=cham_resu, nbret=n1)
        optinit = fieldName
        check = .false.
        do iord = 1, nbordr
            ! SI OPTION CHANGE DANS LE CASE DE EGRU_ELNO
            fieldName = optinit
            !
            call jemarq()
            call jerecu('V')
            iordr = zi(jordr-1+iord)
            call rsexch('F', resultIn, fieldName, iordr, fieldIn, iret)
            call dismoi('NOM_MAILLA', fieldIn(1:19), 'CHAMP', repk=nomma)
            call dismoi('TYPE_CHAMP', fieldIn, 'CHAMP', repk=tych, arret='C', ier=iret)

            !
            if ((fieldName(1:9) .eq. option1(1:9)) .and. (repere(1:15) .eq. 'COQUE_UTIL_INTR')) then
                valk(1) = fieldName(1:9)
                valk(2) = 'COQUE_UTIL_INTR'
                call utmess('F', 'ALGORITH5_87', nk=2, valk=valk)
            end if

            ! TRAITEMENT SPECIFIQUE POUR EFGE_ELNO
            if ((fieldName(1:9) .eq. 'EFGE_ELNO') .and. (n1 .gt. 0)) then
                ASSERT(cham_resu(1:9) .eq. 'EGRU_ELNO')
                ! CREER LE CHAMP EGRU
                call rsexch(' ', resultOut, fieldName, iordr, champ2, iret)
                call copisd('CHAMP_GD', 'G', fieldIn, champ2)
                call rsnoch(resultOut, fieldName, iordr)
                if (lreuse) then
                    call rsexch(' ', resultIn, option1, iordr, champ01, iret)
                    check = .true.
                end if
                fieldName = option1
            end if

            ! CHAMP1 SERA ENSUITE RECREE SUR LA BASE GLOBALE
            call rsexch(' ', resultOut, fieldName, iordr, fieldOut, iret)
            call copisd('CHAMP_GD', 'G', fieldIn, fieldOut)

!           RECUPERATION DU MODELE ASSOCIE AU CHAMP
            model = ''; caraElem = ''
            call rslesd(resultIn(1:8), iordr, model_=model, cara_elem_=caraElem)
            if (iord .eq. 1) then
                modelRefe = model
            else
                if (modelRefe .ne. model) then
                    lModelVariable = ASTER_TRUE
                end if
            end if
            if (model .ne. '') then
                call dismoi('EXI_PLAQUE', model, 'MODELE', repk=exipla)
                call dismoi('EXI_COQUE', model, 'MODELE', repk=exicoq)
                if (((exipla(1:3) .eq. 'OUI') .or. (exicoq(1:3) .eq. 'OUI')) .and. &
                    ((fieldDime .eq. 'TENS_2D') .or. (fieldDime .eq. 'TENS_3D')) .and. &
                    (repere .eq. 'UTILISATEUR')) then
                    call utmess('F', 'ALGORITH3_7')
                end if
            end if
!               Obligatoire : modèle , cara_elem
!                             repere = UTILISATEUR

            if (fieldDime .eq. '1D_GENE') then
                if ((model .eq. '') .or. (caraElem .eq. '') .or. (repere .ne. 'UTILISATEUR')) then
                    call utmess('F', 'ALGORITH2_32')
                end if
            end if
!
!           RECUPERATION DE LA NATURE DES CHAMPS (CHAM_NO OU CHAM_ELEM)
            if (tych(1:4) .eq. 'NOEU') then
                call chrpno(fieldOut, repere, fieldName, fieldDime)
            else if (tych(1:2) .eq. 'EL') then
                if (iord .eq. 1 .or. modelRefe .ne. model) then
                    call exlima('MODI_CHAM', ioccField, 'G', model, ligrelCalc)
                end if
                call chrpel(fieldOut, repere, fieldName, ioccField, fieldDime, &
                            model, caraElem, ligrelCalc, lModelVariable)
            else
                valk(1) = tych
                valk(2) = fieldOut
                call utmess('A', 'ALGORITH9_69', nk=2, valk=valk)
            end if
            call rsnoch(resultOut, fieldName, iordr)
            ! COPIER CHAMP1 SI EGRU_ELNO ET REUSE
            if (check) then
                call copisd('CHAMP_GD', 'G', fieldOut, champ01)
                call rsnoch(resultIn, fieldName, iordr)
            end if
            call jedema()
        end do
    end do

!
    nompar = '&&OP0191.NOMS_PARA'
    call rsnopa(resultIn, 2, nompar, nbac, nbpa)
    nbpara = nbac+nbpa
    call jeveuo(nompar, 'L', jpa)
    do iord = 1, nbordr
        iordr = zi(jordr-1+iord)
        do j = 1, nbpara
            call rsadpa(resultIn, 'L', 1, zk16(jpa+j-1), iordr, 1, sjv=iadin, styp=type, istop=0)
            call rsadpa(resultOut, 'E', 1, zk16(jpa+j-1), iordr, 1, sjv=iadou, styp=type)
            if (type(1:1) .eq. 'I') then
                zi(iadou) = zi(iadin)
            else if (type(1:1) .eq. 'R') then
                zr(iadou) = zr(iadin)
            else if (type(1:1) .eq. 'C') then
                zc(iadou) = zc(iadin)
            else if (type(1:3) .eq. 'K80') then
                zk80(iadou) = zk80(iadin)
            else if (type(1:3) .eq. 'K32') then
                zk32(iadou) = zk32(iadin)
            else if (type(1:3) .eq. 'K24') then
                zk24(iadou) = zk24(iadin)
            else if (type(1:3) .eq. 'K16') then
                zk16(iadou) = zk16(iadin)
            else if (type(1:2) .eq. 'K8') then
                zk8(iadou) = zk8(iadin)
            end if
        end do
    end do
!
    call titre()
    if (niv .eq. 2) call rsinfo(resultOut, ifm)
!
999 continue
!
!   CREATION DE L'OBJET .REFD SI NECESSAIRE :
    call refdcp(resultIn, resultOut)
!
!   Traitement du cas ou il y a reentrance
!       REPERE = 'COQUE_INTR_UTIL' ou 'COQUE_UTIL_INTR' ou 'UTILISATEUR'
    if (lreuse) then
        do ioccField = 1, noccField
            call getvtx('MODI_CHAM', 'NOM_CHAM', iocc=ioccField, scal=fieldName, nbret=n0)
            do iord = 1, nbordr
                call jemarq()
                call jerecu('V')
                iordr = zi(jordr-1+iord)
                call rsexch('F', resultIn, fieldName, iordr, fieldIn, iret)
                call rsexch(' ', resultOut, fieldName, iordr, fieldOut, iret)
                chams0 = '&&CHRPEL.CHAMS0'
                chams1 = '&&CHRPEL.CHAMS1'
                chafus = '&&CHRPEL.CHAFUS'
                chs(1) = chams0
                chs(2) = chams1
                call celces(fieldIn, 'V', chams0)
                call celces(fieldOut, 'V', chams1)
                call cesfus(2, chs, lcumu, lcoer, lcoec, lcoc(1), 'V', chafus)
                call dismoi('NOM_LIGREL', fieldIn, 'CHAM_ELEM', repk=ligrelField)
                call cescel(chafus, ligrelField, fieldName, ' ', 'NAN', nncp, 'G', fieldIn, &
                            'F', ibid)
                call detrsd('CHAMP', fieldOut)
                call jedema()
            end do
        end do
        call detrsd('CHAMP', chams0)
        call detrsd('CHAMP', chams1)
        call detrsd('CHAMP', chafus)
        call detrsd('RESULTAT', resultOut)
    end if
!
    call jedema()
end subroutine
