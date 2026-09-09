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
subroutine thcalr(newcal, tysd, jvListStore, loadNameJv, resultIn, &
                  resultOut, nbStore, model, materField, caraElem, &
                  nbLoad)
!
    use result_module, only: rsCopyPara
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcop.h"
#include "asterfort/callCalcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/erglth.h"
#include "asterfort/exlima.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mecham.h"
#include "asterfort/medom1.h"
#include "asterfort/modopt.h"
#include "asterfort/reslgn.h"
#include "asterfort/resth2.h"
#include "asterfort/resthe.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexc1.h"
#include "asterfort/rsexc2.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8) :: nbStore, nbLoad
    character(len=8) :: resultIn, resultOut, model, caraElem
    character(len=16) :: tysd
    character(len=19) :: jvListStore, loadNameJv
    character(len=24) :: materField
    aster_logical :: newcal
!
! --------------------------------------------------------------------------------------------------
!
! IN  NEWCAL : TRUE POUR UN NOUVEAU CONCEPT RESULTAT, FALSE SINON
! IN  TYSD   : TYPE DU CONCEPT ATTACHE A RESUCO
! IN  KNUM   : NOM D'OBJET DES NUMERO D'ORDRE
! IN  KCHA   : NOM JEVEUX OU SONT STOCKEES LES CHARGES
! IN  RESUCO : NOM DE CONCEPT RESULTAT
! IN  RESUC1 : NOM DE CONCEPT DE LA COMMANDE CALC_ERREUR
! IN  CONCEP : TYPE DU CONCEPT ATTACHE A RESUC1
! IN  NBORDR : NOMBRE DE NUMEROS D'ORDRE
! IN  MODELE : NOM DU MODELE
! IN  MATE   : NOM DU CHAMP MATERIAU
! IN  CARA   : NOM DU CHAMP DES CARACTERISTIQUES ELEMENTAIRES
! IN  NCHAR  : NOMBRE DE CHARGES
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8), parameter :: numeHarm = 0
    integer(kind=8) :: vali
    integer(kind=8) :: iStore, numeStore, jcha, iret1, iret, bufin1, iad, numeStore0
    integer(kind=8) :: ifm, niv, linst, niveau, nbRet
    integer(kind=8) :: iOption, nbOption
    real(kind=8) :: valthe, insold, inst
    character(len=8) :: mesh
    character(len=8) :: psourc
    character(len=16) :: option
    character(len=19) :: cartef, nomgdf, carteh, nomgdh, cartet, nomgdt, cartes
    character(len=19) :: nomgds, jvResultOut
    character(len=24) :: chelem, chtemm, chtemp
    character(len=24) :: chflum, chsour, chflup, cherre, cherrn
    character(len=24) :: chgeom, chharm, materCode
    character(len=24), parameter :: listOptionJv = '&&THCALR.LES_OPTION'
    character(len=24) :: ligrel, modelLigrel
    aster_logical :: evol
    integer(kind=8), pointer :: listStore(:) => null()
    character(len=16), pointer :: listOption(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call jerecu('V')
    call infmaj()
    call infniv(ifm, niv)

! - Initializations
    chgeom = " "
    chtemp = " "
    chharm = " "
    chelem = " "
    jvResultOut = resultOut

! - Create list of options to compute
    call getvtx(' ', 'OPTION', nbval=0, nbret=nbRet)
    nbOption = -nbRet
    call wkvect(listOptionJv, 'V V K16', nbOption, vk16=listOption)
    call getvtx(' ', 'OPTION', nbval=nbOption, vect=listOption, nbret=nbRet)
    call modopt(resultIn, model, listOptionJv, nbOption)
    call jeveuo(listOptionJv, 'L', vk16=listOption)

! - Access to loads
    call jeveuo(loadNameJv//'.LCHA', 'L', jcha)

! - Access to storage
    call jeveuo(jvListStore, 'L', vi=listStore)
    numeStore0 = listStore(1)

! - Create new datastructure
    if (newcal) then
        call rscrsd('G', resultOut, tysd, nbStore)
        call titre()
    end if

! - Copy parameters
    if (newcal) then
        call rsCopyPara(resultIn, jvResultOut, nbStore, listStore)
    end if
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call exlima(' ', 0, 'V', model, ligrel)

! - Process options
    do iOption = 1, nbOption
        option = listOption(iOption)
        if (option .eq. ' ') goto 120

        if (callCalcul(option)) then
            call calcop(option, listOptionJv, resultIn, resultOut, jvListStore, &
                        nbStore, tysd, iret)
            if (iret .eq. 0) goto 120
        end if

! ----- Get parameters
        call medom1(model, materField, materCode, caraElem, loadNameJv, nbLoad, &
                    resultIn, numeStore0)
        call jeveuo(loadNameJv//'.LCHA', 'L', jcha)
!
        call mecham(option, model, numeHarm, &
                    chgeom, chharm, iret)
        if (iret .ne. 0) goto 190
!
!    ------------------------------------------------------------------
!    -- OPTION "ERTH_ELEM"
!    ------------------------------------------------------------------
!
        if (option .eq. 'ERTH_ELEM') then
!
!
! PAR DECRET EDA DU 22/08/01 ON SUPPRIME LE PARAMETRE NIVEAU ET ON LE
! FIXE A 2 (15 VALEURS DE PARAMETRES).
            niveau = 2
!
! RECUPERATION NIVEAU AFFICHAGE
            call infniv(ifm, niv)
!
! BOUCLE SUR LES INSTANTS CHOISIS PAR LE USER
            insold = zero
            chtemm = ' '
            chtemp = ' '
            chflum = ' '
            chflup = ' '
!
! PREPARATION DES CALCULS D'INDICATEUR (CONNECTIVITE INVERSE, CHARGE)
            call jeveuo(loadNameJv//'.LCHA', 'L', jcha)
            call resth2(model, modelLigrel, zk8(jcha), nbLoad, mesh, &
                        cartef, nomgdf, carteh, nomgdh, cartet, &
                        nomgdt, cartes, nomgds, chgeom, chsour, &
                        psourc)
!
            if (niv .ge. 1) then
                write (ifm, *)
                write (ifm, *) '*********************************************'
                write (ifm, *) '  CALCUL DE CARTES D''ERREURS EN RESIDU'
                write (ifm, *) '       POUR LE PROBLEME THERMIQUE'
                write (ifm, *)
                write (ifm, *) '  OPTION DE CALCUL   ERTH_ELEM'
                write (ifm, *) '  MODELE                ', model
                write (ifm, *) '  SD EVOL_THER DONNEE   ', resultIn
                write (ifm, *) '             RESULTAT   ', resultOut
                write (ifm, *)
                write (ifm, *) '* CONTRAIREMENT AUX CALCULS THERMIQUES, POUR *'
                write (ifm, *) '* UN TYPE DE CHARGEMENT DONNE, ON NE RETIENT *'
                write (ifm, *) '* QUE LA DERNIERE OCCURENCE DE AFFE_CHAR_THER*'
                write (ifm, *) '  LISTE DES CHARGEMENTS :'
                do bufin1 = 1, nbLoad
                    write (ifm, *) '                        ', zk8(jcha+bufin1-1)
                end do
                write (ifm, *) '  CL DE FLUX RETENUE      ', nomgdf
                write (ifm, *) '  CL D''ECHANGE RETENUE    ', nomgdh
                write (ifm, *) '  SOURCE RETENUE          ', nomgds
                write (ifm, *) '  MATERIAU PRIS EN COMPTE ', materField(1:8)
                write (ifm, *) '  NOMBRE DE NUMERO D''ORDRE ', nbStore
            end if
!
! BOUCLE SUR LES PAS DE TEMPS
            do iStore = 1, nbStore
                call jemarq()
                call jerecu('V')
                numeStore = listStore(iStore)
                call medom1(model, materField, materCode, caraElem, loadNameJv, nbLoad, &
                            resultIn, numeStore)
! RECUPERATION DU PARM_THETA CORRESPONDANT A IORDR
                call jenonu(jexnom(resultIn//'           .NOVA', 'PARM_THETA'), iad)
                if (iad .eq. 0) then
                    valthe = 0.57d0
                    call utmess('A', 'CALCULEL4_98', sk=resultIn)
                else
                    call rsadpa(resultIn, 'L', 1, 'PARM_THETA', numeStore, 0, sjv=iad)
                    valthe = zr(iad)
                    if ((valthe .gt. 1.d0) .or. (valthe .lt. 0.d0)) then
                        call utmess('F', 'INDICATEUR_5', sk=resultIn)
                    end if
                end if
                if (niv .ge. 1) then
                    write (ifm, *) '   PARAM-THETA/IORDR ', valthe, numeStore
                    if (iStore .eq. nbStore) then
                        write (ifm, *) '**********************************************'
                        write (ifm, *)
                    end if
                end if
!
! CALCUL DU CRITERE D'EVOLUTION LEVOL (TRUE=TRANSITOIRE)
! CAS PARTICULIER DE L'INSTANT INITIAL D'UN CALCUL TRANSITOIRE
! ON ESTIME SON ERREUR COMME EN STATIONNAIRE
                if (iStore .eq. 1) then
                    evol = .false.
                else
                    evol = .true.
                end if
!
! RECUPERATION DU NOM DES CHAMP_GD = RESUCO('FLUX_ELNO',I)
! ET RESUCO('TEMP',I) POUR I=IORDR. POUR IORDR-1 ILS SONT STOCKES
! DANS CHFLUM/CHTEMM DEPUIS LA DERNIERE ITERATION.
! RESUCO = NOM USER DE LA SD DESIGNEE PAR LE MOT-CLE RESULTAT
                call rsexc2(1, 1, resultIn, 'TEMP', numeStore, &
                            chtemp, option, iret)
                if (iret .gt. 0) then
                    vali = numeStore
                    call utmess('F', 'CALCULEL6_46', si=vali)
                end if
                call rsexc2(1, 1, resultIn, 'FLUX_ELNO', numeStore, &
                            chflup, option, iret)
                if (iret .gt. 0) then
                    vali = numeStore
                    call utmess('F', 'CALCULEL6_47', si=vali)
                end if
!
! RECUPERATION DE L'INSTANT CORRESPONDANT A IORDR
                call rsadpa(resultIn, 'L', 1, 'INST', numeStore, 0, sjv=linst)
                inst = zr(linst)
!
! IMPRESSIONS NIVEAU 2 POUR DIAGNOSTIC...
                if (niv .eq. 2) then
                    write (ifm, *) 'THCALR **********'
                    write (ifm, *) 'EVOL/I/IORDR', evol, iStore, numeStore
                    write (ifm, *) 'INST/INSOLD', inst, insold
                    write (ifm, *) 'CHTEMM/CHTEMP', chtemm, ' / ', chtemp
                    write (ifm, *) 'CHFLUM/CHFLUP', chflum, ' / ', chflup
                end if
!
! RECUPERATION DU NOM DU CHAMP_GD = RESUC1('ERTH_ELEM',IORDR)
! RESUC1 = NOM USER DE LA SD CORRESPONDANT AU RESULTAT DE CALC_ERREUR
                call rsexc1(jvResultOut, option, numeStore, chelem)
! PREPARATION DES DONNEES/LANCEMENT DU CALCUL DES INDICATEURS
                call resthe(modelLigrel, evol, chtemm, chtemp, chflum, &
                            chflup, materCode, valthe, insold, inst, &
                            chelem, niveau, ifm, niv, mesh, &
                            cartef, nomgdf, carteh, nomgdh, cartet, &
                            nomgdt, cartes, nomgds, chgeom, chsour, &
                            psourc, iStore)
! CALCUL DE L'ESTIMATEUR GLOBAL
                call erglth(chelem, inst, niveau, numeStore, resultIn)
! NOTATION DE LA SD RESULTAT LERES1
                call rsnoch(jvResultOut, option, numeStore)
!
! INIT. POUR LE NUMERO D'ORDRE SUIVANT
                if (nbStore .ne. 1 .and. iStore .ne. nbStore) then
                    chtemm = chtemp
                    chflum = chflup
                    insold = inst
                end if
                call jedema()
            end do
! DESTRUCTION DES OBJETS JEVEUX VOLATILES
            call jedetr(cartef//'.PTMA')
            call jedetr(carteh//'.PTMA')
            call jedetr(cartet//'.PTMA')
            call jedetr(cartef//'.PTMS')
            call jedetr(carteh//'.PTMS')
            call jedetr(cartet//'.PTMS')
!
!    ------------------------------------------------------------------
!    -- OPTION "ERTH_ELNO"
!    ------------------------------------------------------------------
!
        else if (option .eq. 'ERTH_ELNO') then
!
            do iStore = 1, nbStore
                call jemarq()
                call jerecu('V')
                numeStore = listStore(iStore)
! RECUPERATION DU NOM DU CHAMP_GD = RESUCO('ERTH_ELEM',IORDR)
                call rsexc2(1, 1, resultIn, 'ERTH_ELEM', numeStore, &
                            cherre, option, iret1)
                if (iret1 .gt. 0) goto 40
! RECUPERATION DU NOM DU CHAMP_GD = RESUC1('ERTH_ELNO',IORDR)
! RESUC1 = NOM USER DE LA SD CORRESPONDANT AU RESULTAT DE CALC_ERREUR
                call rsexc1(jvResultOut, option, numeStore, cherrn)
                call reslgn(modelLigrel, option, cherre, cherrn)
! NOTATION DE LA SD RESULTAT LERES1
                call rsnoch(jvResultOut, option, numeStore)
40              continue
                call jedema()
            end do
        else
            call utmess('A', 'CALCULEL3_22', sk=option)
        end if
!
120     continue
    end do
!
190 continue
!
    call jedema()
end subroutine
