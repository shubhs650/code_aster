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
subroutine mecalr(newcal, tysd, jvListStore, loadNameJv, resultIn, &
                  resultOut, nbStore, model, materField, caraElem, &
                  nbLoad)
!
    use result_module, only: rsCopyPara
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcop.h"
#include "asterfort/callCalcul.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesces.h"
#include "asterfort/cetule.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
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
#include "asterfort/meca01.h"
#include "asterfort/mecham.h"
#include "asterfort/medom1.h"
#include "asterfort/modopt.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexc1.h"
#include "asterfort/rsexc2.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/singue.h"
#include "asterfort/singum.h"
#include "asterfort/sinoz1.h"
#include "asterfort/sinoz2.h"
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
! IN  KNUM   : NOM D'OBJET DES NUMEROS D'ORDRE
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
    integer(kind=8), parameter :: numeHarm = 0
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: numeStore0, numeStore
    integer(kind=8) :: iret, iret1, iret2, iret3, iret4, iret5, ireter, nbRet
    integer(kind=8) :: nbOption
    integer(kind=8) :: iStore, ibid, iOption
    integer(kind=8) :: jcha
    integer(kind=8) :: jcoor, ltymo
    integer(kind=8) :: nnoem, nelem, ndim, nncp
    character(len=8) :: mesh
    character(len=19) :: pfchno
    character(len=16) :: option, types
    character(len=19) :: jvResultOut
    character(len=19) :: cherrs, chenes, chsins, chsinn
    character(len=24) :: cheneg, chsing, cherr1, cherr2, cherr3, cherr4, materCode
    character(len=24) :: chamgd, chsig, chsign, chgeom, chharm, chelem
    character(len=24) :: ligrel
    character(len=24), parameter :: listOptionJv = '&&MECALR.LES_OPTION'
    character(len=24) :: modelLigrel
    character(len=19) :: chvarc
!
    real(kind=8) :: prec
    real(kind=8) :: tbgrca(3)
!
    character(len=24) :: valkm(2)
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: meshDime(:) => null()
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
    chamgd = " "
    chgeom = " "
    chharm = " "
    chsig = " "
    chelem = " "
    chvarc = '&&MECALR.CHVARC'
    jvResultOut = resultOut

! - Create list of options to compute
    call getvtx(' ', 'OPTION', nbval=0, nbret=nbRet)
    nbOption = -nbRet
    call wkvect(listOptionJv, 'V V K16', nbOption, vk16=listOption)
    call getvtx(' ', 'OPTION', nbval=nbOption, vect=listOption, nbret=nbRet)
    call modopt(resultIn, model, listOptionJv, nbOption)
    call jeveuo(listOptionJv, 'L', vk16=listOption)

!     ON RECUPERE LE TYPE DE MODE: DYNAMIQUE OU STATIQUE
    if (tysd .eq. 'MODE_MECA') then
        call rsadpa(resultIn, 'L', 1, 'TYPE_MODE', 1, 0, sjv=ltymo)
    end if

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
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call exlima(' ', 0, 'V', model, ligrel)

! - GRANDEURS CARACTERISTIQUES DE L'ETUDE
    call cetule(model, tbgrca, iret)

! - Process options
    do iOption = 1, nbOption
        option = listOption(iOption)
        if (option .eq. ' ') goto 660

        if (callCalcul(option)) then
            call calcop(option, listOptionJv, resultIn, resultOut, jvListStore, &
                        nbStore, tysd, iret)
            if (iret .eq. 0) goto 660
        end if

! ----- Get parameters
        call medom1(model, materField, materCode, caraElem, loadNameJv, nbLoad, &
                    resultIn, numeStore0)
        call jeveuo(loadNameJv//'.LCHA', 'L', jcha)
!
        call mecham(option, model, numeHarm, &
                    chgeom, chharm, iret)
        if (iret .ne. 0) goto 690
!
!    ------------------------------------------------------------------
!    -- OPTIONS "SIZ1_NOEU","SIZ2_NOEU"
!    ------------------------------------------------------------------
        if (option .eq. 'SIZ1_NOEU' .or. option .eq. 'SIZ2_NOEU') then
!
!
            do iStore = 1, nbStore
                call jemarq()
                call jerecu('V')
                numeStore = listStore(iStore)

                call medom1(model, materField, materCode, caraElem, loadNameJv, nbLoad, &
                            resultIn, numeStore)
                call jeveuo(loadNameJv//'.LCHA', 'L', jcha)
                call rsexc2(1, 1, resultIn, 'DEPL', numeStore, &
                            chamgd, option, iret)
                if (iret .gt. 0) goto 150
                call rsexc2(1, 1, resultIn, 'SIEF_ELGA', numeStore, &
                            chsig, option, iret)
                if (iret .gt. 0) then
                    call utmess('A', 'CALCULEL3_7', sk=option)
                    call jedema()
                    goto 660
!
                end if
                call rsexc1(jvResultOut, option, numeStore, chsign)
                if (option .eq. 'SIZ1_NOEU') then
                    call sinoz1(model, chsig, chsign)
                else if (option .eq. 'SIZ2_NOEU') then
                    call dismoi('NUME_EQUA', chamgd, 'CHAM_NO', repk=pfchno)
                    call sinoz2(model, pfchno, chsig, chsign)
                end if
                call rsnoch(jvResultOut, option, numeStore)
150             continue
                call jedema()
            end do
!
!    ------------------------------------------------------------------
!    -- OPTIONS DES INDICATEURS D'ERREURS
!    ------------------------------------------------------------------
        elseif (option .eq. 'ERZ1_ELEM' .or. option .eq. 'ERZ2_ELEM' .or. &
                option .eq. 'ERME_ELEM' .or. option .eq. 'ERME_ELNO' .or. &
                option .eq. 'QIRE_ELEM' .or. option .eq. 'QIRE_ELNO' .or. &
                option .eq. 'QIZ1_ELEM' .or. option .eq. 'QIZ2_ELEM') then
!
            call meca01(option, nbStore, listStore, nbLoad, jcha, &
                        loadNameJv, tbgrca, resultIn, resultOut, &
                        jvResultOut, mesh, model, modelLigrel, materField, &
                        caraElem, chvarc, iret)
!
            if (iret .eq. 1) then
                goto 690
!
            else if (iret .eq. 2) then
                goto 660
!
            end if
!
!    ------------------------------------------------------------------
!    -- OPTION "SING_ELEM"
!    ------------------------------------------------------------------
        else if (option .eq. 'SING_ELEM') then
!
            call getvr8(' ', 'PREC_ERR', scal=prec, nbret=iret1)
            if (iret1 .ne. 1) then
                call utmess('F', 'CALCULEL3_12')
            else
                if (prec .le. 0.d0) then
                    call utmess('F', 'CALCULEL3_13')
                end if
            end if
!
            types = ' '
            call getvtx(' ', 'TYPE_ESTI', scal=types, nbret=ireter)
            if (ireter .gt. 0) then
                call utmess('I', 'CALCULEL3_24', sk=types)
            end if
!
! 1 - RECUPERATION DE :
!  NNOEM : NOMBRE DE NOEUDS
!  NELEM : NOMBRE D ELEMENTS FINIS (EF)
!  NDIM  : DIMENSION
!  JCOOR : ADRESSE DES COORDONNEES
!  JTYPE : ADRESSE DU TYPE D ELEMENTS FINIS
!
            call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
!
            call jeveuo(mesh//'.DIME', 'L', vi=meshDime)
            call jeveuo(mesh//'.COORDO    .VALE', 'L', jcoor)
            call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
!
            nnoem = meshDime(1)
            nelem = meshDime(3)
            ndim = meshDime(6)
!
! 2 - CREATION D OBJETS TEMPORAIRES UTILES POUR LA SUITE
! '&&SINGUM.DIME' (DIM=3) CONTIENT
!   NBRE MAX DE NOEUDS SOMMETS CONNECTES AUX EF (NSOMMX)
!   NBRE MAX D EF CONNECTES AUX NOEUDS (NELCOM)
!   DEGRE DES EF (1 SI LINEAIRE ET 2 SI QUADRATIQUE)
! '&&SINGUM.MESU' (DIM=NELEM) CONTIENT L AIRE OU LE VOLUME DES EFS
! '&&SINGUM.CONN' (DIM=NELEM*(NSOMMX+2)) CONTIENT
!   1ERE VALEUR = NBRE DE NOEUDS SOMMETS CONNECTES A L EF N
!   2EME VALEUR = 1 SI EF EST SURFACIQUE EN 2D ET VOLUMIQUE EN 3D
!                 0 SINON
!   CONNECTIVITE EF N
! '&&SINGUM.CINV' (DIM=NNOEM*(NELCOM+2)) CONTIENT
!   1ERE VALEUR = NBRE D EF CONNECTES AU NOEUD N
!   2EME VALEUR = 0 NOEUD MILIEU OU NON CONNECTE A UN EF UTILE
!                 1 NOEUD SOMMET A L INTERIEUR + LIE A UN EF UTILE
!                 2 NOEUD SOMMET BORD + LIE A UN EF UTILE
!                 EF UTILE = EF SURF EN 2D ET VOL EN 3D
!   CONNECTIVITE INVERSE NOEUD N
!
            call singum(mesh, ndim, nnoem, nelem, typmail, &
                        zr(jcoor))
!
! 3 - BOUCLE SUR LES INSTANTS DEMANDES
!
            do iStore = 1, nbStore
                call jemarq()
                numeStore = listStore(iStore)
!
                if (ireter .gt. 0) then
                    call rsexch(' ', resultIn, types, numeStore, cherr4, iret5)
!
                    if (iret5 .gt. 0) then
                        valkm(1) = types
                        valkm(2) = resultIn
                        call utmess('A', 'CALCULEL3_26', nk=2, valk=valkm)
                        iret = 1
                    end if
!
! 3.1 - RECUPERATION DE LA CARTE D ERREUR ET D ENERGIE
!       SI PLUSIEURS INDICATEURS ON PREND PAR DEFAUT
!       ERME_ELEM SI IL EST PRESENT
!       ERZ2_ELEM PAR RAPPORT A ERZ1_ELEM
!
                else
!
                    iret5 = 1
                    call rsexch(' ', resultIn, 'ERME_ELEM', numeStore, cherr1, &
                                iret1)
                    call rsexch(' ', resultIn, 'ERZ1_ELEM', numeStore, cherr2, &
                                iret2)
                    call rsexch(' ', resultIn, 'ERZ2_ELEM', numeStore, cherr3, &
                                iret3)
!
                    if (iret1 .gt. 0 .and. iret2 .gt. 0 .and. iret3 .gt. 0) then
                        call utmess('A', 'CALCULEL3_14')
                        iret = 1
                    end if
!
                end if
!
                if (tysd .eq. 'EVOL_NOLI') then
                    call rsexch(' ', resultIn, 'ETOT_ELEM', numeStore, cheneg, &
                                iret4)
                else
                    call rsexch(' ', resultIn, 'EPOT_ELEM', numeStore, cheneg, &
                                iret4)
                end if
                if (iret4 .gt. 0) then
                    call utmess('A', 'CALCULEL3_29')
                end if
!
                if ((iret+iret4) .gt. 0) then
                    call utmess('A', 'CALCULEL3_36')
                    goto 250
!
                end if
! 3.2 - TRANSFORMATION DE CES DEUX CARTES EN CHAM_ELEM_S
!
                cherrs = '&&MECALR.ERRE'
!
                if (iret5 .eq. 0) then
                    call celces(cherr4(1:19), 'V', cherrs)
                else if (iret1 .eq. 0) then
                    call celces(cherr1(1:19), 'V', cherrs)
                    if ((iret2 .eq. 0) .or. (iret3 .eq. 0)) then
                        call utmess('A', 'CALCULEL3_15')
                    end if
                else if (iret3 .eq. 0) then
                    call celces(cherr3(1:19), 'V', cherrs)
                    if (iret2 .eq. 0) then
                        call utmess('A', 'CALCULEL3_16')
                    end if
                else if (iret2 .eq. 0) then
                    call celces(cherr2(1:19), 'V', cherrs)
                else
                    ASSERT(.false.)
                end if
!
                chenes = '&&MECALR.ENER'
                call celces(cheneg(1:19), 'V', chenes)
!
! 3.3 - ROUTINE PRINCIPALE QUI CALCULE DANS CHAQUE EF :
!       * LE DEGRE DE LA SINGULARITE
!       * LE RAPPORT ENTRE L ANCIENNE ET LA NOUVELLE TAILLE
!       DE L EF CONSIDERE
!       => CE RESULAT EST STOCKE DANS CHELEM (CHAM_ELEM)
!       CES DEUX COMPOSANTES SONT CONSTANTES PAR ELEMENT
!
                call rsexc1(jvResultOut, option, numeStore, chelem)
!
                call singue(cherrs, chenes, mesh, ndim, nnoem, &
                            nelem, zr(jcoor), prec, modelLigrel, chelem, &
                            types)
!
                call rsnoch(jvResultOut, option, numeStore)
!
! 3.4 - DESTRUCTION DES CHAM_ELEM_S
!
                call detrsd('CHAM_ELEM_S', cherrs)
                call detrsd('CHAM_ELEM_S', chenes)
!
250             continue
                call jedema()
            end do
!
! 4 - DESTRUCTION DES OBJETS TEMPORAIRES
!
            call jedetr('&&SINGUM.DIME           ')
            call jedetr('&&SINGUM.MESU           ')
            call jedetr('&&SINGUM.CONN           ')
            call jedetr('&&SINGUM.CINV           ')
!    ------------------------------------------------------------------
!    -- OPTION "SING_ELNO"
!    ------------------------------------------------------------------
        else if (option .eq. 'SING_ELNO') then
            do iStore = 1, nbStore
                call jemarq()
                numeStore = listStore(iStore)
!
! 1 - RECUPERATION DE LA CARTE DE SINGULARITE
!
                call rsexc2(1, 1, resultIn, 'SING_ELEM', numeStore, &
                            chsing, option, iret1)
!
                if (iret1 .gt. 0) goto 270
!
! 2 - TRANSFORMATION DE CE CHAMP EN CHAM_ELEM_S
!
                chsins = '&&MECALR.SING'
                call celces(chsing(1:19), 'V', chsins)
!
! 3 - TRANSFOMATION DU CHAMP CHSINS ELEM EN ELNO
!
                chsinn = '&&MECALR.SINN'
                call cesces(chsins, 'ELNO', ' ', ' ', ' ', &
                            'V', chsinn)
!
! 4 - STOCKAGE
!
                call rsexc1(jvResultOut, option, numeStore, chelem)
!
                call cescel(chsinn, modelLigrel(1:19), 'SING_ELNO', 'PSINGNO', 'NON', &
                            nncp, 'G', chelem(1:19), 'F', ibid)
!
                call rsnoch(jvResultOut, option, numeStore)
!
! 5 - DESTRUCTION DES CHAM_ELEM_S
!
                call detrsd('CHAM_ELEM_S', chsins)
                call detrsd('CHAM_ELEM_S', chsinn)
!
270             continue
                call jedema()
            end do
!
!      -----------------------------------------------------------------
!
        else
            call utmess('A', 'CALCULEL3_22', sk=option)
        end if
!
660     continue
    end do
!
!============= FIN DE LA BOUCLE SUR LES OPTIONS A CALCULER =============
!
690 continue
    call jedema()
end subroutine
