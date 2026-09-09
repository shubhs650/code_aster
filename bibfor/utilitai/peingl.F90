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
subroutine peingl(tablOutZ, &
                  modelZ, materFieldZ, materCodeZ, caraElemZ, numeHarm, &
                  nbFactorKeyword, factorKeywordZ)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/alchml.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/exisdg.h"
#include "asterfort/exixfe.h"
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
#include "asterfort/jexnum.h"
#include "asterfort/mecham.h"
#include "asterfort/mesomm.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/setStructFields.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/umalma.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/wkvect.h"
#include "asterfort/xajcin.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: tablOutZ
    character(len=*), intent(in) :: modelZ, materFieldZ, materCodeZ, caraElemZ
    integer(kind=8), intent(in) :: numeHarm, nbFactorKeyword
    character(len=*), intent(in) :: factorKeywordZ

    ! character(len=*) :: modele
    ! character(len=*) :: mate, mateco
    ! character(len=*) :: cara
    ! integer(kind=8) :: nh
    ! integer(kind=8) :: nbocc
    ! character(len=*) :: factorKeywordZ
!
! --------------------------------------------------------------------------------------------------
!
!     OPERATEUR   POST_ELEM
!                  TRAITEMENT DU MOT-FACTEUR "INDIC_ENER"
!                          ET DU MOT-FACTEUR "INDIC_SEUIL"
!
!                  CALCUL DES INDICATEURS GLOBAUX DE
!                  DE PERTE DE PROPORTIONNALITE DU CHARGEMENT.
!
! --------------------------------------------------------------------------------------------------
!
!           -POUR LE MOT-CLE INDIC_ENER, ON CALCULE L'INDICATEUR
!            GLOBAL ENERGETIQUE DETERMINE PAR L'EXPRESSION
!            SUIVANTE :
!            IE = (SOMME_DOMAINE((1 - PSI(EPS)/OMEGA(EPS,VARI)).DV)/V
!
!        OU  .OMEGA EST LA DENSITE D'ENERGIE TOTALE
!            (I.E. OMEGA = SOMME_0->T(SIGMA:D(EPS)/DT).DTAU
!            .PSI EST LA DENSITE D'ENERGIE ELASTIQUE 'TOTALE'
!            (I.E. ASSOCIEE A LA COURBE DE TRACTION SI ON
!                  CONSIDERAIT LE MATERIAU ELASTIQUE NON-LINEAIRE)
!            .V EST LE VOLUME DU GROUPE DE MAILLES TRAITE
!
! -----------------------------------------------------------------
!
!           -POUR LE MOT-CLE INDIC_SEUIL, ON CALCULE L'INDICATEUR
!            GLOBAL  DETERMINE PAR L'EXPRESSION SUIVANTE :
!
!   IS = (SOMME_DOMAINE(1 - ((SIG-X):EPS_PLAST)/((SIG_Y+R)*P)).DV)/V
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES
!            .X         EST LE TENSEUR DE RAPPEL
!            .EPS_PLAST EST LE TENSEUR DES DEFORMATIONS PLASTIQUES
!            .SIG_Y     EST LA LIMITE D'ELASTICITE
!            .R         EST LA FONCTION D'ECROUISSAGE
!            .P         EST LA DEFORMATION PLASTIQUE CUMULEE
!            .V EST LE VOLUME DU GROUPE DE MAILLES TRAITE
! -----------------------------------------------------------------
!
!  MOT-CLE ENER_ELAS : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  =================   DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!   ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES
!            .D         EST LE TENSEUR DE HOOKE
!
! -----------------------------------------------------------------
!
!  MOT-CLE ENER_ELTR : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  =================   MODIFIEE DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!   ENELAS =  0.5*Lame*H(tr(EPS))*tr(EPS)**2+mu*SUM(H(Ei)*Ei**2)
!
!        OU  .EPS      EST LE TENSEUR DES DEFORMATIONS ELASTIQUES
!            .Ei       SONT (pour i=1..3) LES DEFORMATIONS PROPRES
!            .H        LA FONCTION D'HEAVISIDE
!
! -----------------------------------------------------------------
!
!  MOT-CLE ENER_TOTALE : CALCUL DE L'ENERGIE DE DEFORMATION TOTALE
!  ===================   DETERMINEE PAR L'EXPRESSION S        character(len=*) :: resu

!
!   ENER_TOTALE =  ENELAS + EPLAS
!
!          AVEC : ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!                 ENELAS EST L'ENERGIE DE DEFORMATION ELASTIQUE
!
!           OU  .SIG       EST LE TENSEUR DES CONTRAINTES
!       !
! -----------------------------------------------------------------
!        .D         EST LE TENSEUR DE HOOKE
!
!          ET   : EPLAS = SOMME_VOLUME((R(P))*D(P))
!                 EPLAS EST L'ENERGIE DE DEFORMATION PLASTIQUE
!
!           OU  .P         EST LA DEFORMATION PLASTIQUE CUMULEE
!           ET   R(P) EST CALCULE POUR LES COMPORTEMENTS SUIVANTS :
!                      .VMIS_ISOT_LINE
!                      .VMIS_ISOT_TRAC
!                      .VMIS_ECMI_LINE
!                      .VMIS_ECMI_TRAC
!                      .VMIS_CINE_LINE
!                      .VISC_CIN1_CHAB
!                      .VISC_CIN2_CHAB
!
!          POUR LES AUTRES COMPORTEMENTS ON S'ARRETE EN ERREUR FATALE
!
! -----------------------------------------------------------------
!
!  MOT-CLE ENER_DISS : CALCUL DE L'ENERGIE DISSIPEE
!  =================   DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!   EDISS =  SOMME_VOLUME(VINT*K0)
!
!        OU  .VINT      EST LA VARIABLE INTERNE ASSOCIEE A LA DISSIPATION
!            .K0        EST LE SEUIL DU DOMAINE ELASTIQUE
!
! -----------------------------------------------------------------
!
!   ARGUMENT        E/S  TYPE         ROLE
!    RESU           VAR    K*      TABLE EN SORTIE DE LA COMMANDE
!    MODELE         IN     K*      NOM DU MODELE SUR-LEQUEL ON FAIT
!                                  LE CALCUL
!    MATE           IN     K*      NOM DU CHAMP MATERIAU
!    CARA           IN     K*      NOM DU CHAMP DES CARA_ELEM
!    NH             IN     I       NUMERO D'HARMONIQUE DE FOURIER
!    NBOCC          IN     I       NOMBRE D'OCCURENCES DU MOT-FACTEUR
!                                  INDIC_ENER
!    MOTFAZ         IN     K*      NOM DU MOT-CLE FACTEUR "INDIC_ENER"
!                                                     OU  "INDIC_SEUIL"
!                                                     OU  "ENER_ELAS"
!                                                     OU  "ENER_ELTR"
!                                                     OU  "ENER_TOTALE"
!                                                     OU  "ENER_DISS"
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    integer(kind=8) :: nbFieldIn, nbFieldOut
    integer(kind=8) :: nbStore, iStore, numeStore, numeStorePrev
    integer(kind=8) :: iret, jvPara, nt, nbRet
    integer(kind=8) :: ng, nbgrma, jgr, ig, nbCell, jad, nbCellil, jma, im, iFactorKeyword, nume
    integer(kind=8) :: ngdmax, ncmpmx, igd, idebgd, dg, iCell, iconex, nbno, nec, ivari, nm
    integer(kind=8) :: nbgrma_tot, deca, nbtot, nbMaiT
    real(kind=8) :: work(5), indic1, volume, inst, valr(6), prec
    real(kind=8) :: energy_tout, energy_ma
    complex(kind=8) :: c16b
    character(len=19), parameter :: chvarc = '&&PEECIN.VARC', chvref = '&&PEINGL.CHVARC.REF'
    character(len=19), parameter :: listStoreJv = '&&PEINGL.NUME_ORDRE'
    integer(kind=8), pointer :: listStore(:) => null()
    character(len=19), parameter :: listTimeJv = '&&PEECIN.INSTANT'
    real(kind=8), pointer :: listTime(:) => null()
    character(len=2) :: codret
    character(len=8) :: result, crit, mesh, nommai, vk8(2), numeStorePrevStr
    character(len=8) :: numeStoreStr, k8b, physQuanName
    character(len=16) :: resultType, factorKeyword, modelLigrel, compt, option
    character(len=19) :: ligrel, compor
    character(len=24) :: chgeom, caraElem, chharm, chvari, chdepl
    character(len=24) :: vk24(2), nomgrm, chsig
    character(len=24) :: chsigm, chdepm, chbid
    aster_logical :: evol
    integer(kind=8), pointer :: ptma(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    character(len=16), pointer :: vale(:) => null()
    real(kind=8), pointer :: energy_grpma(:) => null()
    integer(kind=8), pointer :: v_allma(:) => null()
    aster_logical :: lxfem
    integer(kind=8) :: nbParaResu
    integer(kind=8), parameter :: nbParaResuMax = 9
    character(len=16) :: tablParaName(nbParaResuMax)
    character(len=8), parameter :: tablParaType(nbParaResuMax) = &
                                   (/'I  ', 'R  ', 'K24', &
                                     'K8 ', 'R  ', 'R  ', &
                                     'R  ', 'R  ', 'R  '/)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    c16b = (0.d0, 0.d0)
    ivari = 0
    chbid = '&&PEINGL.VARINUL'
    compor = '&&PEINGL.COMPNUL'
    compt = 'XXXXXXXXXXXXXXXX'
    evol = .false.
    factorKeyword = factorKeywordZ
    option = factorKeywordZ
    caraElem = caraElemZ
    work = zero
    valr = zero
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)
    call exixfe(modelZ, iret)
    lxfem = iret .ne. 0
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Parameters in output table and option to compute
    tablParaName(1) = 'NUME_ORDRE'
    tablParaName(2) = 'INST'
    tablParaName(3) = 'LIEU'
    tablParaName(4) = 'ENTITE'
    tablParaName(5) = factorKeyword
    nbParaResu = 5
    if (factorKeyword(1:4) .eq. 'ENER') then
        tablParaName(5) = 'TOTALE'
        if (factorKeyword .eq. 'ENER_DISS') then
            option = 'DISS_ELEM'
        else if (factorKeyword .eq. 'ENER_ELAS') then
            option = 'ENEL_ELEM'
            tablParaName(6) = 'MEMBRANE'
            tablParaName(7) = 'FLEXION'
            tablParaName(8) = 'CISAILLE'
            tablParaName(9) = 'COUPL_MF'
            nbParaResu = 9
        else if (factorKeyword .eq. 'ENER_ELTR') then
            option = 'ENTR_ELEM'
        else if (factorKeyword .eq. 'ENER_TOTALE') then
            option = 'ENER_TOTALE'
        end if
    end if

!
    nbgrma_tot = 1
    if (factorKeyword(1:5) == 'ENER_') then
        call dismoi('NOM_MAILLA', modelZ, 'MODELE', repk=mesh)
        do iFactorKeyword = 1, nbFactorKeyword
            call getvem(mesh, 'GROUP_MA', factorKeyword, 'GROUP_MA', iFactorKeyword, 0, k8b, ng)
            if (ng < 0) then
                nbgrma_tot = nbgrma_tot-ng
            end if
        end do
    else
        nbgrma_tot = nbgrma_tot+1
    end if
    AS_ALLOCATE(vr=energy_grpma, size=nbgrma_tot)
    energy_tout = 0.d0
    energy_ma = 0.0
    energy_grpma = 0.d0

! - RECUPERATION DU RESULTAT A TRAITER
    call getvid(' ', 'RESULTAT', scal=result, nbret=nbret)
    if (nbret .eq. 0) then
        call utmess('F', 'UTILITAI3_76')
    end if
    call gettco(result, resultType)
    evol = (resultType(1:9) .eq. 'EVOL_NOLI') .or. (resultType(1:9) .eq. 'EVOL_ELAS')
    if (.not. evol) then
        call utmess('F', 'UTILITAI3_77')
    end if

! - Get parameters to create list of time steps
    call getvr8(' ', 'PRECISION', scal=prec, nbret=nbRet)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nbRet)

! - Create list of store index
    call rsutnu(result, ' ', 0, listStoreJv, nbStore, prec, crit, iret)
    if (iret .ne. 0) goto 70
    call jeveuo(listStoreJv, 'L', vi=listStore)

! - Create list of time steps
    call wkvect(listTimeJv, 'V V R', nbStore, vr=listTime)
    call jenonu(jexnom(result//'           .NOVA', 'INST'), iret)
    if (iret .ne. 0) then
        do iStore = 1, nbStore
            numeStore = listStore(iStore)
            call rsadpa(result, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
            listTime(iStore) = zr(jvPara)
        end do
    end if

! - Prepare input fields
    call mecham(option, modelZ, numeHarm, &
                chgeom, chharm, iret)
    if (iret .ne. 0) goto 80
    mesh = chgeom(1:8)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCell)
!
    call exlim3(factorKeyword, 'V', modelZ, ligrel)

! - Create output table
    call tbcrsd(tablOutZ, 'G')
    call tbajpa(tablOutZ, nbParaResu, tablParaName, tablParaType)

    do iStore = 1, nbStore
        call jemarq()
        call jerecu('V')

! ----- Current storing index
        numeStore = listStore(iStore)
        call codent(numeStore, 'G', numeStoreStr)
        inst = listTime(iStore)
        valr(1) = inst

! ----- Get external state variables
        call vrcins(modelZ, materFieldZ, caraElemZ, inst, chvarc, codret)
        call vrcref(modelZ(1:8), materFieldZ(1:8), caraElemZ(1:8), chvref(1:19))
!
        if (resultType(1:9) .eq. 'EVOL_NOLI') then
            call rsexch('F', result, 'COMPORTEMENT', numeStore, compor, iret)
            call etenca(compor, modelLigrel, iret)
            if (iret .ne. 0) then
                call utmess('F', 'UTILITAI2_62')
            end if
            call jeveuo(compor//'.DESC', 'L', vi=desc)
            ngdmax = desc(2)
            physQuanName = 'COMPOR  '
!
            call dismoi('NB_EC', physQuanName, 'GRANDEUR', repi=nec)
            if (nec .gt. 1) then
                call utmess('F', 'UTILITAI2_61')
            end if

! ---    NOMBRE DE COMPOSANTES ASSOCIEES A LA GRANDEUR  ---
            call jelira(jexnom('&CATA.GD.NOMCMP', physQuanName), 'LONMAX', ncmpmx)
!
! ---    TABLEAU DE VALEURS DE LA CARTE COMPO     ---
! ---    (CONTENANT LES VALEURS DU COMPORTEMENT)  ---
!
            call jeveuo(compor//'.VALE', 'L', vk16=vale)
!
! ---    RECUPERATION DU VECTEUR D'ADRESSAGE DANS LA CARTE  ---
! ---    CREE PAR ETENCA                                    ---
!
            call jeveuo(compor//'.PTMA', 'L', vi=ptma)
!
! ---    AFFECTATION DU TABLEAU DES NOEUDS  ---

!bFieldOut,
            do iCell = 1, nbCell
                if (ptma(iCell) .ne. 0) then
                    igd = ptma(iCell)
                    idebgd = (igd-1)*ncmpmx
                    dg = desc(1+3+2*ngdmax+ptma(iCell)-1)
!
! ---        ON S'ASSURE QUE LA PREMIERE COMPOSANTE DE LA GRANDEUR
! ---        QUI EST RELCOM A BIEN ETE AFFECTEE .or. factorKeyword(1:6) .eq. 'INDIC_') then

!
                    if (.not. exisdg([dg], 1)) then
                        call utmess('F', 'UTILITAI2_63')
                    end if
! ---        RECUPERATION DU COMPORTEMENT AFFECTE A LA MAILLE
                    compt = vale(1+idebgd+1-1)
!
! ---        RECUPERATION DES NUMEROS DES NOEUDS DE LA MAILLE
                    call jeveuo(jexnum(mesh//'.CONNEX', iCell), 'L', iconex)
                    call jelira(jexnum(mesh//'.CONNEX', iCell), 'LONMAX', nbno)
!
                end if
            end do
        end if
!
! ---  RECUPERATION DU CHAMP DE CONTRAINTES ASSOCIE AU
! ---  NUMERO D'ORDRE COURANT POUR ENER_ELAS, ENER_ELTR ET ENER_TOTALE:
!      -----------------------------------------------------
        if (factorKeyword .eq. 'ENER_TOTALE' .or. &
            factorKeyword .eq. 'ENER_ELAS' .or. &
            factorKeyword .eq. 'ENER_ELTR' .or. &
            factorKeyword(1:6) .eq. 'INDIC_') then
            call rsexch('F', result, 'SIEF_ELGA', numeStore, chsig, iret)
!
! --- SI LE NUMERO COURANT EST INFERIEUR A NBORDR ON RECUPERE LES
! --- CONTRAINTES DE L INSTANT PRECEDENT
!
            if (iStore .gt. 1) then
                numeStorePrev = listStore(iStore-1)
                call codent(numeStorePrev, 'G', numeStorePrevStr)
                call rsexch('F', result, 'SIEF_ELGA', numeStorePrev, chsigm, iret)
            end if
        end if
!
! ---  RECUPERATION DU CHAMP DES VARIABLES INTERNES ASSOCIE AU
! ---  NUMERO D'ORDRE COURANT DANS LE CAS DES EVOL_NOLI
!      ----------------------
        if (resultType(1:9) .eq. 'EVOL_NOLI') then
            call rsexch(' ', result, 'VARI_ELGA', numeStore, chvari, iret)
            ivari = 1
            if (iret .gt. 0) then
                if (factorKeyword .ne. 'ENER_ELAS' .and. factorKeyword .ne. 'ENER_ELTR') then
                    vk24(1) = result
                    vk24(2) = numeStoreStr
                    call utmess('F', 'UTILITAI3_79', nk=2, valk=vk24)
                else
!                   -- creation d'un champ de variables internes nul
                    ivari = 0
                    chbid = '&&PEINGL.VARINUL'
                    call alchml(modelLigrel, 'TOU_INI_ELGA', 'PVARI_R', 'V', chbid, iret, ' ')
                end if
            end if
        end if
!
! ---  RECUPERATION DU CHAMP DES DEPLACEMENTS ASSOCIE AU
! ---  NUMERO D'ORDRE COURANT POUR ENER_ELAS, ENER_ELTR ET ENER_TOTALE:
!      -----------------------------------------------------
        if (factorKeyword .eq. 'ENER_TOTALE' .or. &
            factorKeyword .eq. 'ENER_ELAS' .or. &
            factorKeyword .eq. 'ENER_ELTR' .or. &
            factorKeyword(1:6) .eq. 'INDIC_') then
            call rsexch('F', result, 'DEPL', numeStore, chdepl, iret)
!
! ---  RECUPERATION DU CHAMP DES DEPLACEMENTS ASSOCIE AU
! ---  NUMERO D'ORDRE PRECEDENT :
!      ----------------------
            if (iStore .gt. 1) then
                call rsexch('F', result, 'DEPL', numeStorePrev, chdepm, iret)
            end if
        end if
!

! ---  CALCUL DE L'INDICATEUR GLOBAL DE PERTE DE RADIALITE
! ---  SUR TOUTES LES MAILLES DU MODELE :
!      --------------------------------
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom(1:19)
        lpain(2) = 'PMATERC'
        lchin(2) = materCodeZ
        lpain(3) = 'PVARIPR'
        if (ivari .eq. 1) then
            lchin(3) = chvari(1:19)
        else
            lchin(3) = ' '
        end if
        lpain(4) = 'PCOMPOR'
        lchin(4) = compor(1:19)
        lpain(5) = 'PVARCPR'
        lchin(5) = chvarc
        lpain(6) = 'PVARCRR'
        lchin(6) = chvref
        nbFieldIn = 6
        if (option .eq. 'ENER_TOTALE') then
            if (iStore .gt. 1) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PCONTMR'
                lchin(nbFieldIn) = chsigm(1:19)
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PDEPLM'
                lchin(nbFieldIn) = chdepm(1:19)
            end if
        end if
        if (option .eq. 'ENER_TOTALE' .or. &
            option .eq. 'ENEL_ELEM' .or. &
            option .eq. 'ENTR_ELEM' .or. &
            factorKeyword(1:6) .eq. 'INDIC_') then
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PDEPLR'
            lchin(nbFieldIn) = chdepl(1:19)
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PCONTPR'
            lchin(nbFieldIn) = chsig(1:19)
        end if

! ----- Add fields for structural elements
        call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElem)

! ----- Add input XFEM fields if required
        if (lXFEM .and. option .eq. "ENEL_ELEM") then
            call xajcin(modelZ, option, nbFieldInMax, lchin, lpain, nbFieldIn)
        end if

        if (option .eq. 'INDIC_ENER' .or. option .eq. 'INDIC_SEUIL') then
            nbFieldOut = 2
            lpaout(1) = 'PENERD1'
            lchout(1) = '&&PEINGL.INDIC'
            lpaout(2) = 'PENERD2'
            lchout(2) = '&&PEINGL.VOLUME'
        else if (option .eq. 'ENEL_ELEM' .or. option .eq. 'ENER_TOTALE') then
            nbFieldOut = 1
            lpaout(1) = 'PENERD1'
            lchout(1) = '&&PEINGL.INDIC'
        else if (option .eq. 'ENTR_ELEM') then
            nbFieldOut = 1
            lpaout(1) = 'PENTRD1'
            lchout(1) = '&&PEINGL.INDIC'
        else if (option .eq. 'DISS_ELEM') then
            nbFieldOut = 1
            lpaout(1) = 'PDISSD1'
            lchout(1) = '&&PEINGL.INDIC'
        end if
!
        call calcul('S', option, ligrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    'V', 'OUI')

        do iFactorKeyword = 1, nbFactorKeyword
            work = 0.d0
            valr(2:6) = 0.d0
            deca = 0
!
! ---   RECUPERATION DES MAILLES POUR LESQUELLES ON VA CALCULER
! ---   L'INDICATEUR :
!       ------------
            call getvtx(factorKeyword, 'TOUT', iocc=iFactorKeyword, nbval=0, nbret=nt)
            call getvem(mesh, 'MAILLE', factorKeyword, 'MAILLE', iFactorKeyword, &
                        0, k8b, nm)
            call getvem(mesh, 'GROUP_MA', factorKeyword, 'GROUP_MA', iFactorKeyword, &
                        0, k8b, ng)
!
! ---   TRAITEMENT DU MOT CLE "TOUT" ,LA QUANTITE EST CALCULEE
! ---   SUR TOUT LE MODELE :
!       ------------------
            if (nt .ne. 0) then
                if (factorKeyword .eq. 'INDIC_ENER' .or. &
                    factorKeyword .eq. 'INDIC_SEUIL') then
!
! ---     SOMMATION DES INTEGRALES SUIVANTES SUR LE MODELE
! ---     LA PREMIERE INTEGRALE CALCULEE EST :
! ---     SOMME_DOMAINE((1 - PSI(EPS)/OMEGA(EPS,VARI)).DV
! ---     LA SECONDE INTEGRALE CALCULEE EST LE VOLUME :
!         -------------------------------------------
                    call mesomm(lchout(1), 1, vr=work(1))
                    call mesomm(lchout(2), 1, vr=work(2))
!
                    indic1 = work(1)
                    volume = work(2)
!
                    if (indic1 .le. 1.0d4*r8prem()) then
                        indic1 = zero
                    end if
!
                    if (volume .le. r8prem()) then
                        call utmess('F', 'UTILITAI3_80')
                    end if
!
                    valr(2) = indic1/volume
                    vk8(1) = mesh
                    vk8(2) = 'TOUT'
!
                else if (factorKeyword .eq. 'ENER_ELAS' .or. &
                         factorKeyword .eq. 'ENER_ELTR' .or. &
                         factorKeyword .eq. 'ENER_TOTALE' .or. &
                         factorKeyword .eq. 'ENER_DISS') then
!
! ---          SOMMATION DE L'ENERGIE ( ELASTIQUE OU TOTALE)
! ---          SUR LE MODELE :
!              -------------
                    if (factorKeyword .eq. 'ENER_TOTALE' .or. &
                        factorKeyword .eq. 'ENER_DISS') then
                        call mesomm(lchout(1), 1, vr=work(1))
                    else
                        call mesomm(lchout(1), 5, vr=work(1))
                    end if
! ---  BOUCLE SUR LES PAS DE TEMPS ON SOMME LES TERMES DE
! ---  L ENERGIE TOTAL
                    if ((compt(1:9) .ne. 'VMIS_ISOT') .and. (compt(1:4) .ne. 'ELAS') .and. &
                        (factorKeyword .ne. 'ENER_ELAS' .and. &
                         factorKeyword .ne. 'ENER_ELTR') .and. &
                        (factorKeyword .ne. 'ENER_DISS')) then
                        energy_tout = energy_tout+work(1)
                    else
                        energy_tout = work(1)
                    end if
!
                    vk8(1) = mesh
                    vk8(2) = 'TOUT'
                    valr(2) = energy_tout
                    if (factorKeyword .eq. 'ENER_ELAS' .or. &
                        factorKeyword .eq. 'ENER_ELTR' .or. &
                        factorKeyword .eq. 'ENER_TOTALE') then
                        valr(3) = work(2)
                        valr(4) = work(3)
                        if (factorKeyword .eq. 'ENER_ELAS') then
! ---    AJOUT INUTILE POUR L INSTANT PUISQUE WORK(4) ET WORK(5)
!        SONT NULS. EN PREVISION DU CALCUL DE L ENERGIE ELASTIQUE
!        DE CISAILLEMENT ET DE COUPLAGE MEMBRANE FLEXION POUR LES
!        PLAQUES EN MECA STATIQUE UNIQUEMENT, SI ON L AUTORISE
!        UN JOUR.
                            valr(5) = work(4)
                            valr(6) = work(5)
                        end if
                    end if
!
                end if
!
! ---    ECRITURE DE L'INDICATEUR OU DE L'ENERGIE DANS LA TABLE :
!        ------------------------------------------------------
                call tbajli(tablOutZ, nbParaResu, tablParaName, [numeStore], valr, &
                            [c16b], vk8, 0)
            end if
!
! ---   TRAITEMENT DU MOT CLE "GROUP_MA" ,LA QUANTITE EST CALCULEE
! ---   SUR LE GROUP_MA COURANT :
!       -----------------------
            if (ng .ne. 0) then
                nbgrma = -ng
                call wkvect('&&PEINGL_GROUPM', 'V V K24', nbgrma, jgr)
                call getvem(mesh, 'GROUP_MA', factorKeyword, 'GROUP_MA', iFactorKeyword, &
                            nbgrma, zk24(jgr), ng)
!
! ---     BOUCLE SUR LES GROUPES DE MAILLES :
!         ---------------------------------
                vk24(2) = 'GROUP_MA'
                do ig = 1, nbgrma
                    nomgrm = zk24(jgr+ig-1)
                    call jeexin(jexnom(mesh//'.GROUPEMA', nomgrm), iret)
                    if (iret .eq. 0) then
                        call utmess('F', 'UTILITAI3_46', sk=nomgrm)
                    end if
                    call jelira(jexnom(mesh//'.GROUPEMA', nomgrm), 'LONUTI', nbCell)
                    if (nbCell .eq. 0) then
                        call utmess('F', 'UTILITAI3_47', sk=nomgrm)
                    end if
                    call jeveuo(jexnom(mesh//'.GROUPEMA', nomgrm), 'L', jad)
!
                    if (factorKeyword .eq. 'INDIC_ENER' .or. &
                        factorKeyword .eq. 'INDIC_SEUIL') then
!
! ---      SOMMATION DES INTEGRALES SUIVANTES SUR LES
! ---      MAILLES DU GROUP_ MA
! ---      LA PREMIERE INTEGRALE CALCULEE EST :
! ---      SOMME_DOMAINE((1 - PSI(EPS)/OMEGA(EPS,VARI)).DV
! ---      LA SECONDE INTEGRALE CALCULEE EST LE VOLUME :
!          -------------------------------------------
                        call mesomm(lchout(1), 1, vr=work(1), nbma=nbCell, linuma=zi(jad))
                        call mesomm(lchout(2), 1, vr=work(2), nbma=nbCell, linuma=zi(jad))
!
                        indic1 = work(1)
                        volume = work(2)
!
                        if (indic1 .le. 1.0d4*r8prem()) then
                            indic1 = zero
                        end if
!
                        if (volume .le. r8prem()) then
                            call utmess('F', 'UTILITAI3_81', sk=nomgrm)
                        end if
!
                        valr(2) = indic1/volume
                        vk24(1) = nomgrm
!
                    else if (factorKeyword .eq. 'ENER_ELAS' .or. &
                             factorKeyword .eq. 'ENER_ELTR' .or. &
                             factorKeyword .eq. 'ENER_TOTALE' .or. &
                             factorKeyword .eq. 'ENER_DISS') then
!
! ---          SOMMATION DE L'ENERGIE ( ELASTIQUE OU TOTALE)
! ---          SUR LE MODELE :
!              -------------
                        deca = deca+1
                        ASSERT(deca < nbgrma_tot)
                        if (factorKeyword .eq. 'ENER_TOTALE' .or. &
                            factorKeyword .eq. 'ENER_DISS') then
                            call mesomm(lchout(1), 1, vr=work(1), nbma=nbCell, linuma=zi(jad))
                        else
                            call mesomm(lchout(1), 5, vr=work, nbma=nbCell, linuma=zi(jad))
                        end if
!
! ---  BOUCLE SUR LES PAS DE TEMPS ON SOMME LES TERMES DE
! ---  L ENERGIE TOTAL
!
                        if ((compt(1:9) .ne. 'VMIS_ISOT') .and. (compt(1:4) .ne. 'ELAS') .and. &
                            (factorKeyword .ne. 'ENER_ELAS' .and. &
                             factorKeyword .ne. 'ENER_ELTR') .and. &
                            (factorKeyword .ne. 'ENER_DISS')) then
!
                            energy_grpma(deca) = energy_grpma(deca)+work(1)
                        else
                            energy_grpma(deca) = work(1)
                        end if
!
                        vk24(1) = nomgrm
                        valr(2) = energy_grpma(deca)
                        if (factorKeyword .eq. 'ENER_ELAS' .or. &
                            factorKeyword .eq. 'ENER_ELTR' .or. &
                            factorKeyword .eq. 'ENER_TOTALE') then
                            valr(3) = work(2)
                            valr(4) = work(3)
                            if (factorKeyword .eq. 'ENER_ELAS') then
! ---    AJOUT INUTILE POUR L INSTANT PUISQUE WORK(4) ET WORK(5)
!        SONT NULS. EN PREVISION DU CALCUL DE L ENERGIE ELASTIQUE
!        DE CISAILLEMENT ET DE COUPLAGE MEMBRANE FLEXION POUR LES
!        PLAQUES EN MECA STATIQUE UNIQUEMENT, SI ON L AUTORISE
!        UN JOUR.
                                valr(5) = work(4)
                                valr(6) = work(5)
                            end if
                        end if
                    end if
!
!
! ---    ECRITURE DE L'INDICATEUR OU DE L'ENERGIE DANS LA TABLE :
!        ------------------------------------------------------
!
! ---      ECRITURE DE L'INDICATEUR DANS LA TABLE :
!          --------------------------------------
                    call tbajli(tablOutZ, nbParaResu, tablParaName, [numeStore], valr, &
                                [c16b], vk24, 0)
                end do
!
! --- UNION
                if (nbgrma > 1) then
                    nomgrm = "UNION_GROUP_MA"
                    call umalma(mesh, zk24(jgr), nbgrma, v_allma, nbtot)
                    ASSERT(nbtot > 0)
                    if (factorKeyword .eq. 'INDIC_ENER' .or. &
                        factorKeyword .eq. 'INDIC_SEUIL') then
!
! ---      SOMMATION DES INTEGRALES SUIVANTES SUR LES
! ---      MAILLES DU GROUP_ MA
! ---      LA PREMIERE INTEGRALE CALCULEE EST :
! ---      SOMME_DOMAINE((1 - PSI(EPS)/OMEGA(EPS,VARI)).DV
! ---      LA SECONDE INTEGRALE CALCULEE EST LE VOLUME :
!          -------------------------------------------
                        call mesomm(lchout(1), 1, vr=work(1), nbma=nbtot, linuma=v_allma)
                        call mesomm(lchout(2), 1, vr=work(2), nbma=nbtot, linuma=v_allma)
!
                        indic1 = work(1)
                        volume = work(2)
!
                        if (indic1 .le. 1.0d4*r8prem()) then
                            indic1 = zero
                        end if
!
                        if (volume .le. r8prem()) then
                            call utmess('F', 'UTILITAI3_81', sk=nomgrm)
                        end if
!
                        valr(2) = indic1/volume
                        vk24(1) = nomgrm
!
                    else if (factorKeyword .eq. 'ENER_ELAS' .or. &
                             factorKeyword .eq. 'ENER_ELTR' .or. &
                             factorKeyword .eq. 'ENER_TOTALE' .or. &
                             factorKeyword .eq. 'ENER_DISS') then
!
! ---          SOMMATION DE L'ENERGIE ( ELASTIQUE OU TOTALE)
! ---          SUR LE MODELE :
!              -------------
                        deca = deca+1
                        ASSERT(deca <= nbgrma_tot)
                        if (factorKeyword .eq. 'ENER_TOTALE' .or. &
                            factorKeyword .eq. 'ENER_DISS') then
                            call mesomm(lchout(1), 1, vr=work(1), nbma=nbtot, linuma=v_allma)
                        else
                            call mesomm(lchout(1), 5, vr=work, nbma=nbtot, linuma=v_allma)
                        end if
!
! ---  BOUCLE SUR LES PAS DE TEMPS ON SOMME LES TERMES DE
! ---  L ENERGIE TOTAL
!
                        if ((compt(1:9) .ne. 'VMIS_ISOT') .and. (compt(1:4) .ne. 'ELAS') .and. &
                            (factorKeyword .ne. 'ENER_ELAS' .and. &
                             factorKeyword .ne. 'ENER_ELTR') .and. &
                            (factorKeyword .ne. 'ENER_DISS')) then
!
                            energy_grpma(deca) = energy_grpma(deca)+work(1)
                        else
                            energy_grpma(deca) = work(1)
                        end if
!
                        vk24(1) = nomgrm
                        valr(2) = energy_grpma(deca)
                        if (factorKeyword .eq. 'ENER_ELAS' .or. &
                            factorKeyword .eq. 'ENER_ELTR' .or. &
                            factorKeyword .eq. 'ENER_TOTALE') then
                            valr(3) = work(2)
                            valr(4) = work(3)
                            if (factorKeyword .eq. 'ENER_ELAS') then
! ---    AJOUT INUTILE POUR L INSTANT PUISQUE WORK(4) ET WORK(5)
!        SONT NULS. EN PREVISION DU CALCUL DE L ENERGIE ELASTIQUE
!        DE CISAILLEMENT ET DE COUPLAGE MEMBRANE FLEXION POUR LES
!        PLAQUES EN MECA STATIQUE UNIQUEMENT, SI ON L AUTORISE
!        UN JOUR.
                                valr(5) = work(4)
                                valr(6) = work(5)
                            end if
                        end if
                    end if
!
!
! ---    ECRITURE DE L'INDICATEUR OU DE L'ENERGIE DANS LA TABLE :
!        ------------------------------------------------------
!
! ---      ECRITURE DE L'INDICATEUR DANS LA TABLE :
!          --------------------------------------
                    call tbajli(tablOutZ, nbParaResu, tablParaName, [numeStore], valr, &
                                [c16b], vk24, 0)
                    AS_DEALLOCATE(vi=v_allma)
                end if
!
                call jedetr('&&PEINGL_GROUPM')
            end if
!
! ---   TRAITEMENT DU MOT CLE "MAILLE" ,L'INDICATEUR EST CALCULE
! ---   SUR LA MAILLE COURANTE :
!       ----------------------
            if (nm .ne. 0) then
                nbCellil = -nm
                call wkvect('&&PEINGL_MAILLE', 'V V K8', nbCellil, jma)
                call getvem(mesh, 'MAILLE', factorKeyword, 'MAILLE', iFactorKeyword, &
                            nbCellil, zk8(jma), nm)
!
! ---    BOUCLE SUR LES MAILLES :
!        ----------------------
                vk8(2) = 'MAILLE'
                call jelira(mesh//'.TYPMAIL', 'LONMAX', nbMaiT)
                do im = 1, nbCellil
                    nommai = zk8(jma+im-1)
                    nume = char8_to_int(nommai)
                    if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                        call utmess('A', 'UTILITAI3_49', sk=zk8(jma+im-1))
                    end if
!
                    if (factorKeyword .eq. 'INDIC_ENER' .or. &
                        factorKeyword .eq. 'INDIC_SEUIL') then
!
! ---      LES INTEGRALES SONT CALCULEES SUR LA MAILLE COURANTE
! ---      LA PREMIERE INTEGRALE CALCULEE EST :
! ---      SOMME_DOMAINE((1 - PSI(EPS)/OMEGA(EPS,VARI)).DV
! ---      LA SECONDE INTEGRALE CALCULEE EST LE VOLUME :
!          -------------------------------------------
                        call mesomm(lchout(1), 1, vr=work(1), nbma=1, linuma=[nume])
                        call mesomm(lchout(2), 1, vr=work(2), nbma=1, linuma=[nume])
!
                        indic1 = work(1)
                        volume = work(2)
!
                        if (indic1 .le. 1.0d4*r8prem()) then
                            indic1 = zero
                        end if
!
                        if (volume .le. r8prem()) then
                            call utmess('F', 'UTILITAI3_82', sk=nommai)
                        end if
!
                        valr(2) = indic1/volume
                        vk8(1) = nommai
!
                    else if (factorKeyword .eq. 'ENER_ELAS' .or. &
                             factorKeyword .eq. 'ENER_ELTR' .or. &
                             factorKeyword .eq. 'ENER_TOTALE' .or. &
                             factorKeyword .eq. 'ENER_DISS') then
!
! ---          SOMMATION DE L'ENERGIE ( ELASTIQUE OU TOTALE)
! ---          SUR LE MODELE :
!              -------------
                        if (factorKeyword .eq. 'ENER_TOTALE' .or. &
                            factorKeyword .eq. 'ENER_DISS') then
                            call mesomm(lchout(1), 1, vr=work(1), nbma=1, linuma=[nume])
                        else
                            call mesomm(lchout(1), 5, vr=work, nbma=1, linuma=[nume])
                        end if
!
                        if ((compt(1:9) .ne. 'VMIS_ISOT') .and. (compt(1:4) .ne. 'ELAS') .and. &
                            (factorKeyword .ne. 'ENER_ELAS' .and. &
                             factorKeyword .ne. 'ENER_ELTR')) then
!
                            ASSERT(nbCellil == 1)
                            energy_ma = energy_ma+work(1)
                        else
                            energy_ma = work(1)
                        end if
!
                        valr(2) = energy_ma
                        vk8(1) = nommai
                        if (factorKeyword .eq. 'ENER_ELAS' .or. &
                            factorKeyword .eq. 'ENER_ELTR' .or. &
                            factorKeyword .eq. 'ENER_TOTALE') then
                            valr(3) = work(2)
                            valr(4) = work(3)
                            if (factorKeyword .eq. 'ENER_ELAS') then
! ---    AJOUT INUTILE POUR L INSTANT PUISQUE WORK(4) ET WORK(5)
!        SONT NULS. EN PREVISION DU CALCUL DE L ENERGIE ELASTIQUE
!        DE CISAILLEMENT ET DE COUPLAGE MEMBRANE FLEXION POUR LES
!        PLAQUES EN MECA STATIQUE UNIQUEMENT, SI ON L AUTORISE
!        UN JOUR.
                                valr(5) = work(4)
                                valr(6) = work(5)
                            end if
                        end if
!
                    end if
!
! ---      ECRITURE DE L'INDICATEUR DANS LA TABLE :
!          --------------------------------------
                    call tbajli(tablOutZ, nbParaResu, tablParaName, [numeStore], valr, &
                                [c16b], vk8, 0)
                end do
!
                call jedetr('&&PEINGL_MAILLE')
            end if
        end do
        call jedetr('&&MECHTI.CH_INST_R')
        call detrsd('CHAM_ELEM', chvarc)
        call detrsd('CHAM_ELEM', chvref)
        call jedetr(compor//'.PTMA')
        call jedema()
    end do
70  continue
    call jedetr(listStoreJv)
    call jedetr(listTimeJv)
    call jedetr('&&PEINGL.INDIC')
    call jedetr('&&PEINGL.VOLUME')
    call jedetr('&&MEHARM.NUME_HARM')
    if (ivari .eq. 0) then
        call jedetr(chbid)
    end if
!
80  continue
    AS_DEALLOCATE(vr=energy_grpma)
    call jedema()
end subroutine
