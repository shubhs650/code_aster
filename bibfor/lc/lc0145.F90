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
! aslint: disable=W0104
!
subroutine lc0145(fami, kpg, ksp, ndim, imate, &
                  carcri, instam, instap, epsm, &
                  deps, sigm, nvi, vim, option, &
                  sigp, vip, typmod, &
                  dsidep, codret)

!
! ----------------------------------------------------------------------
!
!              LOI BETON_RAG
!
! ----------------------------------------------------------------------
!
    use beton_rag_module
    use tenseur_meca_module
    implicit none
!
#include "asterc/r8t0.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/get_varc.h"
#include "asterfort/rcexistvarc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: imate, ndim, kpg, ksp, nvi
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam, instap
    real(kind=8), intent(in) :: epsm(6), deps(6), sigm(6), vim(nvi)
    character(len=16), intent(in) ::  option
    character(len=8), intent(in) :: typmod(*)
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(out) :: codret
    real(kind=8), intent(out) :: sigp(6)
    real(kind=8), intent(out) :: vip(nvi)
    real(kind=8), intent(out) :: dsidep(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
    type(beton_rag_materiau)   :: mater_bet_rag
    type(beton_rag_parametres) :: param_bet_rag

    real(kind=8)        ::  depsldc(6), epsmldc(6), sigmldc(6)

    integer(kind=8), parameter   :: nbval = 12
    integer(kind=8)             :: icodre(nbval), erre_compat_hypotheses
    character(len=16)   :: nomres(nbval)
    real(kind=8)        :: valres(nbval)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)       :: ii, jj, iret, jcret, nbres
    real(kind=8)  :: epsmeca(6), perturb, vperturb(6), NormSigm, numerateur
    real(kind=8)  :: sigptb(6), viptb(1), dsideptb(6, 6)
    aster_logical :: rigi, resi, elas, isnogood
!
    real(kind=8)        :: valr(8)
    character(len=30)   :: valk(5)
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
!   RIGI_MECA_TANG ->        DSIDEP        -->  RIGI
!   FULL_MECA      ->  SIGP  DSIDEP  VARP  -->  RIGI  RESI
!   RAPH_MECA      ->  SIGP          VARP  -->        RESI
!   RIGI_MECA_ELAS ->        ELAS          -->  RIGI        ELAS
!   FULL_MECA_ELAS ->  SIGP  ELAS    VARP  -->  RIGI  RESI  ELAS

    rigi = (option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RIGI')
    resi = (option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RAPH')
    elas = option(11:14) .eq. 'ELAS'
    !
    ! Incrément de temps
    param_bet_rag%instap = instap
    param_bet_rag%instam = instam
    param_bet_rag%dtemps = instap-instam
    ! Nombre d'itérations maxi (ITER_INTE_MAXI)
    param_bet_rag%nbdecp = int(carcri(1))
    ! Tolérance de convergence (RESI_INTE)
    param_bet_rag%errmax = carcri(3)
    !
    param_bet_rag%rigi = rigi
    param_bet_rag%resi = resi
    !
    ! Récupération de la TEMPÉRATURE
    call get_varc(fami, kpg, ksp, 'T', &
                  param_bet_rag%temperm, param_bet_rag%temperp, &
                  param_bet_rag%temperref, param_bet_rag%istemper)
    !
    if (param_bet_rag%istemper) then
        param_bet_rag%dtemper = param_bet_rag%temperp-param_bet_rag%temperm
    else
        param_bet_rag%dtemper = 0.0d0
        param_bet_rag%temperm = 0.0d0
        param_bet_rag%temperp = 0.0d0
        param_bet_rag%temperref = 0.0d0
    end if
    !
    param_bet_rag%ishydrat = rcexistvarc('HYDR')
    param_bet_rag%issechag = rcexistvarc('SECH')
    !
    ! Caractéristiques élastiques
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nbres = 2
    if (param_bet_rag%istemper) then
        nbres = nbres+1
        nomres(nbres) = 'ALPHA'
    end if
    if (param_bet_rag%ishydrat) then
        nbres = nbres+1
        nomres(nbres) = 'B_ENDOGE'
    end if
    call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.0D0], &
                nbres, nomres, valres, icodre, 1)
    ! Mécanique
    mater_bet_rag%young = valres(1)
    mater_bet_rag%nu = valres(2)
    nbres = 2
    if (param_bet_rag%istemper) then
        nbres = nbres+1
        mater_bet_rag%alpha = valres(nbres)
    end if
    if (param_bet_rag%ishydrat) then
        nbres = nbres+1
        mater_bet_rag%bendo = valres(nbres)
    end if
    ! Variables de commandes
    if (param_bet_rag%ishydrat) then
        call rcvarc('F', 'HYDR', '-', fami, kpg, ksp, param_bet_rag%hydratm, iret)
        call rcvarc('F', 'HYDR', '+', fami, kpg, ksp, param_bet_rag%hydratp, iret)
        param_bet_rag%dhydrat = param_bet_rag%hydratp-param_bet_rag%hydratm
    end if
    if (param_bet_rag%issechag) then
        call rcvarc('F', 'SECH', '-', fami, kpg, ksp, param_bet_rag%sechagm, iret)
        call rcvarc('F', 'SECH', '+', fami, kpg, ksp, param_bet_rag%sechagp, iret)
        ! En dessous de BR_SECHAGE_MINI cela ne veut plus rien dire
        if ((param_bet_rag%sechagm < BR_SECHAGE_MINI) .or. &
            (param_bet_rag%sechagp < BR_SECHAGE_MINI)) then
            valr(1) = instap
            valr(2) = 0.10
            valr(3) = param_bet_rag%sechagm
            valr(4) = param_bet_rag%sechagp
            call utmess('F', 'COMPOR3_50', nr=4, valr=valr)
        end if
        param_bet_rag%dsechag = param_bet_rag%sechagp-param_bet_rag%sechagm
    end if
    !
    nomres(1) = 'COMP_BETON'
    nomres(2) = 'ENDO_MC'
    nomres(3) = 'ENDO_SIGUC'
    nomres(4) = 'ENDO_MT'
    nomres(5) = 'ENDO_SIGUT'
    nomres(6) = 'ENDO_DRUPRA'
    call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'BETON_RAG', 0, ' ', [0.0D0], &
                6, nomres, valres, icodre, 1)
    !
    param_bet_rag%loi_integre = nint(valres(1))
    mater_bet_rag%mc = valres(2)
    mater_bet_rag%siguc = valres(3)
    mater_bet_rag%mt = valres(4)
    mater_bet_rag%sigut = valres(5)
    mater_bet_rag%dhom = valres(6)
    ! vérification des données
    isnogood = (valres(2) <= 0.0) .or. (valres(4) <= 0.0) .or. (valres(6) < 0.0)
    if (isnogood) then
        valk(1) = 'ENDO_MC ENDO_MT ENDO_DRUPRA'
        call utmess('F', 'COMPOR3_54', nk=1, valk=valk)
    end if

    ! si on fait au moins du fluage
    if (param_bet_rag%loi_integre .ge. BR_LOI_FLUA_LIN) then
        nomres(1) = 'FLUA_KAPPAI'
        nomres(2) = 'FLUA_SPH_KR'
        nomres(3) = 'FLUA_SPH_NR'
        nomres(4) = 'FLUA_SPH_NI'
        nomres(5) = 'FLUA_DEV_KR'
        nomres(6) = 'FLUA_DEV_NR'
        nomres(7) = 'FLUA_DEV_NI'
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'BETON_RAG', 0, ' ', [0.0D0], &
                    7, nomres, valres, icodre, 0)
        ! On vérifie que l'utilisateur a précisé une valeur de kappa_i
        if (icodre(1) .eq. 0) then
            mater_bet_rag%kappa_fluage = valres(1)
            param_bet_rag%loi_integre = param_bet_rag%loi_integre+1
        else
            mater_bet_rag%kappa_fluage = 1.
        end if
        ! Fluage sphérique
        mater_bet_rag%fluage_sph%kr = valres(2)
        mater_bet_rag%fluage_sph%nr = valres(3)
        mater_bet_rag%fluage_sph%ni = valres(4)
        ! Fluage déviatorique
        mater_bet_rag%fluage_dev%kr = valres(5)
        mater_bet_rag%fluage_dev%nr = valres(6)
        mater_bet_rag%fluage_dev%ni = valres(7)
        ! vérification des données : soit ils sont tous positifs (cas normal)
        !                            soit tous négatifs (cas spécifique de test)
        isnogood = XOR((mater_bet_rag%kappa_fluage <= 0.0) .or. (valres(2) <= 0.0) &
                       .or. (valres(3) <= 0.0) .or. (valres(4) <= 0.0) &
                       .or. (valres(5) <= 0.0) .or. (valres(6) <= 0.0) &
                       .or. (valres(7) <= 0.0) &
                       , &
                       (valres(2) < 0.0) &
                       .and. (valres(3) < 0.0) .and. (valres(4) < 0.0) &
                       .and. (valres(5) < 0.0) .and. (valres(6) < 0.0) &
                       .and. (valres(7) < 0.0))
        ! isnogood = (mater_bet_rag%kappa_fluage <= 0.0) .or. (valres(2) <= 0.0) &
        !             .or. (valres(3) <= 0.0) .or. (valres(4) <= 0.0) &
        !             .or. (valres(5) <= 0.0) .or. (valres(6) <= 0.0) &
        !             .or. (valres(7) <= 0.0)
        if (isnogood) then
            valk(1) = 'FLUA_SPH_* FLUA_DEV_*'
            call utmess('F', 'COMPOR3_54', nk=1, valk=valk)
        end if
    end if

    ! Si l'utilisateur a demandé à faire de la RAG, on récupère les paramètres
    if (param_bet_rag%loi_integre .ge. BR_LOI_RAG_FLIN) then
        nomres(1) = 'GEL_ALPHA0'
        nomres(2) = 'GEL_TREF'
        nomres(3) = 'GEL_EAR'
        nomres(4) = 'GEL_SR0'
        nomres(5) = 'GEL_VG'
        nomres(6) = 'GEL_MG'
        nomres(7) = 'GEL_BG'
        nomres(8) = 'GEL_A0'
        nomres(9) = 'RAG_EPSI0'
        nomres(10) = 'PW_A'
        nomres(11) = 'PW_B'
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'BETON_RAG', 0, ' ', [0.0D0], &
                    11, nomres, valres, icodre, 1)
        !  Avancement du gel
        mater_bet_rag%gel%alpha0 = valres(1)
        mater_bet_rag%gel%tref = valres(2)
        mater_bet_rag%gel%ear = valres(3)
        mater_bet_rag%gel%sr0 = valres(4)
        ! vérification des données
        isnogood = (valres(1) < 0.0) .or. (valres(2) <= -r8t0()) .or. (valres(3) <= 0.0)
        if (isnogood) then
            valk(1) = 'GEL_ALPHA0 GEL_TREF GEL_EAR'
            call utmess('F', 'COMPOR3_54', nk=1, valk=valk)
        end if
        ! Protection sur les valeurs de séchage initial
        isnogood = (valres(4) < BR_SECHAGE_MINI) .or. (valres(4) >= BR_SECHAGE_MAXI)
        if (isnogood) then
            valr(1) = valres(4)
            valr(2) = BR_SECHAGE_MINI
            valr(3) = BR_SECHAGE_MAXI
            valk(1) = 'GEL_SR0'
            call utmess('F', 'COMPOR3_53', nr=3, valr=valr, nk=1, valk=valk)
        end if
        ! Pression du gel
        mater_bet_rag%gel%vg = valres(5)
        mater_bet_rag%gel%mg = valres(6)
        mater_bet_rag%gel%bg = valres(7)
        mater_bet_rag%gel%a0 = valres(8)
        ! vérification des données
        isnogood = (valres(5) <= 0.0) .or. (valres(6) <= 0.0) .or. (valres(7) < 0.0)
        if (isnogood) then
            valk(1) = 'GEL_VG GEL_MG GEL_BG'
            call utmess('F', 'COMPOR3_54', nk=1, valk=valk)
        end if
        ! Protection sur les valeurs du seuil d'avancement
        isnogood = (valres(8) < 0.0) .or. (valres(8) >= 0.9999)
        if (isnogood) then
            valr(1) = valres(8)
            valr(2) = 0.0
            valr(3) = 0.9999
            valk(1) = 'GEL_A0'
            call utmess('F', 'COMPOR3_53', nr=3, valr=valr, nk=1, valk=valk)
        end if
        ! Déformation visqueuse RAG
        mater_bet_rag%gel%epsi0 = valres(9)
        ! vérification des données
        isnogood = (valres(9) < 0.0)
        if (isnogood) then
            valk(1) = 'RAG_EPSI0'
            call utmess('F', 'COMPOR3_54', nk=1, valk=valk)
        end if
        ! Coefficients Van Genuchten
        mater_bet_rag%pw%a = valres(10)
        mater_bet_rag%pw%b = valres(11)
        ! vérification des données
        isnogood = (valres(10) < 0.0) .or. (valres(11) <= 1.0)
        if (isnogood) then
            valk(1) = 'PW_A PW_B'
            call utmess('F', 'COMPOR3_54', nk=1, valk=valk)
        end if

    end if

    ! protection développeur : loi_integre doit se trouver dans les valeurs acceptées
    isnogood = (param_bet_rag%loi_integre >= BR_LOI_MECA) .and. &
               (param_bet_rag%loi_integre <= BR_LOI_RAG_FNL)
    ASSERT(isnogood)

    ! erreur = 0 : tout va bien, pas d'erreur
    ! erreur = 1 : si on a commencé à faire de la RAG on doit continuer à faire de la RAG
    ! erreur = 5 : si on a commencé un calcul avec du fluage linéaire/non-linéaire, alors
    !                   il faut continuer avec la même hypothèse
    erre_compat_hypotheses = BR_Compat_Hypotheses_Calcul(vim, param_bet_rag%loi_integre)
    if (erre_compat_hypotheses .eq. 1) then
        call utmess('F', 'COMPOR3_51')
    else if (erre_compat_hypotheses .eq. 5) then
        call utmess('F', 'COMPOR3_55')
    end if
    !
    ! Si on fait de la RAG il faut les champs Temper et Sech
    if ((param_bet_rag%loi_integre == BR_LOI_RAG_FLIN) .or. &
        (param_bet_rag%loi_integre == BR_LOI_RAG_FNL)) then
        ! ASSERT( param_bet_rag%issechag .and. param_bet_rag%istemper )
        if (.not. (param_bet_rag%issechag .and. param_bet_rag%istemper)) then
            call utmess('F', 'COMPOR3_52')
        end if
    end if
    !
    !
    if (resi) then
        depsldc = VecteurAsterVecteur(deps)
        vip(1:BR_VARI_NOMBRE) = vim(1:BR_VARI_NOMBRE)
    else
        depsldc = 0.0d0
    end if
    !
    param_bet_rag%perturbation = .false.
    !
    epsmldc = VecteurAsterVecteur(epsm)
    sigmldc = VecteurAsterVecteur(sigm)
    call ldc_beton_rag(epsmldc, depsldc, sigmldc, vim, mater_bet_rag, param_bet_rag, &
                       sigp, vip, dsidep, iret)
    !
    if (iret .ne. 0) then
        goto 900
    end if
    if (elas) then
        goto 900
    end if
    !
    perturb = maxval(abs(depsldc))
    if (rigi .and. perturb > 1.0D-07) then
        ! Le calcul de la matrice tangente par perturbation
        param_bet_rag%perturbation = ASTER_TRUE
        ! Le calcul ne se fait que sur la mécanique
        param_bet_rag%loi_integre = BR_LOI_MECA
        sigmldc = VecteurAsterVecteur(sigp)
        ! Perturbation des déformations
        perturb = 1.0D-06
        !
        NormSigm = 10.0**(nint(log10(abs(mater_bet_rag%siguc))-3))
        ! Déformations Mécanique
        epsmeca = epsmldc+depsldc
        do ii = 1, 6
            vperturb(:) = 0.0d0
            ! La perturbation est dans la direction de l'incrément de déformation
            vperturb(ii) = sign(perturb, depsldc(ii))
            call ldc_beton_rag(epsmeca, vperturb, sigmldc, vip, &
                               mater_bet_rag, param_bet_rag, &
                               sigptb, viptb, dsideptb, iret)
            do jj = 1, 6
                numerateur = abs(sigptb(jj)-sigp(jj))
                if (numerateur/NormSigm > 1.0D-03) then
                    dsidep(jj, ii) = numerateur/perturb
                end if
            end do
        end do
    end if
    !
900 continue
    if (resi) then
        call tecach('NNO', 'PCODRET', 'E', jj, iad=jcret)
        if (jj .ne. 0) then
            zi(jcret) = iret
        end if
    end if
end subroutine lc0145
