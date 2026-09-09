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
! aslint: disable=W0413,W1501,W1306
!
subroutine lcejdm(BEHInteg, fami, kpg, ksp, ndim, &
                  mate, option, epsm, deps, sigma, &
                  dsidep, vim, vip, typmod, instam, &
                  instap)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/matinv.h"
#include "asterfort/pmavec.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "asterfort/finlf1.h"
#include "asterfort/finlf2.h"
    type(Behaviour_Integ), intent(in) :: BEHInteg
    integer(kind=8) :: mate, ndim, kpg, ksp
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: sigma(6), dsidep(6, 6)
    real(kind=8) :: vim(*), vip(*), instam, instap
    character(len=8) :: typmod(*)
    character(len=16) :: option
    character(len=*) :: fami
!
!-----------------------------------------------------------------------
!     LOI DE COMPORTEMENT DES JOINTS DE BARRAGE : JOINT_MECA_ENDO
!     POUR LES ELEMENTS DE JOINT ET JOINT_HYME 2D ET 3D
!
! DEUX TYPES DE CALCUL SONT POSSIBLES:
!
! 1) MODELISATIONS *_JOINT_HYME, AVEC PARAMETRE HYDRO PRESENTS
!    CALCUL COUPLE HYDRO-MECANIQUE SUR UN MAILLAGE QUADRATIQUE
!
! 2) MODELISATIONS *_JOINT, PAS DE PARAMETRES HYDRO POSSIBLE
!    CALCUL MECA SUR UN MAILLAGE LINEAIRE OU QUADRATIQUE
!
! IN : EPSM - SAUT INSTANT MOINS ET GRAD PRESSION ET PRES FLUIDE SI HYME
! IN : MATE, OPTION, VIM, COOROT, INSTAM, INSTAP
! IN : SIGMO - SIGMA INSTANT MOINS ET FLUX HYDRO SI HYME
! OUT : SIGMA , DSIDEP , VIP
!-----------------------------------------------------------------------
    integer(kind=8) :: nbpa
    parameter(nbpa=11)
    integer(kind=8) :: cod(nbpa), i, j, n, kronec
! Indicateur de plastification
    integer(kind=8) :: ifplas
! Parametres meca de la loi
    real(kind=8) :: kn, kt, mu, cbar
    real(kind=8) :: m1, m2, D1, Bn, Bt, CC, TT
! Parametres code_aster
    real(kind=8) :: val(nbpa), valpar(ndim+1), inst, valr(3)
! Parametres hydro de la loi
    real(kind=8) :: rhof, visf, amin
    real(kind=8) :: presfl, presg
! Transfert repere global vers local
    real(kind=8) :: gp(ndim-1), gploc(ndim), gpglo(ndim), fhloc(ndim), fhglo(ndim)
    real(kind=8) :: coorot(ndim+ndim*ndim), invrot(ndim, ndim)
! Parametre d'ouverture
    real(kind=8) :: delta(ndim)
! Parametres plastique
    real(kind=8) :: dlam, plasti(ndim), dplasti(ndim), criter
! Force thermodynamique generalisee et ses normes
    real(kind=8) :: xvec(ndim), xvecpl(ndim), xnormp, xnprov, xnormt
! Parametre d'endommagement
    real(kind=8) :: alpha
! Parametres de la perte de rigidite adimensionnee et ses derivees
    real(kind=8) :: S_alpha
! Parametres S(alpha) calcule pour alpha = alpha_max et ses derivees
    real(kind=8) :: S_alpha_max, Sp_alpha_max
! Parametres de verification du domaine de validite de la loi
    real(kind=8) :: R_max, D1_min
! Parametres temporaire pour le calcul iteratif de alpha
    real(kind=8) :: alpha1, alpha2, alpha_dicho, alpha_corde, alpha_pro, alpha_pred
    real(kind=8) :: res1, res2, res_dicho, res_corde, res_inter, crit_conv
    character(len=8) :: methode
! Critere de regularisation de la phase d'initilisation d'endommagement
    real(kind=8) :: alpha_min
! Autres parametres de la loi
    real(kind=8) :: alpha_max
! Parametres numeriques de precision
    real(kind=8) :: r8bid
    character(len=8) :: nompar(ndim+1)
    character(len=16) :: nom(nbpa)
    character(len=1) :: poum
! Parametre numerique pour le calcul de la matrice tangente
    real(kind=8) :: An, Apn, A2pn, At, Apt, A2pt
    real(kind=8) :: W_vect(2), W_norm, pt_norm, ptxw, disc, d_lambda
    real(kind=8) :: ddldalp, dfddn, dfddt(2), dfdalp, ddlddn, ddlddt(2)
    real(kind=8) :: dgddn, dgddt(2), dgdalp, ddpnddn, ddpnddt(2), ddptddn(2), ddptddt(2, 2)
    integer(kind=8) :: identity
!
! Parametre numerique aster
    aster_logical :: resi, rigi, elas, ifpahm, ifhyme
    blas_int :: b_incx, b_incy, b_n
!
! OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
!---------------------------------------------------------
! TODO : activer l'option RIGI_MECA_TANG dans l'element fini
! CALCUL DE CONTRAINTE (RESIDU)
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
! CALCUL DE LA MATRICE TANGEANTE (RIGIDITE)
    rigi = option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA'
! CALCUL DE LA MATRICE ELASTIQUE A LA PLACE DE LA MATRICE TANGENTE
    elas = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'RIGI_MECA_ELAS'
!
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if
!
! INDICATEUR AVEC/SANS HYDRO
    if (typmod(2) .eq. 'EJ_HYME') ifhyme = .true.
    if (typmod(2) .eq. 'ELEMJOIN') ifhyme = .false.
!
! SAUT DE DEPLACEMENT EN T- OU T+
!   Recuperation de delta-
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm, b_incx, delta, b_incy)
!   Calcul de delta+
    if (resi) then
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps, b_incx, delta, &
                   b_incy)
    end if
!
! GRADIENT DE PRESSION ET PRESSION EN T- OU T+
    if (ifhyme) then
        do n = 1, ndim-1
            gp(n) = epsm(ndim+n)
            if (resi) gp(n) = gp(n)+deps(ndim+n)
        end do
        presg = epsm(2*ndim)
        if (resi) presg = presg+deps(2*ndim)
    end if
!
! INSTANT DE CALCUL T- OU T+
    inst = instam
    if (resi) inst = instap
!
! RECUPERATION DES PARAMETRES PHYSIQUES DE LA LDC
!------------------------------------------------
    nom(1) = 'K_N'
    nom(2) = 'MU'
    nom(3) = 'K_T'
    nom(4) = 'PENA_RUPTURE'
    nom(5) = 'Bn'
    nom(6) = 'Bt'
    nom(7) = 'ALPHA'
    nom(8) = 'PRES_FLUIDE'
    nom(9) = 'RHO_FLUIDE'
    nom(10) = 'VISC_FLUIDE'
    nom(11) = 'OUV_MIN'
!
! LECTURE DES PARAMETRES MECA DE LA LDC
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_ENDO', 0, ' ', [0.d0], &
                7, nom, val, cod, 2)
!
! Coefficient de rigidite normale
    kn = val(1)
! Coefficient de frottement
    mu = val(2)
! Coefficient de rigidite tangentielle (si pas renseignee : K_T = K_N)
    if (cod(3) .eq. 0) then
        kt = val(3)
    else
        kt = kn
    end if
! Coefficient pour la fonction d'endommagement
    D1 = val(4)
    m1 = 3.0
    m2 = 0.5
! Coefficients d'endommagement qui fixent la cohesion et la resistance a la traction
    Bn = val(5)
    Bt = val(6)
! Calcul du max de R(alpha) (EQ 2.33) en dur car depend de m1 et m2
    R_max = 0.763244531975
! Cohesion residuelle (EVOLUTION : a faire remonter au catalogue de la loi)
    cbar = 0.d0
! Calcul de alpha_max en fonction de m1 et m2
    alpha_max = (-m2+sqrt(m1*m2/(m1-m2+1.d0)))/(m1-m2)
! Critere de regularisation de la phase d'initilisation d'endommagement (alpha_max/1000) Eq 2.42
    if (cod(7) .eq. 0) then
        alpha_min = val(7)
    else
        alpha_min = alpha_max/1.e5
    end if
! Critere de convergence pour F(alpha) = 0 ou G(alpha) = 0
    crit_conv = alpha_min
! Calcul de S(alpha) pour alpha = alpha_max
    S_alpha_max = ((1.d0-alpha_max)**m1)/(alpha_max**m2)
! Calcul de la derivee S(alpha) pour alpha = alpha_max
    Sp_alpha_max = ((m2-m1)*alpha_max-m2)*((1.d0-alpha_max)**(m1-1))/(alpha_max**(m2+1))
! Calcul de la cohesion
    CC = sqrt(2.d0*D1*(Bn*(mu**2)+Bt))*S_alpha_max/sqrt(-Sp_alpha_max)
! Calcul de la resistance a la traction
    TT = sqrt(2.d0*D1*Bn)*S_alpha_max/sqrt(-Sp_alpha_max)
! Calcul de la valeur minimale de D1
    D1_min = max(TT**2/(2.d0*kn), CC**2/(2.d0*kt))*(-Sp_alpha_max/(S_alpha_max**2))*R_max
!
! VERIFICATION DU DOMAINE DE VALIDITE DES PARAMETRES
!---------------------------------------------------
! Condition sur D1 defini en page 47 de la these
    if (D1 .le. D1_min) then
        valr(1) = D1
        valr(2) = D1_min
        call utmess('F', 'ALGORITH17_33', nr=2, valr=valr)
    end if
! Condition sur CC defini par l'EQ 2.50
    if (CC .lt. mu*TT .or. CC .gt. sqrt(2.d0)*mu*TT) then
        valr(1) = mu*TT
        valr(2) = CC
        valr(3) = sqrt(2.d0)*mu*TT
        call utmess('F', 'ALGORITH17_34', nr=3, valr=valr)
    end if
! Condition sur kn defini par l'EQ 2.46
    if (kn .le. Bn*R_max) then
        valr(1) = kn
        valr(2) = Bn*R_max
        valr(3) = R_max
        call utmess('F', 'ALGORITH17_35', nr=3, valr=valr)
    end if
! Condition sur kt defini par l'EQ 2.41
    if (kt .le. (mu**2*Bn+Bt)*R_max) then
        valr(1) = kt
        valr(2) = (mu**2*Bn+Bt)*R_max
        valr(3) = R_max
        call utmess('F', 'ALGORITH17_36', nr=3, valr=valr)
    end if
! Condition sur Bn et Bt qui decoule de l'EQ 2.19
    if (Bn .le. 0.d0 .or. Bt .le. 0.d0) then
        valr(1) = Bn
        Valr(2) = Bt
        call utmess('F', 'ALGORITH17_37', nr=2, valr=valr)
    end if
!
! DEFINITION DES PARAMETRES POUR LA RECUPERATION DES FONCTIONS
    coorot = 0.d0
    do i = 1, ndim
        coorot(i) = BEHInteg%behavESVA%behavESVAGeom%coorElga(kpg, i)
    end do
    do i = 1, ndim*ndim
        coorot(ndim+i) = BEHInteg%behavESVA%behavESVAOther%rotpg(i)
    end do
!
    nompar(1) = 'INST'
    nompar(2) = 'X'
    nompar(3) = 'Y'
    valpar(1) = inst
    valpar(2) = coorot(1)
    valpar(3) = coorot(2)
    if (ndim .eq. 3) then
        nompar(4) = 'Z'
        valpar(4) = coorot(3)
    end if
!
! RECUPERATION DE LA PRESS FLUIDE (MODELISATION MECA PURE)
!---------------------------------------------------------
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_ENDO', ndim+1, nompar, valpar, &
                1, nom(8), val(8), cod(8), 0)
!
    if (cod(8) .eq. 0) then
        presfl = val(8)
    else
        presfl = 0.d0
    end if
!
! RECUPERATION DE LA MASSE VOL ET DE LA VISCO (MODELISATION JOINT HM)
!--------------------------------------------------------------------
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_ENDO', 0, ' ', [0.d0], &
                1, nom(9), val(9), cod(9), 0)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_ENDO', 0, ' ', [0.d0], &
                1, nom(10), val(10), cod(10), 0)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_ENDO', 0, ' ', [0.d0], &
                1, nom(11), val(11), cod(11), 0)
!
    if (cod(9) .eq. 0) rhof = val(9)
    if (cod(10) .eq. 0) visf = val(10)
    if (cod(11) .eq. 0) amin = val(11)
!
! INDICATEUR SI LES PARAMETRES HYDRO SONT RENSEIGNES
    ifpahm = (cod(9) .eq. 0) .and. (cod(10) .eq. 0) .and. (cod(11) .eq. 0)
!
! VERIFICATION DE LA PRESENCE/ABSENCE DE PARAMETRES
! EN FONCTION DE LA MODELISATION MECA PUR OU HYDRO MECA
!------------------------------------------------------
    if (ifhyme) then
!       POUR LE CALCUL HYDRO => PAS DE PRES_FLUIDE
        if (cod(8) .eq. 0) then
            call utmess('F', 'ALGORITH17_14')
        end if
!       POUR LE CALCUL HYDRO => PRESENCE DE PARA_HM
        if (.not. ifpahm) then
            call utmess('F', 'ALGORITH17_15')
        end if
    else
!       POUR LE CALCUL MECA => PAS DE PARAMETRE HYDRO
        if (ifpahm) then
            call utmess('F', 'ALGORITH17_16')
        end if
    end if
!
! RECUPERATION DES VARIABLES INTERNES
!     PLASTI = vecteur de deplacement plastique normal et tangentiel
!     ALPHA = variable d'endommagement
!-------------------------------------------------------------------
! IDICATEUR DE PLASTIFICATION A L'INSTANT ACTUEL
    ifplas = nint(vim(2))
! Vecteur de plasticite : plasti(1) = vim(3); plasti(2) = vim(4); plasti(3) = vim(5)
    do i = 1, ndim
        plasti(i) = vim(i+2)
    end do
! Parametre d'endommagement
    alpha = vim(6)
! EVOLUTION : Introduire le calcul de plasti(1) par l'eqution fin de page 37
! Afin de travailler sur une structure pré-endo artificiellement
!
! PHASE DE PREDICTION POUR NEWTON
    if (.not. resi) goto 500
!
! PHASE D'INITIALISATION DES VARAIBLES TEMPORAIRES
    dlam = 0.d0
    do i = 1, ndim
        dplasti(i) = 0.d0
    end do
    S_alpha = 0.d0
!   Bornes de recherche pour F(alpha) = 0 ou G(alpha) = 0
    alpha1 = 0.d0
    alpha2 = 0.d0
!   alpha provisoire pour selection entre F(alpha) = 0 ou G(alpha) = 0
    alpha_pro = alpha
!   alpha prediction pour la boucle while locale F(alpha) = 0 ou G(alpha) = 0
    alpha_pred = alpha
!
! INITIALISATION DE LA CONTRAINTE
    call r8inir(6, 0.d0, sigma, 1)
!
! CALCUL DE LA CONTRAINTE HYDRO : DEBIT (LOI CUBIQUE)
! RECUPERATION DE LA PRESSION AU POINT DE GAUSS
!----------------------------------------------------
    if (ifhyme) then
        do n = 1, ndim-1
            sigma(ndim+n) = -rhof*gp(n)*(max(amin, delta(1)+amin))**3/(12*visf)
        end do
    end if
!
! CALCUL DE LA CONTRAINTE MECANIQUE
!----------------------------------
! CONTRAINTE DE PREDICTION ELASTIQUE
!
! CONTRAINTE NORMALE MECANIQUE
    sigma(1) = kn*(delta(1)-plasti(1))
!
! CONTRAINTE TANGENTIELLE
    do i = 2, ndim
        sigma(i) = kt*(delta(i)-plasti(i))
    end do
!
! CALCUL DE LA FORCE THERMODINAMIQUE PLASTIQUE X
! PHASE ELASTIQUE INITIALE A ALPHA=0 xvec=sigma
! PHASE D'EVOLUTION APRES L'ENDO INITIAL
!-----------------------------------------------
    if (alpha .ge. alpha_min) then
        S_alpha = ((1.d0-alpha)**m1)/(alpha**m2)
    end if
    xvec(1) = sigma(1)-Bn*S_alpha*plasti(1)
    do i = 2, ndim
        xvec(i) = sigma(i)-Bt*S_alpha*plasti(i)
    end do
!  Calcul de la norme de la partie tangente de X
    xnormt = 0.d0
    do i = 2, ndim
        xnormt = xnormt+xvec(i)**2
    end do
    xnormt = sqrt(xnormt)
!  Critere dans la prediction elastique
    criter = xnormt+mu*xvec(1)-cbar
    if (criter .le. 0.d0) then
!       La solution elastique est la solution
        ifplas = 0
        dlam = 0.d0
        do i = 1, ndim
            dplasti(i) = 0.
        end do
    else
!       Il y a endommagement, il faut calculer alpha
        ifplas = 1
!       Resolution de F(alpha) = 0 pour trouver la valeur de alpha
        if (alpha .le. alpha_min) then
            alpha1 = alpha_min
        else
            alpha1 = alpha
        end if
        alpha2 = 1.d0
        call finlf1(ndim, delta, plasti, alpha1, kn, &
                    kt, mu, Bn, Bt, m1, &
                    m2, cbar, D1, res1)
        call finlf1(ndim, delta, plasti, alpha2, kn, &
                    kt, mu, Bn, Bt, m1, &
                    m2, cbar, D1, res2)
!       Alpha ne bouge pas si la solution est hors intervalle (probleme de presicion num)
        if (res1*res2 .gt. 0.) then
            alpha_pro = alpha
!       Activation de la recherche de alpha pour F(alpha)=0 par methode mixte
!       On ne rentre pas dans le while si F(alpha-) < crit_conv
        else
!           METHODE DE BISECTION POUR TROUVER LA VALEUR DE ALPHA
            methode = 'corde'
            res_inter = res1
            alpha_pred = alpha
            do while (abs(alpha2-alpha1) .gt. alpha_min .and. abs(res_inter) .gt. crit_conv)
!               Recherche de alpha par methode mixte (dichotomie + corde)
                alpha_corde = (alpha1*res2-alpha2*res1)/(res2-res1)
                alpha_dicho = (alpha1+alpha2)/2.
                call finlf1(ndim, delta, plasti, alpha_corde, kn, &
                            kt, mu, Bn, Bt, m1, &
                            m2, cbar, D1, res_corde)
                call finlf1(ndim, delta, plasti, alpha_dicho, kn, &
                            kt, mu, Bn, Bt, m1, &
                            m2, cbar, D1, res_dicho)
                if (methode .eq. 'corde') then
                    res_inter = res_dicho
                    alpha_pred = alpha_dicho
                    methode = 'dicho'
                else
                    res_inter = res_corde
                    alpha_pred = alpha_corde
                    methode = 'corde'
                end if
                if (res_inter .gt. 0.) then
                    alpha1 = alpha_pred
                    res1 = res_inter
                else
                    alpha2 = alpha_pred
                    res2 = res_inter
                end if
            end do
!           Valeur provisoire de alpha
            alpha_pro = alpha_pred
        end if
!
!       CONSTRUCTION DU VECTEUR X(delta+,pl-,alpha+)
        if (alpha_pro .le. alpha_min) then
            ifplas = 0
            dlam = 0.
            do i = 1, ndim
                dplasti(i) = 0.
            end do
        else
!           Toutes les équations sont regularisees par alpha**m2
            xvecpl(1) = sigma(1)*alpha_pro**m2-Bn*plasti(1)*(1-alpha_pro)**m1-cbar*alpha**m2/mu
            do i = 2, ndim
                xvecpl(i) = sigma(i)*alpha_pro**m2-Bt*plasti(i)*(1-alpha_pro)**m1
            end do
            xnormp = 0.d0
            do i = 2, ndim
                xnormp = xnormp+xvecpl(i)**2
            end do
            xnormp = sqrt(xnormp)
!           Dlambda de l'equation 2.63
            dlam = xnormp+mu*xvecpl(1)
            dlam = dlam/((mu**2*kn+kt)*alpha_pro**m2+(mu**2*Bn+Bt)*(1-alpha_pro)**m1)
!           Dplasti_n de l'equation 2.54a
            dplasti(1) = mu*dlam
!           Xn de l'equation 2.51c multiplie par alpha**m2
            xnprov = xvecpl(1)-(kn*alpha_pro**m2+Bn*(1-alpha_pro)**m1)*dplasti(1)
!
            if (xnprov .le. cbar*alpha_pro**m2/mu) then
!               La solution a été trouvée
                alpha = alpha_pro
                do i = 2, ndim
                    dplasti(i) = dlam*xvecpl(i)/xnormp
                end do
            else
!               METHODE DE BISECTION POUR TROUVER LA VALEUR DE ALPHA
!               Même logique de resolution que dans F(alpha)=0
!               Selection de la matrice tangente 2 de la fonction G(alpha)
                ifplas = 2
                if (alpha .le. alpha_min) then
                    alpha1 = alpha_min
                else
                    alpha1 = alpha
                end if
                alpha2 = 1.d0
                call finlf2(ndim, delta, alpha1, kn, kt, &
                            mu, Bn, Bt, m1, m2, &
                            cbar, D1, res1)
                call finlf2(ndim, delta, alpha2, kn, kt, &
                            mu, Bn, Bt, m1, m2, &
                            cbar, D1, res2)
                if (res1*res2 .gt. 0.) then
                    alpha_pred = alpha
                else
                    methode = 'corde'
                    res_inter = 1.d0
                    alpha_pred = alpha
                    do while (abs(alpha2-alpha1) .gt. alpha_min .and. abs(res_inter) .gt. alpha_min)
!                       Recherche de alpha par methode mixte (dichotomie + corde)
                        alpha_corde = (alpha1*res2-alpha2*res1)/(res2-res1)
                        alpha_dicho = (alpha1+alpha2)/2.
                        call finlf2(ndim, delta, alpha_corde, kn, kt, &
                                    mu, Bn, Bt, m1, m2, &
                                    cbar, D1, res_corde)
                        call finlf2(ndim, delta, alpha_dicho, kn, kt, &
                                    mu, Bn, Bt, m1, m2, &
                                    cbar, D1, res_dicho)
                        if (methode .eq. 'corde') then
                            res_inter = res_dicho
                            alpha_pred = alpha_dicho
                            methode = 'dicho'
                        else
                            res_inter = res_corde
                            alpha_pred = alpha_corde
                            methode = 'corde'
                        end if
                        if (res_inter .gt. 0.) then
                            alpha1 = alpha_pred
                            res1 = res_inter
                        else
                            alpha2 = alpha_pred
                            res2 = res_inter
                        end if
                    end do
                end if
                alpha = alpha_pred
!
                if (alpha .le. alpha_min) then
                    ifplas = 0
                    dlam = 0.
                    do i = 1, ndim
                        dplasti(i) = 0.
                    end do
                else
                    xvecpl(1) = sigma(1)*alpha**m2-Bn*plasti(1)*(1-alpha)**m1-cbar*alpha**m2/mu
                    do i = 2, ndim
                        xvecpl(i) = sigma(i)*alpha**m2-Bn*plasti(i)*(1-alpha)**m1
                    end do
                    dplasti(1) = xvecpl(1)/(kn*alpha**m2+Bn*(1-alpha)**m1)
                    do i = 2, ndim
                        dplasti(i) = xvecpl(i)/(kt*alpha**m2+Bt*(1-alpha)**m1)
                    end do
                    dlam = dplasti(1)/mu
                end if
            end if
        end if
!
!       Actualisation des contraintes
        sigma(1) = sigma(1)-kn*dplasti(1)
        do i = 2, ndim
            sigma(i) = sigma(i)-kt*dplasti(i)
        end do
    end if
!
! PRISE EN COMPTE DE LA PRESSION DE FLUIDE EVENTUELLE
! PRESFL : IMPOSEE, PRESG : CALCULEE (MODELISATION HYME)
    if (ifhyme) then
        sigma(1) = sigma(1)-presg
    else
        sigma(1) = sigma(1)-presfl
    end if
!
!
! ACTUALISATION DES VARIABLES INTERNES
!-------------------------------------
! V1 : INDICATEUR DE PLASTICITE CUMULEE AU COURS DU GLISSEMENT
! V2 : INDICATEUR DE GLISSEMENT (0 : NON, 1 : OUI F, 2 : OUI G)
! V3 : DEPLACEMENT NORM PLASTIQUE PAR RAPPORT AU POINT DE DEPART
! V4-V5 : VECTEUR DE DEPLACEMENT TANG PLASTIQUE PAR RAPPORT AU POINT DE DEPART
!            (INDIQUE LA POSITION D'EQUILIBRE ACTUELLE)
! V6 : VARIABLE D'ENDOMMAGEMENT
! V7 A V9 : VALEURS DU DEPLACEMENT DANS LE REPERE LOCAL
! V10 : PAS UTILISE
! V11 : CONTRAINTE MECANIQUE NORMALE (SANS INFLUENCE PRESSION DE FLUIDE)
! V12 A V14 : COMPOSANTES DU GRADIENT DE PRESSION DANS LE REPERE GLOBAL
! V15 A V17 : COMPOSANTES DU FLUX HYDRO DANS LE REPERE GLOBAL
! V18 : PRESSION DE FLUIDE IMPOSEE OU CALCULEE ET INTERPOLEE (EN HYME)
! V19-V20 : CONTRAINTES MECANIQUE TANGENTIELLES
!
    vip(1) = vim(1)+dlam
    vip(2) = ifplas
    do i = 1, ndim
        vip(i+2) = vim(i+2)+dplasti(i)
    end do
    vip(6) = alpha
    vip(7) = delta(1)
    vip(8) = delta(2)
    if (ndim .eq. 3) then
        vip(9) = delta(3)
    else
        vip(9) = 0.d0
    end if
    vip(10) = 0.d0
!
!   FLUX, GRAD DE PRESSION ET PRESSION DANS LE REPERE GLOBAL
    if (ifhyme) then
        gploc(1) = 0.d0
        gploc(2) = gp(1)
        if (ndim .eq. 3) then
            gploc(3) = gp(2)
        end if
!
        fhloc(1) = 0.d0
        fhloc(2) = sigma(ndim+1)
        if (ndim .eq. 3) then
            fhloc(3) = sigma(2*ndim-1)
        end if
!
        call matinv('S', ndim, coorot(ndim+1), invrot, r8bid)
        call pmavec('ZERO', ndim, invrot, gploc, gpglo)
        call pmavec('ZERO', ndim, invrot, fhloc, fhglo)
!
!       CONTRAINTE MECANIQUE NORMALE SANS PRESSION DE FLUIDE CALCULEE
!       ON ANNULE SON INFLUENCE
        vip(11) = sigma(1)+presg
        vip(12) = gpglo(1)
        vip(13) = gpglo(2)
        vip(15) = fhglo(1)
        vip(16) = fhglo(2)
        if (ndim .eq. 3) then
            vip(14) = gpglo(3)
            vip(17) = fhglo(3)
        else
            vip(14) = 0.d0
            vip(17) = 0.d0
        end if
!
!       PRESSION DE FLUIDE CALCULEE AUX NOEUDS (DDL) ET INTERPOL AU PG
        vip(18) = presg
    else
!       CONTRAINTE MECANIQUE NORMALE SANS PRESSION DE FLUIDE IMPOSEE
!       ON ANNULE SON INFLUENCE
        vip(11) = sigma(1)+presfl
!       VI PAS UTILISEES EN MODELISATION NON HYME
        vip(12) = 0.d0
        vip(13) = 0.d0
        vip(14) = 0.d0
        vip(15) = 0.d0
        vip(16) = 0.d0
        vip(17) = 0.d0
!       PRESSION DE FLUIDE IMPOSEE AU PG :
        vip(18) = presfl
    end if
    vip(19) = sigma(2)
    if (ndim .eq. 3) then
        vip(20) = sigma(3)
    else
        vip(20) = 0.d0
    end if
!
500 continue
!
! INITIALISATION DE LA MATRICE TANGENTE
    call r8inir(6*6, 0.d0, dsidep, 1)
!
    if (.not. rigi) goto 999
!
! CALCUL DE LA MATRICE TANGENTE HYDRO
!------------------------------------
    if (ifhyme) then
!       TERME : DW/DGP  (POUR KTAN P P)
        do n = 1, ndim-1
            dsidep(ndim+n, ndim+n) = -rhof*(max(amin, delta(1)+amin))**3/(12*visf)
        end do
!       TERME : DW/DDELTA_N  (POUR KTAN P U)
        do n = 1, ndim-1
            if (delta(1) .lt. 0.d0) then
                dsidep(ndim+n, 1) = 0.d0
            else
                dsidep(ndim+n, 1) = -3*rhof*gp(n)*(delta(1)+amin)**2/(12*visf)
            end if
        end do
    end if
!
! CALCUL DE LA MATRICE TANGENTE MECA
!------------------------------------
    if (ifplas .eq. 0 .or. alpha .lt. alpha_min*10.d0) then
!     Utilisation de la matrice elastique si :
!       - alpha est plus petit que alpha_min*10 (transition douce)
!       - ifplas = 0 cas de la decharge
!     DSIGMA_N/DDELTA_N
        dsidep(1, 1) = kn
!     DSIGMA_N/DDELTA_T
        do i = 2, ndim
            dsidep(1, i) = 0.d0
        end do
!     DSIGMA_T/DDELTA_N
        do i = 2, ndim
            dsidep(i, 1) = 0.d0
        end do
!     DSIGMA_T/DDELTA_T
        do j = 2, ndim
            do i = j, ndim
                if (i .eq. j) then
                    kronec = 1
                else
                    kronec = 0
                end if
                dsidep(j, i) = kt*kronec
                dsidep(i, j) = kt*kronec
            end do
        end do
!
    else
!     Utilisation de la matrice tangente A ou B :
!       - Matrice tangente 1 : si resolution de F(alpha) = 0
!       - Matrice tangente 2 : si resolution de G(alpha) = 0
!
!       Calcul des termes communs aux deux matrices
        An = Bn*((1.d0-alpha)**m1)/(alpha**m2)
        Apn = Bn*((m2-m1)*alpha-m2)*((1-alpha)**(m1-1))/(alpha**(m2+1))
        A2pn = Bn*( &
               (m1-m2)*(m1-m2-1.)*alpha**2+2.*m2*(m1-m2-1.)*alpha+m2*(m2+1.))*(1.-alpha**(m1-2))/&
               &alpha**(m2+2 &
               )
        At = Bt*((1.d0-alpha)**m1)/(alpha**m2)
        Apt = Bt*((m2-m1)*alpha-m2)*((1-alpha)**(m1-1))/(alpha**(m2+1))
        A2pt = Bt*( &
               (m1-m2)*(m1-m2-1.)*alpha**2+2.*m2*(m1-m2-1.)*alpha+m2*(m2+1.))*(1.-alpha**(m1-2))/&
               &alpha**(m2+2 &
               )
!
        if (ifplas .eq. 1) then
!           Utilisation de la matrice tangente 1 de la fonction F(alpha)
!           Initialisation des vecteurs
            do i = 1, 2
                W_vect(i) = 0.d0
                dfddt(i) = 0.d0
                ddlddt(i) = 0.d0
            end do
!
!           Calcul des expressions intermediaires de la matrice tangente 1
!           Expression du vecteur W (equation 2.55)
            do i = 2, ndim
                W_vect(i-1) = kt*(delta(i)-plasti(i))-At*plasti(i)
            end do
!           Calcul de la norme de W
            W_norm = sqrt(DOT_PRODUCT(W_vect(:ndim-1), W_vect(:ndim-1)))
            if (W_norm .eq. 0.d0) then
                W_norm = 1.d0
            end if
!           Calcul de la norme de plasti_t selon ndim
            pt_norm = sqrt(DOT_PRODUCT(plasti(2:ndim), plasti(2:ndim)))
!           Calcul du produit scalaire plasti*W_vect
            ptxw = DOT_PRODUCT(plasti(2:ndim), W_vect(:ndim-1))
!           Calcul du discriminant de l'equation 2.60
            disc = ( &
                   mu*Apn*plasti(1)+Apt*ptxw/W_norm)**2-(mu**2*Apn+Apt)*(Apn*plasti(1)**2+Apt*pt_&
                   &norm**2+2.d0*D1 &
                   )
!           Calcul de la derivee de delta lambda par rapport a alpha
            ddldalp = ( &
                      -( &
                      Apt*ptxw/W_norm+mu*Apn*plasti(1))*(mu**2*kn+kt+mu**2*An+At)-(W_norm+mu*kn*(&
                      &delta(1)-plasti(1))-mu*An*plasti(1)-cbar)*(mu**2*Apn+Apt) &
                      )/(mu**2*kn+kt+mu**2*An+At &
                      )**2
!           Calcul de la derivee de F par rapport a delta_n
            dfddn = mu*kn*(mu**2*Apn+Apt)
!           Calcul de la derivee de F par rapport a delta_t
            do i = 2, ndim
                dfddt(i-1) = ( &
                             kt*W_vect(i-1)/W_norm)*(mu**2*Apn+Apt)+(1.d0+(mu*Apn*plasti(1)+Apt*p&
                             &txw/W_norm)/sqrt(disc))*(mu**2*kn+kt+mu**2*An+At)*Apt*kt/W_norm*(p&
                             &lasti(i)-DOT_PRODUCT(plasti(2:ndim), W_vect(:ndim-1))*W_vect(i-1)/W&
                             &_norm**2 &
                             )
            end do
!           Calcul de la derivee de F par rapport a alpha
            dfdalp = ( &
                     W_norm+mu*kn*( &
                     delta(1)-plasti(1))-mu*An*plasti(1))*(mu**2*A2pn+A2Pt)+((1.d0+(mu*Apn*plasti&
                     &(1)+Apt*ptxw/W_norm**mu)/sqrt(disc))*(mu*A2pn*plasti(1)+A2pt*ptxw/W_norm-Ap&
                     &t**2/W_norm*(DOT_PRODUCT(plasti(2:ndim), plasti(2:ndim))-DOT_PRODUCT(plasti&
                     &(2:ndim), W_vect(:ndim-1))**2/W_norm**2))-((mu*2*A2pn+A2pt)*(Apn*plasti(1)*&
                     &*2+Apt*pt_norm**2+2.d0*D1)+(mu**2*Apn+Apt)*(A2pn*plasti(1)**2+A2pt*pt_norm*&
                     &*2))/(2.d0*sqrt(disc)))*(kt+mu**2*kn+At+mu**2*An)+sqrt(disc &
                     )*(Apt+mu**2*Apn &
                     )
!           Calcul de la derivee de delta lambda par rapport a delta_n
            ddlddn = (mu*kn/(kt+mu**2*kn+At+mu**2*An))-ddldalp*dfddn*(1.d0/dfdalp)
!           Calcul de la derivee de delta lambda par rapport a delta_t
            do i = 2, ndim
                ddlddt(i-1) = W_vect(i-1)*kt/((kt+mu**2*kn+At+mu**2*An)*W_norm)-ddldalp*dfddt(i-1&
                              &)*(1.d0/dfdalp)
            end do
!           Expression de delta lambda
            d_lambda = ( &
                       W_norm+mu*kn*(delta(1)-plasti(1))-mu*An*plasti(1)-cbar &
                       )/(mu**2*kn+kt+mu**2*An+At &
                          )
!
!           Ecriture des termes de la matrice tangente 1
!           DSIGMA_N/DDELTA_N
            dsidep(1, 1) = kn*(1.d0-mu*ddlddn)
!           DSIGMA_N/DDELTA_T
            do i = 2, ndim
                dsidep(1, i) = -mu*kn*ddlddt(i-1)
            end do
!           DSIGMA_T/DDELTA_N
            do i = 2, ndim
                dsidep(i, 1) = -kt*( &
                               ddlddn*W_vect(i-1)/W_norm+d_lambda*Apt*dfddn/(W_norm*dfdalp)*(plas&
                               &ti(i)-DOT_PRODUCT(plasti(2:ndim), W_vect(:ndim-1))*W_vect(i-1)/W_&
                               &norm**2) &
                               )
            end do
!           DSIGMA_T/DDELTA_T
            do j = 2, ndim
                do i = 2, ndim
                    if (i .eq. j) then
                        identity = 1
                    else
                        identity = 0
                    end if
                    dsidep(i, j) = kt*( &
                                   identity-ddlddt(i-1)*W_vect(j-1)/W_norm-d_lambda/W_norm*(ident&
                                   &ity*kt+plasti(i)*dfddt(j-1)*Apt/dfdalp-W_vect(i-1)*W_vect(j-1&
                                   &)*kt/W_norm**2-W_vect(i-1)*dfddt(j-1)*DOT_PRODUCT(plasti(2:nd&
                                   &im), W_vect(:ndim-1))*Apt/(dfdalp*W_norm**2)) &
                                   )
                end do
            end do
!
        else if (ifplas .eq. 2) then
!           Utilisation de la matrice tangente 2 de la fonction G(alpha)
!           Initialisation des vecteurs
            do i = 1, 2
                dgddt(i) = 0.d0
                ddpnddt(i) = 0.d0
                ddptddn(i) = 0.d0
            end do
            do i = 1, 2
                do j = 1, 2
                    ddptddt(i, j) = 0.d0
                end do
            end do
!
!           Calcul des expressions intermediaires de la matrice tangente 2
!           Calcul de la derivee de G par rapport a delta_n
            dgddn = 2.d0*Apn*kn*(mu*kn*delta(1)-cbar)/(mu*(kn+An)**2)
!           Calcul de la derivee de G par rapport a delta_t
            do i = 2, ndim
                dgddt(i-1) = 2.d0*Apt*kt**2*delta(i)/((kt+At)**2)
            end do
!           Calcul de la derivee de G par rapport a alpha
            dgdalp = ( &
                     A2pn*(kn+An)-2.d0*Apn**2)*(mu*kn*delta(1)-cbar)**2/(mu**2*(kn+An)**3)+(A2pt*&
                     &(kt+At)-2.d0*Apt**2)*(kt*sqrt(DOT_PRODUCT(plasti(2:ndim), plasti(2:ndim))))&
                     &**2/((kt+At)**3 &
                     )
!           Calcul de la derivee de delta_pn par rapport a delta_n
            ddpnddn = kn/(kn+An)+(mu*kn*delta(1)-cbar)*Apn/(mu*(kn+An)**2)*dgddn/dgdalp
!           Calcul de la derivee de delta_pn par rapport a delta_t
            do i = 2, ndim
                ddpnddt(i-1) = ((mu*kn*delta(1)-cbar)*Apn)/(mu*(kn+An)**2)*dgddt(i-1)/dgdalp
            end do
!           Calcul de la derivee de delta_pt par rapport a delta_n
            do i = 2, ndim
                ddptddn(i-1) = kt*Apt*delta(i)/((kt+At)**2)*dgddn/dgdalp
            end do
!           Calcul de la derivee de delta_pt par rapport a delta_t
            do i = 2, ndim
                do j = 2, ndim
                    ddptddt(i-1, j-1) = kt*Apt*delta(i)*dgddt(j-1)/(dgdalp*(kt+At)**2)
                end do
            end do
!
!           Ecriture des termes de la matrice tangente 2
!           DSIGMA_N/DDELTA_N
            dsidep(1, 1) = kn*(1.d0-ddpnddn)
!           DSIGMA_N/DDELTA_T
            do i = 2, ndim
                dsidep(1, i) = -kn*ddpnddt(i-1)
            end do
!           DSIGMA_T/DDELTA_N
            do i = 2, ndim
                dsidep(i, 1) = -kt*ddptddn(i-1)
            end do
!           DSIGMA_T/DDELTA_T
            do j = 2, ndim
                do i = 2, ndim
                    if (i .eq. j) then
                        identity = 1
                    else
                        identity = 0
                    end if
                    dsidep(i, j) = kt*(identity-(identity*kt/(kt+At)+ddptddt(i-1, j-1)))
                end do
            end do
        end if
    end if
!
999 continue
end subroutine
