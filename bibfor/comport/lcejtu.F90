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
subroutine lcejtu(BEHInteg, &
                  ndim, option, epsm, deps, sigm, &
                  sigp, dsidep, vim, vip, typmod, &
                  instam, instap)
!
    use Behaviour_type
    implicit none
!
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    integer(kind=8), intent(in) :: ndim
    real(kind=8), intent(in) :: epsm(6), deps(6), sigm(6), vim(*)
    real(kind=8), intent(in) :: instam, instap
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: vip(*), sigp(6), dsidep(6, 6)
!-----------------------------------------------------------------------
!     LOI DE COMPORTEMENT CZM_TURON D'UN JOINT ANISOTROPE
!     COMPORTEMENTS NORMAL ET TANGENTIEL DIFFERENTS
!     CRITERE D'INITIATION DE L'ENDO DE TYPE BENZEGGAGH-KENANE (BK) OU YE
!     CRITERE DE PROPAGATION DE LA FISSURE DE TYPE BK
!     MODELE DE TURON 2006 (university of Girona, Spain)
!     Référence : Turon et al. / Mechanics of Materials 38 (2006) 1072-1089
!
! IN : EPSM - SAUT INSTANT MOINS
! IN : DEPS - INCREMENT DE SAUT
! IN : IMATE, OPTION, VIM, COOROT, INSTAM, INSTAP
! IN : SIGM - CONTRAINTES A L'INSTANT MOINS
! OUT : SIGMA , DSIDEP , VIP
!-----------------------------------------------------------------------
    integer(kind=8) :: nbpa
    parameter(nbpa=8)
    integer(kind=8) :: cod(nbpa)
    character(len=16) :: nom(nbpa)
    real(kind=8) :: val(nbpa)
    integer(kind=8) :: i, j, diss, cass
    real(kind=8) :: inst, delta(3), ddelta(3)
    real(kind=8) :: k, c, eta, crit
    real(kind=8) :: delta_N_0, delta_T_0
    real(kind=8) :: delta_N_f, delta_T_f
    real(kind=8) :: delta_N_pos, delta_T, lambda
    real(kind=8) :: beta, b, t, delta_0, delta_f, r, d, g
    real(kind=8) :: delta_f_N, delta_f_T
    real(kind=8) :: quot, a(6, 6), aa(6, 6)
    real(kind=8) :: zero, pi
    character(len=16) :: type_comp
    character(len=1) :: poum
    aster_logical :: resi, rigi, elas, pred, l_lambda0
    blas_int :: b_incx, b_incy, b_n
!
    data nom/'K', 'SIGM_C_N', 'SIGM_C_T', 'GC_N', 'GC_T', 'C_RUPT', &
        'ETA_BK', 'CRIT_INIT'/
!
!-----------------------------------------------------------------------
! RQ SUR LES NOTATIONS :
!   - delta_N_0 et delta_N_f : sauts norm pour init et prop en mode N pur
!   - delta_T_0 et delta_T_f : sauts tang pour init et prop en mode T pur
!   - delta_f_N et delta_f_T : sauts norm et tang lorsque le saut équivalent
!     (lambda) atteint son seuil de prop delta_f en mode mixte
!-----------------------------------------------------------------------
!
    ASSERT(typmod(2) .eq. 'ELEMJOIN')
!
! INITIALISATIONS
    type_comp = 'COMP_ELAS'
    zero = r8prem()
    pi = r8pi()
    resi = .false.
    rigi = .false.
    elas = .false.
    l_lambda0 = .false.
    inst = 0.d0
    delta = 0.d0
    ddelta = 0.d0
    delta_N_pos = 0.d0
    delta_T = 0.d0
    delta_0 = 0.d0
    delta_f = 0.d0
    delta_f_N = 0.d0
    delta_f_T = 0.d0
    diss = 0.d0
    cass = 0.d0
    sigp = 0.d0
    dsidep = 0.d0
    quot = 0.d0
    a = 0.d0
    aa = 0.d0
!
! OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
    rigi = option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA'
    elas = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'RIGI_MECA_ELAS'
    pred = option .eq. 'RIGI_MECA_TANG'
!
! INSTANT DE CALCUL T- (INSTANT PRECEDENT) OU T+ (INSTANT ACTUEL)
    inst = instam
    if (resi) inst = instap
!
! SAUT DE DEPLACEMENT A L'INSTANT ACTUEL
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm, b_incx, delta, b_incy)
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, deps, b_incx, ddelta, b_incy)
    if (resi) then
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, ddelta, b_incx, delta, &
                   b_incy)
    end if
!
! RECUPERATION DES PARAMETRES PHYSIQUES
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if
!
    call rcvalb(BEHInteg%materPara%schemePara%fami, &
                BEHInteg%materPara%schemePara%kpg, &
                BEHInteg%materPara%schemePara%ksp, &
                poum, &
                BEHInteg%materPara%jvMaterCode, &
                ' ', 'RUPT_TURON', &
                0, ' ', [0.d0], &
                nbpa, nom, val, cod, 2)
!
!   * VERIFICATION
    if (val(1)*val(2)*val(3)*val(4)*val(5) .lt. r8prem()) then
        call utmess('F', 'COMPOR4_74', nk=5, valk=nom(1:5))
    end if
!
!   * RIGIDITE DE PENALISATION (IDENTIQUE POUR TOUS LES MODES)
    k = val(1)
!
!   * CALCUL DES SEUILS D'INITIATION DE L'ENDOMMAGEMENT POUR LES MODES PURS
    delta_N_0 = val(2)/k
    delta_T_0 = val(3)/k
! METTRE UN MESSAGE D'ERREUR SI LES VALEURS SONT NULLES
!
!   * CALCUL DES SEUILS DE PROPAGATION DE LA FISSURE POUR LES MODES PURS
    delta_N_f = 2*val(4)/(k*delta_N_0)
    delta_T_f = 2*val(5)/(k*delta_T_0)
!
!   * COEFFICIENT DE RIGIDITE POST-RUPTURE
    c = val(6)
!
!   * PARAMETRE PUISSANCE DE LA FORMULE DE BENZEGGAGH ET KENANE
    eta = val(7)
!
!   * CRITERE D'INITIATION DE L'ENDOMMAGEMENT
!     CRIT = 0 : TURON, CRIT = 1 :YE
    crit = val(8)
!
! PARTIE POSITIVE DU SAUT NORMAL
    delta_N_pos = max(0.d0, delta(1))
!
! SAUT EQUIVALENT TANGENTIEL
    delta_T = sqrt(delta(2)**2+delta(3)**2)
!
! SAUT EQUIVALENT TOTAL
    lambda = sqrt(delta_N_pos**2+delta_T**2)
    l_lambda0 = (lambda .lt. r8prem())
!
    if (.not. l_lambda0) then
! TAUX DE MIXITE A T+ (BETA A PARTIR DES SAUTS, B A PARTIR DES TAUX DE RESTIT G)
        beta = delta_T/(delta_T+delta_N_pos)
        b = beta**2/(1-2*beta+2*beta**2)
! SEUILS D'INITIATION DE L'ENDOMMAGEMENT A T+
! * CRITERE DE TURON (DE TYPE BK, FONCTION DU TAUX DE MIXITE)
        if (crit .eq. 0) then
            delta_0 = sqrt(delta_N_0**2+(delta_T_0**2-delta_N_0**2)*b**eta)
            ASSERT(delta_0 .gt. r8prem())
! * CRITERE DE YE (DE TYPE ELLIPTIQUE)
        else if (crit .eq. 1) then
            if (delta_N_pos .lt. r8prem()) then
                t = pi/2
            else
                t = atan(delta_N_0/delta_T_0*delta_T/delta_N_pos)
            end if
            delta_0 = sqrt((delta_N_0*cos(t))**2+(delta_T_0*sin(t))**2)
            ASSERT(delta_0 .gt. r8prem())
        else
            ASSERT(.False.)
        end if
! SEUIL DE PROPAGATION DE LA FISSURE A T+ (DE TYPE BK)
        delta_f = 1/delta_0*(delta_N_0*delta_N_f+(delta_T_0*delta_T_f-delta_N_0*delta_N_f)*b**eta &
                             )
        ASSERT(delta_f .gt. r8prem())
    else
! CAS PARTICULIERS SI ON (RE)PASSE PAR UN ETAT DE SAUT NUL :
! ON RESTE AVEC LES VALEURS PRECEDENTES
        beta = vim(13)
        b = vim(14)
        delta_0 = vim(15)
        delta_f = vim(16)
! POUR LE PREMIER PAS DE TEMPS SANS ETAT INIT
        if (delta_0 .lt. r8prem()) then
            call utmess('A', 'COMPOR4_73')
            delta_0 = min(delta_N_0, delta_T_0)
            delta_f = min(delta_N_f, delta_T_f)
        end if
    end if
!
! ESTIMATION DE LA VARIABLE SEUIL ACTUELLE
    r = max(delta_0, vim(1))
!
! INITIALISATION DES INDICATEURS DE DISS ET D'ENDO POUR RIGI_MECA_TANG (+ SECANTE PENALISEE)
    if (.not. resi) then
        if (elas) then
            diss = 0
        else
            diss = nint(vim(4))
        end if
        cass = nint(vim(5))
    end if
!
! CALCUL DE LA CONTRAINTE
    if (resi) then
        if (((lambda-delta_f) .ge. zero) .or. ((r-delta_f) .ge. zero)) then
            diss = 0
            cass = 2
            d = 1.0
            g = 1.0
        else
            if ((delta_0-lambda) .ge. zero) then
                g = lambda*lambda/delta_0/delta_f
            else
                g = 1.-(delta_f-lambda)*(delta_f-lambda)/(delta_f-delta_0)/delta_f
            end if
!
            if ((r-lambda) .gt. zero) then
                diss = 0
                if ((r-delta_0) .ge. zero) then
                    cass = 1
                else
                    cass = 0
                end if
                d = delta_f*(r-delta_0)/(r*(delta_f-delta_0))
!                 g = (r-delta_0)/(delta_f-delta_0)
!                 g = (delta_f-(delta_f-r)*(delta_f-r)/(delta_f-delta_0))/delta_f
            else
                diss = 1
                cass = 1
                if (.not. l_lambda0) then
                    d = delta_f*(lambda-delta_0)/(lambda*(delta_f-delta_0))
!                     g = (lambda-delta_0)/(delta_f-delta_0)
!                     g = (delta_f-(delta_f-lambda)*(delta_f-lambda)/(delta_f-delta_0))/delta_f
                else
                    d = vim(3)
!                     g = vim(11)
                end if
            end if
        end if
        ASSERT((d-1.d0) .le. r8prem())
        ASSERT((g-1.d0) .le. r8prem())
!
        if (type_comp .eq. 'COMP_ELAS') then
            sigp(1) = -d*k*max(0.d0, -delta(1))+(1-d)*k*delta(1)
            do i = 2, ndim
                sigp(i) = (1-d)*k*delta(i)
            end do
        else if (type_comp .eq. 'COMP_INCR') then
            sigp(1) = sigm(1)-(d*k*max(0.d0, -delta(1))*delta(1)-vim(3)*k*max(0.d0, -(delta(1)-dd&
                      &elta(1)))*(delta(1)-ddelta(1)))+(1-d)*k*delta(1)-(1-vim(3))*k*(delta(1)-dd&
                      &elta(1))
            do i = 2, ndim
                sigp(i) = sigm(i)+(1-d)*k*delta(i)-(1-vim(3))*k*(delta(i)-ddelta(i))
            end do
        else
            ASSERT(.false.)
        end if
!         write(6,*) ''
!         write(6,*) 'sigm',sigm
!         write(6,*) 'sigp',sigp
!         write(6,*) 'd',d
!         write(6,*) 'k',k
!         write(6,*) 'ddelta',ddelta
!         write(6,*) 'delta', delta
!         write(6,*) 'epsm', epsm
!         write(6,*) 'deps', deps
!         write(6,*) 'option', option
!
!
! ACTUALISATION DES VARIABLES INTERNES
!-------------------------------------
! V1 : PLUS GRANDE NORME DU SAUT EQUIVALENT TOTAL
! V2 : VARIABLE SEUIL
! V3 : VARIABLE D'ENDOMMAGEMENT
! V4 : INDICATEUR DE DISSIPATION (0 : SI REGIME LIN, 1 : SI REGIME DISS)
! V5 : INDICATEUR D'ENDOMMAGEMENT (0 : SAIN, 1: ENDOMMAGE, 2: CASSE)
! V6 A V8 : VALEURS DU SAUT DANS LE REPERE LOCAL
! V9 : VALEUR DU SAUT EQUIVALENT TOTAL
! V10 : VALEUR DU SAUT EQUIVALENT TANGENTIEL
! V11 : POURCENTAGE D'ENERGIE DISSIPEE
! V12 : VALEUR DE L'ENERGIE DISSIPEE
! V13 : TAUX DE MIXITE BETA A T+ (CALCULE PAR LES SAUTS)
! V14 : TAUX DE MIXITE B A T+ (CALCULE PAR LES TAUX DE RESTITUTION D'ENERGIE)
! V15 : SEUIL D'INITIATION DE L'ENDOMMAGEMENT EN MODE MIXTE A T+
! V16 : SEUIL DE PROPAGATION DE LA FISSURE EN MODE MIXTE A T+
!
        vip(1) = max(lambda, vim(1))
        vip(2) = r
        vip(3) = d
        vip(4) = diss
        vip(5) = cass
        vip(6) = delta(1)
        vip(7) = delta(2)
        if (ndim .eq. 3) then
            vip(8) = delta(3)
        else
            vip(8) = 0.d0
        end if
        vip(9) = lambda
        vip(10) = delta_T
        vip(11) = g
        vip(12) = g*1/2*k*delta_0*delta_f
        vip(13) = beta
        vip(14) = b
        vip(15) = delta_0
        vip(16) = delta_f
    end if
!
! CALCUL DE LA MATRICE TANGENTE
! NB : EN TOUTE RIGUEUR, ELLE N'EST PAS DEFINIE EN DELTA(1)=0
! ON LA PROLONGE PAR SA VALEUR A GAUCHE.
    if (rigi) then
!
!   * RIGIDITE ARTIFICIELLE POST-RUPTURE
        if (cass .eq. 2) then
            if (delta(1) .gt. r8prem()) then
                dsidep(1, 1) = c*val(2)/delta_N_f
            else
                dsidep(1, 1) = k
            end if
            do i = 2, ndim
                dsidep(i, i) = c*val(3)/delta_T_f
            end do
        else
!
            if (abs(delta(1)) .gt. r8prem()) then
                quot = max(0.d0, -delta(1))/abs(delta(1))
!   * VAUT 1 si delta(1) < 0 et 0 si delta(1) > 0
            else
                quot = 1.d0
!   * VAUT 1 si delta(1) = 0 (prolongement à gauche)
            end if
!
            if ((diss .eq. 0) .or. elas) then
                d = delta_f*(r-delta_0)/(r*(delta_f-delta_0))
            else
                d = delta_f*(lambda-delta_0)/(lambda*(delta_f-delta_0))
            end if
!
            do i = 1, ndim
                a(i, i) = 1-d
            end do
            a(1, 1) = a(1, 1)+d*quot
            dsidep = a*k
!
            if (diss .eq. 1) then
                do i = 1, ndim
                    do j = 1, ndim
                        aa(i, j) = delta(i)*delta(j)
                    end do
                    aa(i, 1) = aa(i, 1)+delta(i)*max(0.d0, -delta(1))
                end do
                do j = 1, ndim
                    aa(1, j) = aa(1, j)+delta(j)*max(0.d0, -delta(1))
                end do
                aa(1, 1) = aa(1, 1)+(max(0.d0, -delta(1)))**2
                aa = aa*k*delta_f*delta_0/(delta_f-delta_0)*1/lambda**3
                dsidep = dsidep-aa
            end if
        end if
    end if
!
end subroutine
