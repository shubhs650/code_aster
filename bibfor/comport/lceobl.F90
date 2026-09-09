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
! aslint: disable=W1501,W0413
!
subroutine lceobl(ndim, imate, carcri, epsm, &
                  deps, vim, option, sigm, vip, &
                  dsidep, iret)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dfmdf.h"
#include "asterfort/dfpdf.h"
#include "asterfort/diago3.h"
#include "asterfort/lceob1.h"
#include "asterfort/lceob2.h"
#include "asterfort/lceob3.h"
#include "asterfort/lceobb.h"
#include "asterfort/meobl1.h"
#include "asterfort/meobl2.h"
#include "asterfort/meobl3.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvala.h"
#include "asterfort/sigeob.h"
!
    character(len=16) :: option
    integer(kind=8) :: ndim, imate, iret
    real(kind=8) :: epsm(6), deps(6), vim(7)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: sigm(6), vip(7), dsidep(6, 6)
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT DU MODELE D'ENDOMMAGEMENT ANISOTROPE
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : NATURE DU MATERIAU
! IN  CRIT   : CRITERES DE CONVERGENCE LOCAUX
! IN  EPSM    : DEFORMATION EN T- REPERE GLOBAL
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  VIM     : VARIABLES INTERNES EN T-
! IN  OPTION  : OPTION DEMANDEE
!                 RIGI_MECA_TANG ->     DSIDEP
!                 FULL_MECA      -> SIG DSIDEP VIP
!                 RAPH_MECA      -> SIG        VIP
! OUT SIGM     : CONTRAINTE
! OUT VIP     : VARIABLES INTERNES
!                 1 A 6   -> TENSEUR D'ENDOMMAGEMENT DE TRACTION
!                 7       -> ENDOMMAGEMENT DE COMPRESSION
! OUT DSIDEP  : MATRICE TANGENTE DEFO
! OUT         : IRET CODE RETOUR
! ----------------------------------------------------------------------
!
    aster_logical :: rigi, resi, elas, reinit
    aster_logical :: total, tot1, tot2, tot3, dbloq
    integer(kind=8) :: ndimsi, i, j, k, l, p, q, m, n, t(3, 3)
    integer(kind=8) :: bdim, r1(6), r2(6)
!
    integer(kind=8) :: intmax
    real(kind=8) :: toler
!
    real(kind=8) :: treps, treb, eps(6), cc(6)
    real(kind=8) :: kron(6), mult, seuil, rk, rk1, rk2, un, deux
    real(kind=8) :: rac2, bobo(6, 6), zozo(6, 6), zaza(6, 6)
    real(kind=8) :: rtemp1, rtemp2, rtemp3, rtemp4, rtemp5, rtemp6
    real(kind=8) :: e, nu, alpha, lambda, mu, ecrob, ecrod
    real(kind=8) :: b(6), bm(6), br(6)
    real(kind=8) :: dm, d, vecb(3, 3), valb(3), tolb, trepsm
    real(kind=8) :: interm(3, 3), epi(6), bmr(6), epsr(6)
    real(kind=8) :: vecbr(3, 3), valbr(3), binter(6), ad
    real(kind=8) :: deltab(6), deltad, dsisup(6, 6)
    real(kind=8) :: interb(3, 3), epib(6)
    real(kind=8) :: dsiint(6, 6), dsimed(6, 6)
    real(kind=8) :: vecepm(3, 3), valepm(3)
!
    integer(kind=8) :: icodre(6)
    character(len=16) :: nomres(6)
    real(kind=8) :: valres(6)
!
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
    un = 1.d0
    deux = 2.d0
    t(1, 1) = 1
    t(1, 2) = 4
    t(1, 3) = 5
    t(2, 1) = 4
    t(2, 2) = 2
    t(2, 3) = 6
    t(3, 1) = 5
    t(3, 2) = 6
    t(3, 3) = 3
!
    r1(1) = 2
    r1(2) = 3
    r1(3) = 1
    r1(4) = 6
    r1(5) = 4
    r1(6) = 5
!
    r2(1) = 3
    r2(2) = 1
    r2(3) = 2
    r2(4) = 5
    r2(5) = 6
    r2(6) = 4
!
    rac2 = sqrt(deux)
    tolb = 1.d-2
    mult = 0.d0
!=====================================================================
!                            INITIALISATION
! ====================================================================
!
!---------------------------------------------------
! -- OPTION ET MODELISATION
!---------------------------------------------------
    rigi = (option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL')
    resi = (option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL')
    ndimsi = 2*ndim
    total = .false.
    reinit = .false.
!
!---------------------------------------------------
! -- LECTURE DES CARACTERISTIQUES THERMOELASTIQUES
!---------------------------------------------------
    nomres(1) = 'E'
    nomres(2) = 'NU'
    call rcvala(imate, ' ', 'ELAS', 0, ' ', &
                [0.d0], 2, nomres, valres, icodre, &
                1)
    e = valres(1)
    nu = valres(2)
    lambda = e*nu/(un+nu)/(un-deux*nu)
    mu = e/(deux*(un+nu))
!-------------------------------------------------
! -- LECTURE DES CARACTERISTIQUES D'ENDOMMAGEMENT
!-------------------------------------------------
    nomres(1) = 'ALPHA'
    nomres(2) = 'K0'
    nomres(3) = 'K1'
    nomres(4) = 'K2'
    nomres(5) = 'ECROB'
    nomres(6) = 'ECROD'
    call rcvala(imate, ' ', 'ENDO_ORTH_BETON', 0, ' ', &
                [0.d0], 6, nomres, valres, icodre, &
                1)
    alpha = valres(1)
    rk = valres(2)
    rk1 = valres(3)
    rk2 = valres(4)
    ecrob = valres(5)
    ecrod = valres(6)
!
    toler = carcri(3)
    intmax = int(carcri(1))
!
    elas = .true.
!
!---------------------------------------------
!---INITIALISATION MATRICE TANGENTE
!---------------------------------------------
    call r8inir(36, 0.d0, dsidep, 1)
    call r8inir(36, 0.d0, dsisup, 1)
!
!-------------------------------------------------
! -- DEFORMATIONS
!-------------------------------------------------
    if (resi) then
        do k = 1, ndimsi
            eps(k) = epsm(k)+deps(k)
        end do
    else
        do k = 1, ndimsi
            eps(k) = epsm(k)
        end do
        elas = .true.
    end if
!
    if (ndim .lt. 3) then
        do i = ndimsi+1, 6
            eps(i) = 0.d0
            epsm(i) = 0.d0
            deps(i) = 0.d0
        end do
    end if
!
    do i = 4, 6
        eps(i) = eps(i)/rac2
        epsm(i) = epsm(i)/rac2
        deps(i) = deps(i)/rac2
    end do
!-------------------------------------------------
! -- ENDOMMAGEMENT DANS LE REPERE GLOBAL
!-------------------------------------------------
    do i = 1, 3
        bm(i) = un-vim(i)
        b(i) = bm(i)
    end do
    do i = 4, 6
        bm(i) = -vim(i)
        b(i) = bm(i)
    end do
    dm = vim(7)
    d = dm
!
!------------------------------------------------------------------
!-- VERIFICATION SUR LES VALEURS PROPRES DE L ENDOMMAGEMENT
!-- DE TRACTION ET SUR L ENDOMMAGEMENT DE COMPRESSION POUR BLOQUAGE
!-- EVENTUEL DE L EVOLUTION
!------------------------------------------------------------------
    call diago3(b, vecb, valb)
    bdim = 3
    do i = 1, 3
        if (valb(i)-tolb .le. 0.d0) then
            bdim = bdim-1
        end if
    end do
    dbloq = .false.
    if (d-(un-tolb) .ge. 0.d0) then
        dbloq = .true.
    end if
!
!-------------------------------------------------
!--DEFINITION DU SEUIL----------------------------
!-------------------------------------------------
!
    call diago3(epsm, vecepm, valepm)
    trepsm = eps(1)+eps(2)+eps(3)
    if (trepsm .gt. 0.d0) then
        trepsm = 0.d0
    end if
    seuil = rk-rk1*trepsm*(atan2(-trepsm/rk2, un))
!
    if (resi) then
!----------------------------------------------------------------
!----CAS OU LES 3 VALEURS PROPRES D ENDO TRACTION SONT NON NULLES
!----------------------------------------------------------------
!
        if (bdim .eq. 3) then
            call lceob3(intmax, toler, eps, bm, dm, &
                        lambda, mu, alpha, ecrob, ecrod, &
                        seuil, bdim, b, d, mult, &
                        elas, dbloq, iret)
!--VERIF SUR ENDO FINAL POUR VOIR SI ENDO DEPASSE 1 OU PAS
!-- SI ENDO DEPASSE 1 PLUS QUE TOLERANCE ON PASSE DANS LCEOBB
!-- QUI DECOUPE L INCREMENT DE CHARGE POUR ALLER DOUCEMENT A ENDO=1
!
            call diago3(b, vecb, valb)
!
            do i = 1, 3
                if (valb(i) .lt. 0) then
                    reinit = .true.
                end if
                if (valb(i)-tolb .le. 0.d0) then
                    valb(i) = tolb-r8prem()
                    total = .true.
                end if
            end do
            if (d .gt. 1.d0) then
                reinit = .true.
            end if
            if (un-d-tolb .le. 0.d0) then
                d = un-tolb+r8prem()
                dbloq = .true.
            end if
!
            if (.not. reinit) then
                if (total) then
                    call r8inir(6, 0.d0, b, 1)
                    do i = 1, 3
                        do j = i, 3
                            do k = 1, 3
                                b(t(i, j)) = b(t(i, j))+vecb(i, k)*valb(k)* &
                                             vecb(j, k)
                            end do
                        end do
                    end do
                end if
            else
                call lceobb(intmax, toler, epsm, deps, bm, &
                            dm, lambda, mu, alpha, ecrob, &
                            ecrod, rk, rk1, rk2, b, &
                            d, mult, elas, dbloq, iret)
            end if
!
!----------------------------------------------------------------
!----CAS OU 1 VALEUR PROPRE EST NULLE----------------------------
!----------------------------------------------------------------
!
        else if (bdim .eq. 2) then
!
!-- ON RESTREINT L ESPACE CAR L ENDO N EVOLUE PLUS DANS 2 DIRECTIONS
!
            call r8inir(9, 0.d0, interm, 1)
            call r8inir(6, 0.d0, epi, 1)
            do i = 1, 3
                do l = 1, 3
                    do k = 1, 3
                        interm(i, l) = interm(i, l)+vecb(k, i)*eps(t(k, l))
                    end do
                    do j = i, 3
                        epi(t(i, j)) = epi(t(i, j))+interm(i, l)*vecb(l, j)
                    end do
                end do
            end do
            tot1 = .false.
            tot2 = .false.
            tot3 = .false.
            call r8inir(6, 0.d0, bmr, 1)
            if (valb(1)-tolb .le. 0.d0) then
                bmr(1) = valb(2)
                bmr(2) = valb(3)
                do i = 1, 6
                    epsr(i) = epi(r1(i))
                end do
                tot1 = .true.
            else if (valb(2)-tolb .le. 0.d0) then
                bmr(1) = valb(3)
                bmr(2) = valb(1)
                do i = 1, 6
                    epsr(i) = epi(r2(i))
                end do
                tot2 = .true.
            else if (valb(3)-tolb .le. 0.d0) then
                bmr(1) = valb(1)
                bmr(2) = valb(2)
                do i = 1, 6
                    epsr(i) = epi(i)
                end do
                tot3 = .true.
            end if
!
            call lceob2(intmax, toler, epsr, bmr, dm, &
                        lambda, mu, alpha, ecrob, ecrod, &
                        seuil, bdim, br, d, mult, &
                        elas, dbloq, iret)
!--VERIF SUR ENDO FINAL POUR VOIR SI ENDO DEPASSE 1 OU PAS
!-- SI ENDO DEPASSE 1 PLUS QUE TOLERANCE ON PASSE DANS LCEOBB
!-- QUI DECOUPE L INCREMENT DE CHARGE POUR ALLER DOUCEMENT A ENDO=1
!-- ENSUITE ON REVIENT AU 3D DANS REPERE INITIAL
!
            call diago3(br, vecbr, valbr)
            do i = 1, 2
                if (valbr(i) .lt. 0) then
                    reinit = .true.
                end if
                if (valbr(i)-tolb .le. 0.d0) then
                    valbr(i) = tolb-r8prem()
                    total = .true.
                end if
            end do
            if (d .gt. 1.d0) then
                reinit = .true.
            end if
            if (un-d-tolb .le. 0.d0) then
                d = un-tolb+r8prem()
                dbloq = .true.
            end if
            if (.not. reinit) then
                if (total) then
                    call r8inir(6, 0.d0, br, 1)
                    do i = 1, 3
                        do j = i, 3
                            do k = 1, 3
                                br(t(i, j)) = br(t(i, j))+vecbr(i, k)* &
                                              valbr(k)*vecbr(j, k)
                            end do
                        end do
                    end do
                end if
                call r8inir(6, 0.d0, binter, 1)
                if (tot1) then
                    binter(1) = tolb-r8prem()
                    binter(2) = br(1)
                    binter(3) = br(2)
                    binter(6) = br(4)
                else if (tot2) then
                    binter(1) = br(2)
                    binter(2) = tolb-r8prem()
                    binter(3) = br(1)
                    binter(5) = br(4)
                else if (tot3) then
                    binter(1) = br(1)
                    binter(2) = br(2)
                    binter(3) = tolb-r8prem()
                    binter(4) = br(4)
                end if
                call r8inir(9, 0.d0, interm, 1)
                call r8inir(6, 0.d0, b, 1)
                do i = 1, 3
                    do l = 1, 3
                        do k = 1, 3
                            interm(i, l) = interm(i, l)+vecb(i, k)*binter( &
                                           t(k, l))
                        end do
                        do j = i, 3
                            b(t(i, j)) = b(t(i, j))+interm(i, l)*vecb(j, l)
                        end do
                    end do
                end do
            else
                call lceobb(intmax, toler, epsm, deps, bm, &
                            dm, lambda, mu, alpha, ecrob, &
                            ecrod, rk, rk1, rk2, b, &
                            d, mult, elas, dbloq, iret)
            end if
!
!----------------------------------------------------------------
!----CAS OU 2 VALEURS PROPRES SONT NULLES------------------------
!----------------------------------------------------------------
!
        else if (bdim .eq. 1) then
!
!-- ON RESTREINT L ESPACE CAR L ENDO N EVOLUE PLUS DANS UNE DIRECTION
!
            call r8inir(9, 0.d0, interm, 1)
            call r8inir(6, 0.d0, epi, 1)
            do i = 1, 3
                do l = 1, 3
                    do k = 1, 3
                        interm(i, l) = interm(i, l)+vecb(k, i)*eps(t(k, l))
                    end do
                    do j = i, 3
                        epi(t(i, j)) = epi(t(i, j))+interm(i, l)*vecb(l, j)
                    end do
                end do
            end do
            tot1 = .false.
            tot2 = .false.
            tot3 = .false.
            call r8inir(6, 0.d0, bmr, 1)
            if (valb(1)-tolb .gt. 0.d0) then
                bmr(1) = valb(1)
                do i = 1, 6
                    epsr(i) = epi(i)
                end do
                tot1 = .true.
            else if (valb(2)-tolb .gt. 0.d0) then
                bmr(1) = valb(2)
                do i = 1, 6
                    epsr(i) = epi(r1(i))
                end do
                tot2 = .true.
            else if (valb(3)-tolb .gt. 0.d0) then
                bmr(1) = valb(3)
                do i = 1, 6
                    epsr(i) = epi(r2(i))
                end do
                tot3 = .true.
            end if
!
            call lceob1(intmax, toler, epsr, bmr, dm, &
                        lambda, mu, alpha, ecrob, ecrod, &
                        seuil, bdim, br, d, mult, &
                        elas, dbloq, iret)
!--VERIF SUR ENDO FINAL POUR VOIR SI ENDO DEPASSE 1 OU PAS
!-- SI ENDO DEPASSE 1 PLUS QUE TOLERANCE ON PASSE DANS LCEOBB
!-- QUI DECOUPE L INCREMENT DE CHARGE POUR ALLER DOUCEMENT A ENDO=1
!-- ENSUITE ON REVIENT AU 3D DANS REPERE INITIAL
!
            if (br(1) .lt. 0) then
                reinit = .true.
            end if
            if (br(1)-tolb .le. 0.d0) then
                br(1) = tolb-r8prem()
            end if
            if (d .gt. 1.d0) then
                reinit = .true.
            end if
            if (un-d-tolb .le. 0.d0) then
                d = un-tolb+r8prem()
                dbloq = .true.
            end if
            if (.not. reinit) then
                valb(1) = tolb-r8prem()
                valb(2) = tolb-r8prem()
                valb(3) = tolb-r8prem()
                if (tot1) valb(1) = br(1)
                if (tot2) valb(2) = br(1)
                if (tot3) valb(3) = br(1)
                call r8inir(6, 0.d0, b, 1)
                do i = 1, 3
                    do j = i, 3
                        do k = 1, 3
                            b(t(i, j)) = b(t(i, j))+vecb(i, k)*valb(k)* &
                                         vecb(j, k)
                        end do
                    end do
                end do
            else
                call lceobb(intmax, toler, epsm, deps, bm, &
                            dm, lambda, mu, alpha, ecrob, &
                            ecrod, rk, rk1, rk2, b, &
                            d, mult, elas, dbloq, iret)
            end if
!
        end if
!
!-----------------------------------------------------
!--ON RECUPERE L ENDOMMAGEMENT ET LA CONTRAINTE
!-----------------------------------------------------
        do i = 1, 3
            vip(i) = un-b(i)
        end do
        do i = 4, 6
            vip(i) = -b(i)
        end do
        vip(7) = d
        call sigeob(eps, b, d, 3, lambda, &
                    mu, sigm)
    end if
!
! ================================================================
!                            MATRICE TANGENTE
! ================================================================
!
!-on verifie l etat de B et D finaux pour calcul matrice tangente
!
    if (rigi) then
!
        call r8inir(36, 0.d0, dsidep, 1)
        if ((b(1) .eq. 1.d0) .and. (b(2) .eq. 1.d0) .and. &
            (b(3) .eq. 1.d0) .and. (d .eq. 0.d0)) then
            do i = 1, 6
                dsidep(i, i) = dsidep(i, i)+deux*mu
            end do
            do i = 1, 3
                do j = 1, 3
                    dsidep(i, j) = dsidep(i, j)+lambda
                end do
            end do
!
        else
!
            ad = (un-d)**deux
!
!
! -- DSIGMA/DEPS A ENDOMMAGEMENT CONSTANT = MATRICE SECANTE----
!
            call r8inir(6, 0.d0, cc, 1)
            do i = 1, 3
                do j = i, 3
                    do k = 1, 3
                        cc(t(i, j)) = cc(t(i, j))+b(t(i, k))*eps(t(k, j))+ &
                                      b(t(j, k))*eps(t(k, i))
                    end do
                end do
            end do
            treps = 0.d0
            treb = 0.d0
            do i = 1, 3
                treb = treb+cc(i)/deux
                treps = treps+eps(i)
            end do
            if (treb .ge. 0.d0) then
                do i = 1, 6
                    if (i .gt. 3) then
                        rtemp2 = rac2
                    else
                        rtemp2 = 1.d0
                    end if
                    do j = 1, 6
                        if (j .gt. 3) then
                            rtemp3 = rac2
                        else
                            rtemp3 = 1.d0
                        end if
                        dsidep(i, j) = dsidep(i, j)+lambda*b(i)*b(j)* &
                                       rtemp2*rtemp3
                    end do
                end do
            end if
            if (treps .lt. 0.d0) then
                do i = 1, 6
                    if (i .gt. 3) then
                        rtemp2 = rac2
                    else
                        rtemp2 = 1.d0
                    end if
                    do j = 1, 6
                        if (j .gt. 3) then
                            rtemp3 = rac2
                        else
                            rtemp3 = 1.d0
                        end if
                        dsidep(i, j) = dsidep(i, j)+lambda*ad*kron(i)* &
                                       kron(j)*rtemp2*rtemp3
                    end do
                end do
            end if
            call dfmdf(6, eps, zozo)
            do i = 1, 6
                do j = 1, 6
                    dsidep(i, j) = dsidep(i, j)+deux*ad*mu*zozo(i, j)
                end do
            end do
            call dfpdf(6, cc, zaza)
            call r8inir(36, 0.d0, bobo, 1)
            do i = 1, 3
                do j = i, 3
                    if (i .eq. j) then
                        rtemp2 = 1.d0
                    else
                        rtemp2 = rac2
                    end if
                    do p = 1, 3
                        do q = 1, 3
                            if (p .eq. q) then
                                rtemp3 = 1.d0
                            else
                                rtemp3 = un/rac2
                            end if
                            do k = 1, 3
                                if (k .eq. i) then
                                    rtemp4 = 1.d0
                                else
                                    rtemp4 = un/rac2
                                end if
                                if (k .eq. j) then
                                    rtemp5 = 1.d0
                                else
                                    rtemp5 = un/rac2
                                end if
                                do m = 1, 3
                                    do n = 1, 3
                                        if (m .eq. n) then
                                            rtemp6 = 1.d0
                                        else
                                            rtemp6 = 1/rac2
                                        end if
                                        bobo(t(i, j), t(p, q)) = bobo(t(i, j), t(p, q))+ &
                                                                 (zaza(t(i, k), t(m, n))* &
                                                                  (b(t(m, p))*kron(t(q, n))+ &
                                                                   kron(t(m, p))*b(t(q, n)))* &
                                                                  rtemp4*b(t(k, j))+ &
                                                                  (zaza(t(k, j), t(m, n))* &
                                                                   rtemp5*(b(t(m, p))* &
                                                                           kron(t(q, n))+ &
                                                                           kron(t(m, p))* &
                                                                           b(t(q, n)))* &
                                                                   b(t(i, k))))* &
                                                                 rtemp2*rtemp3*rtemp6
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            do i = 1, 6
                do j = 1, 6
                    dsidep(i, j) = dsidep(i, j)+mu/deux*bobo(i, j)
                end do
            end do
            call r8inir(36, 0.d0, dsisup, 1)
            if (option(10:14) .ne. '_ELAS') then
                if (.not. elas) then
!
                    if (bdim .eq. 3) then
!
                        call r8inir(6, 0.d0, deltab, 1)
                        do i = 1, 6
                            deltab(i) = b(i)-bm(i)
                        end do
                        deltad = d-dm
!
                        call meobl3(eps, b, d, deltab, deltad, &
                                    mult, lambda, mu, ecrob, ecrod, &
                                    alpha, rk1, rk2, bdim, dsisup)
!
                    else if (bdim .eq. 2) then
!
                        call diago3(b, vecb, valb)
!
                        call r8inir(9, 0.d0, interm, 1)
                        call r8inir(9, 0.d0, interb, 1)
                        call r8inir(6, 0.d0, epi, 1)
                        call r8inir(6, 0.d0, epib, 1)
                        do i = 1, 3
                            do l = 1, 3
                                do k = 1, 3
                                    interm(i, l) = interm(i, l)+vecb(k, i)* &
                                                   eps(t(k, l))
                                    interb(i, l) = interb(i, l)+vecb(k, i)* &
                                                   bm(t(k, l))
                                end do
                                do j = i, 3
                                    epi(t(i, j)) = epi(t(i, j))+interm(i, &
                                                                       l)*vecb(l, j)
                                    epib(t(i, j)) = epib(t(i, j))+interb( &
                                                    i, l)*vecb(l, j)
                                end do
                            end do
                        end do
                        tot1 = .false.
                        tot2 = .false.
                        tot3 = .false.
                        call r8inir(6, 0.d0, br, 1)
                        call r8inir(6, 0.d0, deltab, 1)
                        if (abs(valb(1))-tolb .le. 0.d0) then
!
                            br(1) = valb(2)
                            br(2) = valb(3)
                            do i = 1, 6
                                epsr(i) = epi(r1(i))
                            end do
                            deltab(1) = valb(2)-epib(2)
                            deltab(2) = valb(3)-epib(3)
                            deltab(4) = -epib(6)
                            deltad = d-dm
                            tot1 = .true.
                        else if (abs(valb(2))-tolb .le. 0.d0) then
!
                            br(1) = valb(3)
                            br(2) = valb(1)
                            do i = 1, 6
                                epsr(i) = epi(r2(i))
                            end do
                            deltab(1) = valb(3)-epib(3)
                            deltab(2) = valb(1)-epib(1)
                            deltab(4) = -epib(5)
                            deltad = d-dm
                            tot2 = .true.
!
                        else if (abs(valb(3))-tolb .le. 0.d0) then
!
                            br(1) = valb(1)
                            br(2) = valb(2)
                            do i = 1, 6
                                epsr(i) = epi(i)
                            end do
                            deltab(1) = valb(1)-epib(1)
                            deltab(2) = valb(2)-epib(2)
                            deltab(4) = -epib(4)
                            deltad = d-dm
                            tot3 = .true.
!
                        end if
!
                        call meobl2(epsr, br, d, deltab, deltad, &
                                    mult, lambda, mu, ecrob, ecrod, &
                                    alpha, rk1, rk2, bdim, dsiint)
!
                        if (tot1) then
                            do i = 1, 6
                                do j = 1, 6
                                    if (i .ge. 4) then
                                        rtemp1 = un/rac2
                                    else
                                        rtemp1 = 1.d0
                                    end if
                                    if (j .ge. 4) then
                                        rtemp2 = un/rac2
                                    else
                                        rtemp2 = 1.d0
                                    end if
                                    dsimed(r1(i), r1(j)) = dsiint(i, j)* &
                                                           rtemp1*rtemp2
                                end do
                            end do
                        else if (tot2) then
                            do i = 1, 6
                                do j = 1, 6
                                    if (i .ge. 4) then
                                        rtemp1 = un/rac2
                                    else
                                        rtemp1 = 1.d0
                                    end if
                                    if (j .ge. 4) then
                                        rtemp2 = un/rac2
                                    else
                                        rtemp2 = 1.d0
                                    end if
                                    dsimed(r2(i), r2(j)) = dsiint(i, j)* &
                                                           rtemp1*rtemp2
                                end do
                            end do
                        else if (tot3) then
                            do i = 1, 6
                                do j = 1, 6
                                    if (i .ge. 4) then
                                        rtemp1 = un/rac2
                                    else
                                        rtemp1 = 1.d0
                                    end if
                                    if (j .ge. 4) then
                                        rtemp2 = un/rac2
                                    else
                                        rtemp2 = 1.d0
                                    end if
                                    dsimed(i, j) = dsiint(i, j)*rtemp1* &
                                                   rtemp2
                                end do
                            end do
                        end if
!
                        call r8inir(36, 0.d0, dsisup, 1)
!
                        do i = 1, 3
                            do j = i, 3
                                if (i .eq. j) then
                                    rtemp1 = 1.d0
                                else
                                    rtemp1 = rac2
                                end if
                                do p = 1, 3
                                    do q = p, 3
                                        if (p .eq. q) then
                                            rtemp2 = 1.d0
                                        else
                                            rtemp2 = rac2
                                        end if
                                        do k = 1, 3
                                            do l = 1, 3
                                                do m = 1, 3
                                                    do n = 1, 3
                                                        dsisup(t(i, j), t(p, q)) = &
                                                            dsisup(t(i, j), t(p, q))+ &
                                                            vecb(i, k)* &
                                                            vecb(j, l)*vecb(p, m)* &
                                                            vecb(q, n) &
                                                            *dsimed(t(k, l), t(m, n))* &
                                                            rtemp1*rtemp2
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
!
!--------------------------------------------------------
!--------------------------------------------------------
!
                    else if (bdim .eq. 1) then
!
                        call diago3(b, vecb, valb)
!
                        call r8inir(9, 0.d0, interm, 1)
                        call r8inir(9, 0.d0, interb, 1)
                        call r8inir(6, 0.d0, epi, 1)
                        call r8inir(6, 0.d0, epib, 1)
                        do i = 1, 3
                            do l = 1, 3
                                do k = 1, 3
                                    interm(i, l) = interm(i, l)+vecb(k, i)* &
                                                   eps(t(k, l))
                                    interb(i, l) = interb(i, l)+vecb(k, i)* &
                                                   bm(t(k, l))
                                end do
                                do j = i, 3
                                    epi(t(i, j)) = epi(t(i, j))+interm(i, &
                                                                       l)*vecb(l, j)
                                    epib(t(i, j)) = epib(t(i, j))+interb( &
                                                    i, l)*vecb(l, j)
                                end do
                            end do
                        end do
!
                        tot1 = .false.
                        tot2 = .false.
                        tot3 = .false.
                        call r8inir(6, 0.d0, br, 1)
                        call r8inir(6, 0.d0, deltab, 1)
!
                        if (abs(valb(1))-tolb .gt. 0.d0) then
!
                            br(1) = valb(1)
                            do i = 1, 6
                                epsr(i) = epi(i)
                            end do
                            deltab(1) = valb(1)-epib(1)
                            deltad = d-dm
                            tot1 = .true.
!
                        else if (abs(valb(2))-tolb .gt. 0.d0) then
!
                            br(1) = valb(2)
!
                            do i = 1, 6
                                epsr(i) = epi(r1(i))
                            end do
                            deltab(1) = valb(2)-epib(2)
                            deltad = d-dm
                            tot2 = .true.
!
                        else if (abs(valb(3))-tolb .gt. 0.d0) then
!
                            br(1) = valb(3)
!
                            do i = 1, 6
                                epsr(i) = epi(r2(i))
                            end do
                            deltab(1) = valb(3)-epib(3)
                            deltad = d-dm
                            tot3 = .true.
!
                        end if
!
                        call meobl1(epsr, b, d, deltab, deltad, &
                                    mult, lambda, mu, ecrob, ecrod, &
                                    alpha, rk1, rk2, bdim, dsiint)
!
                        if (tot1) then
                            do i = 1, 6
                                do j = 1, 6
                                    if (i .ge. 4) then
                                        rtemp1 = un/rac2
                                    else
                                        rtemp1 = 1.d0
                                    end if
                                    if (j .ge. 4) then
                                        rtemp2 = un/rac2
                                    else
                                        rtemp2 = 1.d0
                                    end if
                                    dsimed(i, j) = dsiint(i, j)*rtemp1* &
                                                   rtemp2
                                end do
                            end do
                        else if (tot2) then
                            do i = 1, 6
                                do j = 1, 6
                                    if (i .ge. 4) then
                                        rtemp1 = un/rac2
                                    else
                                        rtemp1 = 1.d0
                                    end if
                                    if (j .ge. 4) then
                                        rtemp2 = un/rac2
                                    else
                                        rtemp2 = 1.d0
                                    end if
                                    dsimed(r1(i), r1(j)) = dsiint(i, j)* &
                                                           rtemp1*rtemp2
                                end do
                            end do
                        else if (tot3) then
                            do i = 1, 6
                                do j = 1, 6
                                    if (i .ge. 4) then
                                        rtemp1 = un/rac2
                                    else
                                        rtemp1 = 1.d0
                                    end if
                                    if (j .ge. 4) then
                                        rtemp2 = un/rac2
                                    else
                                        rtemp2 = 1.d0
                                    end if
                                    dsimed(r2(i), r2(j)) = dsiint(i, j)* &
                                                           rtemp1*rtemp2
                                end do
                            end do
                        end if
!
                        call r8inir(36, 0.d0, dsisup, 1)
!
                        do i = 1, 3
                            do j = i, 3
                                if (i .eq. j) then
                                    rtemp1 = 1.d0
                                else
                                    rtemp1 = rac2
                                end if
                                do p = 1, 3
                                    do q = p, 3
                                        if (p .eq. q) then
                                            rtemp2 = 1.d0
                                        else
                                            rtemp2 = rac2
                                        end if
                                        do k = 1, 3
                                            do l = 1, 3
                                                do m = 1, 3
                                                    do n = 1, 3
                                                        dsisup(t(i, j), t(p, q)) = &
                                                            dsisup(t(i, j), t(p, q))+ &
                                                            vecb(i, k)* &
                                                            vecb(j, l)*vecb(p, m)*vecb(q, n) &
                                                            *dsimed(t(k, l), t(m, n))*rtemp1* &
                                                            rtemp2
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
!
                    end if
!
                end if
            end if
!
            do i = 1, 6
                do j = 1, 6
                    dsidep(i, j) = dsidep(i, j)+dsisup(i, j)
                end do
            end do
!
        end if
    end if
end subroutine
