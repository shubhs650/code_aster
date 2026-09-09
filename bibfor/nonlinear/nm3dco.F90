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

subroutine nm3dco(fami, kpg, ksp, ndim, option, &
                  imate, sigm, deps, vim, sigp, &
                  vip, dsidep, carcri, codret)
!
    implicit none
! ----------------------------------------------------------------------
!          LOI DE L'ACIER SOUMIS A LA CORROSION 3D
!
! IN  NDIM    : DIMENSION
! IN  OPTION  : OPTION DE CALCUL
! IN  IMATE   : POINTEUR MATERIAU
! IN  SIGM    : CONTRAINTE AU TEMPS MOINS
! IN  EPSM    : DEFORMATION TOTALE AU TEMPS MOINS
! IN  DEPS    : DEFORMATION  TOTALE PLUS - DEFORMATION TOTALE MOINS
! IN VIM      : VARIABLE INTERNES AU TEMPS MOINS
! IN CRILDC   : 1 ITERMAX, 3 RESI_INTE
!
! OUT SIGP     : CONTRAINTES PLUS
! OUT VIP       : VARIABLE INTERNES PLUS
! OUT DSIDEP    : DSIG/DEPS
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"

    real(kind=8) :: sigm(6), deps(6), vim(*)
    real(kind=8) :: sigp(6), vip(*), dsidep(6, 6)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16) :: option
    character(len=*) :: fami
    integer(kind=8) :: ndim, imate, codret, kpg, ksp
!
    real(kind=8) :: young, nu, kcoef, mcoef, coefdc, limit
    real(kind=8) :: valres(4)
    integer(kind=8) :: codres(4)
    character(len=16) :: nomres(4)
!
    integer(kind=8) :: iter, itemax, ndimsi, i, k, l, m, itd, ibid
!
    real(kind=8) :: resi, ecum, ecumm, dcoef, plas
    real(kind=8) :: defe, defc, nuetoi, defpc(3), ecumc, ecumd
    real(kind=8) :: coef1, coef2, j2, rini, crit, crit0
    aster_logical :: dconv, pconv, premd, mtang, melas
    real(kind=8) :: terme1, terme2(6), terme4(6), terme5, ter11
    real(kind=8) :: deltap, dp
    real(kind=8) :: criten, crit2, crit0d, treps
!
    real(kind=8) :: sigfi(6), sigd(6), rbid, drdp, hp
    real(kind=8) :: kci, lamda, deumu, acoef, bcoef, ccoef, corrm
!
    ndimsi = 2*ndim
    nomres(1) = 'D_CORR'
    nomres(2) = 'ECRO_K'
    nomres(3) = 'ECRO_M'
    nomres(4) = 'SY'
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'CORR_ACIER', 0, ' ', [0.d0], &
                4, nomres, valres, codres, 1)
!
    coefdc = valres(1)
    kcoef = valres(2)
    mcoef = valres(3)
    limit = valres(4)
    nomres(1) = 'NU'
    nomres(2) = 'E'
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres, valres, codres, 1)
!
    nu = valres(1)
    young = valres(2)
!
    ecumm = vim(1)
    dcoef = vim(2)
    plas = vim(3)
!
! --- PARAMETRES DE CONVERGENCE
    resi = carcri(3)
    itemax = nint(carcri(1))
!
!
!     CALCUL DE LA DEFORMATION CRITIQUE
    defe = limit/young
    call rcvarc('F', 'CORR', '-', fami, kpg, &
                ksp, corrm, ibid)
    if (corrm .le. 15.d0) then
        defc = 2.345d-01-(1.11d-02*corrm)
    else
        defc = 5.1d-02-(6.d-04*corrm)
    end if
    nuetoi = 0.5d0-(defe/defc)*(0.5d0-nu)
    defpc(1) = defc
    defpc(2) = -1.d0*nuetoi*defc
    defpc(3) = defpc(2)
!
!     CALCUL DE LA DEFORMATION PLAST EQUIV CRITIQUE
    ecumc = defpc(1)*defpc(1)
    ecumc = ecumc+defpc(2)*defpc(2)
    ecumc = ecumc+defpc(3)*defpc(3)
    ecumc = (2.d0/3.d0)*ecumc
    ecumc = ecumc**0.5d0
!
!     CALCUL DE DEFORMATION PLASTIQUE EQUIV DE DEBUT D'ENDOMMAGMENT
    ecumd = 0.8d0*ecumc
!
!     CALCUL DE TRIAXIALITE
!
!     PARAMETRES COEF2-LAMBDA ET COEF1-2MU
    coef1 = (young/(1.d0+nu))
    coef2 = (nu*young)/((1.d0+nu)*(1.d0-(2.d0*nu)))
!
!     DES INITIALISATIONS POUR MATRICE TANGENTE ?
    do i = 1, ndimsi
        sigp(i) = sigm(i)
        sigfi(i) = sigp(i)
    end do
!
!     DEFORMATION PLASTIQUE EQUIV A L'INSTANT M
    ecum = ecumm
!
!     CALCUL DES CONTRAINTES ELASTIQUES
    treps = deps(1)+deps(2)+deps(3)
    do i = 1, ndimsi
        sigp(i) = sigp(i)+coef1*deps(i)
    end do
    do i = 1, 3
        sigp(i) = sigp(i)+coef2*treps
    end do
!
    dp = 0.d0
    dconv = .false.
    iter = 0
    itd = 0
    premd = .true.
!
!
!   999=RETOUR ENDO
999 continue
!
    if (.not. dconv) then
!
!      CALCUL DE J2(SIG)
        j2 = 0.d0
        do i = 1, ndimsi
            j2 = j2+(sigp(i)**2)
        end do
        j2 = j2-((1.d0/3.d0)*((sigp(1)+sigp(2)+sigp(3))**2))
        j2 = (3.d0/2.d0*j2)**0.5d0
!
!      CALCUL D'ECROUISSAGE
        rini = kcoef*(ecum**(1.d0/mcoef))
!
!       SURFACE SEUIL
        crit0 = ((j2/(1.d0-dcoef))-rini-limit)
        crit = crit0
!
        if (option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RAPH_MECA') then
!
889         continue
!
            if (crit0 .lt. 0.d0) then
                plas = 0.d0
                vip(3) = plas
                dconv = .true.
                dp = 0.d0
! ON SORT COMPLETEMENT
                goto 999
            else
                plas = 1.d0
                vip(3) = plas
                pconv = .false.
!
!     PLASTICITE  (888 = RETOUR PLASTICITE)
888             continue
                if (.not. pconv) then
                    iter = iter+1
                    if (iter .eq. itemax) then
                        call utmess('A', 'MODELISA5_42')
                        codret = 1
                        do i = 1, ndimsi
                            sigp(i) = 0.d0
                        end do
                        vip(1) = vim(1)
                        vip(2) = vim(2)
                        vip(3) = vim(3)
                        goto 998
                    end if
!     TERME1 : F(SIG,R)
                    terme1 = (j2/(1.d0-dcoef))-rini-limit
!
!     TERME2(*) : DF(SIG,X,R) / DSIG
                    rbid = 1.d0/(j2*(1.d0-dcoef))
                    terme2(1) = rbid*(sigp(1)-0.5d0*sigp(2)-0.5d0*sigp(3))
                    terme2(2) = rbid*(sigp(2)-0.5d0*sigp(1)-0.5d0*sigp(3))
                    terme2(3) = rbid*(sigp(3)-0.5d0*sigp(1)-0.5d0*sigp(2))
                    do i = 4, ndimsi
                        terme2(i) = rbid*1.5d0*sigp(i)
                    end do
!
!     TERME3(*) : DF(SIG,X,R) / DSIG = DF(SIG,X,R) / DSIG
!
!     TERME4(*) : KE * TERME2
                    rbid = coef2*(terme2(1)+terme2(2)+terme2(3))
                    do i = 1, ndimsi
                        terme4(i) = coef1*terme2(i)
                    end do
                    do i = 1, 3
                        terme4(i) = terme4(i)+rbid
                    end do
!
!     TERME5 = TERME2 : TERME4
                    terme5 = 0.d0
                    do i = 1, ndimsi
                        terme5 = terme5+terme2(i)*terme4(i)
                    end do
!
!     TER11 : DF/DR*COEFFIC
                    ter11 = limit/kcoef
                    ter11 = j2/(kcoef*(1.d0-dcoef))-ter11
                    ter11 = ter11**(1.d0-mcoef)
                    ter11 = kcoef/mcoef*ter11
                    ter11 = -1.d0*ter11
!
!     DETERMINATION DE DELTAP
                    deltap = terme1/(terme5-ter11)
!
!      CALCUL  DE TOUTES LES VARIABLES INTERNES :
                    do i = 1, ndimsi
                        sigp(i) = sigp(i)-(deltap*terme4(i))
                    end do
!
!     DETERMINATION DE LA DEFORMATION PLASTIQUE ET P
                    dp = deltap/(1.d0-dcoef)
                    ecum = ecum+dp
!
!     CALCUL DE J2(SIG)
                    j2 = 0.d0
                    do i = 1, ndimsi
                        j2 = j2+(sigp(i)**2)
                    end do
                    j2 = j2-((1.d0/3.d0)*((sigp(1)+sigp(2)+sigp(3))**2))
                    j2 = (3.d0/2.d0*j2)**0.5d0
!
!     DETERMINATION DE L'ECROUISSAGE
                    rini = kcoef*(ecum**(1.d0/mcoef))
!
!     SURFACE SEUIL
                    crit = (j2/(1.d0-dcoef))-rini-limit
                    pconv = (abs(crit/crit0) .le. resi)
! FIN IF PCONV
                    goto 888
                end if
            end if
!
!     CRITERE D'ENDOMMAGEMENT
            criten = ecum-ecumd
            if (criten .le. 0.d0) then
                dconv = .true.
            else
!     COEFFICIENT D'ENDOMMAGEMENT
                dcoef = coefdc*(ecum-ecumd)/(ecumc-ecumd)
                if (dcoef .gt. 0.99d0) then
                    dconv = .true.
                    dcoef = 0.99d0
                    do i = 1, ndimsi
                        sigp(i) = 0.d0
                    end do
                end if
!
                crit2 = (j2/(1.d0-dcoef))-rini-limit
                if (premd) then
                    crit0d = crit2
                    premd = .false.
                end if
                dconv = (abs(crit2/crit0d) .le. resi)
!
! SI PAS CONVERGENCE EN ENDO, RETOUR A LA PLASTICITE
                if (.not. dconv) then
                    pconv = .false.
                    itd = itd+1
                    goto 889
                end if
            end if
!
! FIN IF OPTION
        end if
! FIN IF NOT DCONV
    end if
    vip(1) = ecum
    vip(2) = dcoef
!
!     CALCUL DE LA MATRICE TANGENTE OU ELAS OU SECANTE DECHARGE
    mtang = (option .eq. 'RIGI_MECA_TANG') .or.&
         &      (option .eq. 'FULL_MECA')
    melas = (option .eq. 'RIGI_MECA_ELAS') .or.&
         &      (option .eq. 'FULL_MECA_ELAS')
!
!     ELASTIQUE
    if ((option .eq. 'RIGI_MECA') .or. melas .or. (mtang .and. (plas .lt. 0.5d0))) then
!
        dsidep(:, :) = 0.d0
!
        do k = 1, 6
            dsidep(k, k) = coef1
        end do
        do k = 1, 3
            do l = 1, 3
                dsidep(k, l) = dsidep(k, l)+coef2
            end do
        end do
    end if
!
!     PLASTICITE
    if (mtang .and. (plas .ge. 0.5d0)) then
        kci = 1.d0
        if (option(1:14) .eq. 'RIGI_MECA_TANG') then
            rbid = sigfi(1)+sigfi(2)+sigfi(3)
            do k = 1, 3
                sigd(k) = sigfi(k)-rbid*(1.d0/3.d0)
            end do
            do k = 4, ndimsi
                sigd(k) = sigfi(k)
            end do
        else
            rbid = sigp(1)+sigp(2)+sigp(3)
            do k = 1, 3
                sigd(k) = sigp(k)-rbid*(1.d0/3.d0)
            end do
            do k = 4, ndimsi
                sigd(k) = sigp(k)
            end do
        end if
        drdp = (ecum**((1.d0/mcoef)-1.d0))
        drdp = (kcoef/mcoef)*drdp
!
        if (dcoef .le. 0.d0) then
!          PLASTICITE SANS ENDOMMAGEMENT
            hp = 1.d0+((3.d0/2.d0)*coef1*kci*dp)/((rini+limit))
            lamda = coef2+((coef1/3.d0)*(1.d0-(1.d0/hp)))
            deumu = coef1/hp
            bcoef = 1.d0-(((drdp*dp)/(rini+limit)))
            bcoef = kci*((9.d0*(coef1**2))/(4.d0*hp))*bcoef
            ccoef = drdp+((3.d0/2.d0)*coef1)
            dsidep(:, :) = 0.d0
            do k = 1, ndimsi
                dsidep(k, k) = deumu
            end do
            do k = 1, 3
                do m = 1, 3
                    dsidep(k, m) = dsidep(k, m)+lamda
                end do
            end do
            do k = 1, ndimsi
                do m = 1, ndimsi
                    dsidep(k, m) = dsidep(k, m)- &
                                   ((bcoef/ccoef)*((sigd(k)/(rini+limit))*(sigd(m)/(rini+limit))))
                end do
            end do
!
        else
!          PLASTICITE ET ENDOMMAGEMENT
            hp = (1.d0+((3.d0/2.d0)*coef1*kci*dp)/((1.d0-dcoef)*(rini+limit)))
            lamda = coef2+((coef1/3.d0)*(1.d0-(1.d0/hp)))
            deumu = coef1/hp
            acoef = (coefdc/(ecumc-ecumd))
            bcoef = ( &
                    1.d0-( &
                    dp*(((1.d0-dcoef)*drdp)-(rini*acoef))/((1.d0-dcoef)*(rini+limit))) &
                    )
            bcoef = kci*((9.d0*(coef1**2))/(4.d0*hp))*bcoef
            ccoef = (((1.d0-dcoef)*drdp)+((3.d0/2.d0)*coef1)-(rini*acoef))
!
            dsidep(:, :) = 0.d0
!
            do k = 1, ndimsi
                dsidep(k, k) = deumu
            end do
            do k = 1, 3
                do m = 1, 3
                    dsidep(k, m) = dsidep(k, m)+lamda
                end do
            end do
            do k = 1, ndimsi
                do m = 1, ndimsi
                    dsidep(k, m) = (dsidep(k, m)- &
                                    ((bcoef/ccoef)*((sigd(k)/((1.d0-dcoef)*(rini+limit))) &
                                                    *(sigd(m)/((1.d0-dcoef)*(rini+limit))))))
                end do
            end do
        end if
    end if
!
!     CAS RIGI_MECA_ELAS ET FULL_MECA_ELAS AVEC ENDOMMAGEMENT
    if (melas .and. (dcoef .ge. 0.d0)) then
        do k = 1, ndimsi
            do m = 1, ndimsi
                dsidep(k, m) = (1.d0-dcoef)*dsidep(k, m)
            end do
        end do
    end if
!
998 continue
end subroutine
