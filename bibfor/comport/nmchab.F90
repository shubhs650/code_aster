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
subroutine nmchab(fami, kpg, ksp, ndim, &
                  imate, relaComp, carcri, instam, instap, &
                  deps, sigm, vim, option, sigp, &
                  vip, dsidep, iret)
!
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmcham.h"
#include "asterfort/nmchat.h"
#include "asterfort/nmchdp.h"
#include "asterfort/r8inir.h"
#include "asterfort/radial.h"
#include "asterfort/trace.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
!
    character(len=*) :: fami
    integer(kind=8) :: kpg
    integer(kind=8) :: ksp
    integer(kind=8) :: ndim
    integer(kind=8) :: imate
    character(len=16), intent(in) :: relaComp, option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: instam
    real(kind=8) :: instap
    real(kind=8) :: deps(6)
    real(kind=8) :: sigm(6)
    real(kind=8) :: vim(*)
    real(kind=8) :: sigp(6)
    real(kind=8) :: vip(*)
    real(kind=8) :: dsidep(6, 6)
    integer(kind=8) :: iret
!     INTEGRATION LOCALE DES LOIS DE COMPORTEMENT DE CHABOCHE
!          RELATIONS : 'VMIS_CIN1_CHAB' 'VMIS_CIN2_CHAB'
!          RELATIONS : 'VISC_CIN1_CHAB' 'VISC_CIN2_CHAB'
!          RELATIONS : 'VMIS_CIN2_MEMO' 'VISC_CIN2_MEMO'
!
!     ARGUMENTS :
!       IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!       IN      KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!       IN      NDIM    DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!       IN      TYPMOD  TYPE DE MODELISATION
!               IMATE   ADRESSE DU MATERIAU CODE
!               COMPOR  COMPORTEMENT DE L ELEMENT
!                       COMPOR(1) = RELATION DE COMPORTEMENT
!                       COMPOR(2) = NB DE VARIABLES INTERNES
!                       COMPOR(3) = TYPE DE DEFORMATION
!               CRIT    CRITERES  LOCAUX, EN PARTICULIER
!                       CRIT(1) = NOMBRE D ITERATIONS MAXI A CONVERGENCE
!                                 (ITER_INTE_MAXI == ITECREL)
!                       CRIT(3) = VALEUR DE LA TOLERANCE DE CONVERGENCE
!                                 (RESI_INTE == RESCREL)
!               INSTAM  INSTANT T
!               INSTAP  INSTANT T+DT
!               DEPS    INCREMENT DE DEFORMATION TOTALE
!               SIGM    CONTRAINTE A T
!               VIM     VARIABLES INTERNES A T    + INDICATEUR ETAT T
!    ATTENTION  VIM     VARIABLES INTERNES A T MODIFIEES SI REDECOUPAGE
!               OPTION     OPTION DE CALCUL A FAIRE
!                             'RIGI_MECA_TANG'> DSIDEP(T)
!                             'FULL_MECA'     > DSIDEP(T+DT) , SIG(T+DT)
!                             'RAPH_MECA'     > SIG(T+DT)
!       OUT     SIGP    CONTRAINTE A T+DT
!               VIP     VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSIDEP  MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!               IRET    CODE RETOUR DE  L'INTEGRATION DE LA LDC
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ABSENCE DE CONVERGENCE
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
!               -----------------------------------------------------
!

    integer(kind=8) ::  nbvar
    real(kind=8) :: depsth(6), pm, c2inf, gamm20
    real(kind=8) :: plast, depsmo, sigmmo, e, nu, troisk, deuxmu
    real(kind=8) :: rpm, sieleq, seuil, dp, cm, ainf, cp, rpvm, rpvp
    real(kind=8) :: coef, sigedv(6), kron(6), depsdv(6), epspm(6)
    real(kind=8) :: sigmdv(6), sigpdv(6), em, num, ksip(6), qp
    real(kind=8) :: troikm, deumum, sigmp(6), sigel(6), pp, epspp(6)
    real(kind=8) :: un, rac2, c2p, radi, gamma0, gammap, delta1, delta2
    real(kind=8) :: r0, rinf, b, cinf, k, w, mat(18), c2m, gamm2p, unrac2
    real(kind=8) :: depsp(6), alfam(6), alfa(6), dalfa(6), n1, n2
    real(kind=8) :: alfa2m(6), alfa2(6), dalfa2(6), matel(4), xm(6), xp(6)
    real(kind=8) :: dt, ksim(6), qm, beta1, beta2
    common/fchab/mat, pm, sigedv, epspm, alfam, alfa2m, deuxmu, rpvm, rpvp,&
     &    qm, qp, ksim, ksip, dt, n1, n2, depsp,&
     &    beta1, beta2, ndimsi, nbvar, visc, memo, idelta
    integer(kind=8) :: ndimsi, i, niter, visc, memo, idelta
    blas_int :: b_incx, b_incy, b_n
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
    iret = 0
    call nmcham(fami, kpg, ksp, imate, relaComp, &
                matel, mat, nbvar, memo, visc, &
                idelta, coef)
!
!     NBVARI=2+6*NBVAR+MEMO*14
    un = 1.0d0
    rac2 = sqrt(2.d0)
    unrac2 = 1.d0/sqrt(2.d0)
    em = matel(1)
    num = matel(2)
    deumum = em/(un+num)
    troikm = em/(un-2.d0*num)
    e = matel(3)
    nu = matel(4)
    deuxmu = e/(un+nu)
    troisk = e/(un-2.d0*nu)
    r0 = mat(1)
    rinf = mat(2)
    b = mat(3)
    cinf = mat(4)
    k = mat(5)
    w = mat(6)
    gamma0 = mat(7)
    ainf = mat(8)
    c2inf = mat(9)
    gamm20 = mat(10)
    if (idelta .gt. 0.d0) then
        delta1 = mat(17)
        delta2 = mat(18)
    else
        delta1 = 1.d0
        delta2 = 1.d0
    end if
    n1 = 1.d0
    n2 = 1.d0
!
! --- INITIALISATIONS :
!
    ndimsi = 2*ndim
    dp = 0.d0
    call r8inir(6, 0.d0, depsp, 1)
    seuil = 0.d0
    dt = instap-instam
    pm = vim(1)
    plast = vim(2)
    niter = 0
    if (memo .eq. 0) then
        rpm = rinf+(r0-rinf)*exp(-b*pm)
    else if (memo .eq. 1) then
        rpvm = vim(15)
        rpm = rpvm+r0
        qm = vim(16)
        call r8inir(6, 0.d0, ksim, 1)
        call r8inir(6, 0.d0, ksip, 1)
        call r8inir(6, 0.d0, epspm, 1)
        call r8inir(6, 0.d0, epspp, 1)
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vim(17), b_incx, ksim, b_incy)
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vim(23), b_incx, epspm, b_incy)
        qp = qm
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, ksim, b_incx, ksip, b_incy)
        rpvp = rpvm
    end if
    cm = cinf*(un+(k-un)*exp(-w*pm))
    c2m = c2inf*(un+(k-un)*exp(-w*pm))
!
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vim(3), b_incx, alfam, b_incy)
    b_n = to_blas_int(ndimsi-3)
    b_incx = to_blas_int(1)
    call dscal(b_n, rac2, alfam(4), b_incx)
    if (nbvar .eq. 2) then
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vim(9), b_incx, alfa2m, b_incy)
        b_n = to_blas_int(ndimsi-3)
        b_incx = to_blas_int(1)
        call dscal(b_n, rac2, alfa2m(4), b_incx)
    else
        call r8inir(6, 0.d0, alfa2m, 1)
    end if
!
! --- CALCUL DE DEPSMO ET DEPSDV :
    depsmo = 0.d0
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, deps, b_incx, depsth, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -coef, kron, b_incx, depsth, &
               b_incy)
    depsmo = trace(3, depsth)/3.d0
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, depsth, b_incx, depsdv, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -depsmo, kron, b_incx, depsdv, &
               b_incy)
!
!       -------------------------------------------------
    sigmmo = trace(3, sigm)/3.d0
    do i = 1, ndimsi
        sigmp(i) = deuxmu/deumum*(sigm(i)-sigmmo*kron(i))+troisk/troikm*sigmmo*kron(i)
    end do
!     SIGMMO A PU CHANGER A CAUSE DE TROISK/TROIKM
    sigmmo = trace(3, sigmp)/3.d0
!
! --- CALCUL DU SEUIL :
    sieleq = 0.d0
    do i = 1, ndimsi
        sigmdv(i) = sigmp(i)-sigmmo*kron(i)
        sigedv(i) = sigmdv(i)+deuxmu*depsdv(i)
        sigel(i) = sigedv(i)-cm*alfam(i)/1.5d0-c2m*alfa2m(i)/1.5d0
        sieleq = sieleq+sigel(i)*sigel(i)
    end do
    sieleq = sqrt(1.5d0*sieleq)
    seuil = sieleq-rpm
!
! --- CALCUL DE SIGP,VIP
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
!
!       CALCUL DE DP :
!       ------------
        if (seuil .le. 0.d0) then
            plast = 0.d0
            dp = 0.d0
        else
! ---       DETERMINATION DE DP SOLUTION D'UNE EQUATION NON LINEAIRE
            call nmchdp(carcri, seuil, dp, iret, niter)
            if (iret .gt. 0) goto 999
            plast = un
        end if
!
        pp = pm+dp
        gammap = gamma0*(ainf+(un-ainf)*exp(-b*pp))
        gamm2p = gamm20*(ainf+(un-ainf)*exp(-b*pp))
!
        if (memo .eq. 1) then
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, epspm, b_incx, epspp, b_incy)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, depsp, b_incx, epspp, &
                       b_incy)
        end if
!
        do i = 1, ndimsi
            if (cinf .ne. 0.d0) then
                dalfa(i) = (n1*depsp(i)-gammap*delta1*dp*alfam(i))
                dalfa(i) = dalfa(i)/(un+gammap*delta1*dp)
            else
                dalfa(i) = 0.d0
            end if
            if (nbvar .eq. 2) then
                if (c2inf .ne. 0.d0) then
                    dalfa2(i) = (n2*depsp(i)-gamm2p*delta2*dp*alfa2m(i))
                    dalfa2(i) = dalfa2(i)/(un+gamm2p*delta2*dp)
                else
                    dalfa2(i) = 0.d0
                end if
            end if
        end do
!
        do i = 1, ndimsi
            alfa(i) = alfam(i)+dalfa(i)
            if (nbvar .eq. 2) then
                alfa2(i) = alfa2m(i)+dalfa2(i)
            end if
            sigpdv(i) = sigedv(i)-deuxmu*depsp(i)
            sigp(i) = sigpdv(i)+(sigmmo+troisk*depsmo)*kron(i)
        end do
!
    end if
!
! ---- CALCUL DE LA MATRICE DE COMPORTEMENT TANGENTE COHERENTE DSIDEP
    if (option(1:14) .eq. 'RIGI_MECA_TANG' .or. option(1:9) .eq. 'FULL_MECA') then
        call nmchat(matel, mat, nbvar, memo, visc, &
                    plast, sigmdv, depsdv, pm, dp, &
                    ndimsi, dt, rpvp, qp, vim, &
                    idelta, n1, n2, beta1, beta2, &
                    dsidep)
!
    end if
!
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        vip(1) = pp
        vip(2) = niter
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, alfa, b_incx, vip(3), b_incy)
        b_n = to_blas_int(ndimsi-3)
        b_incx = to_blas_int(1)
        call dscal(b_n, unrac2, vip(6), b_incx)
        if (nbvar .eq. 2) then
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, alfa2, b_incx, vip(9), b_incy)
            b_n = to_blas_int(ndimsi-3)
            b_incx = to_blas_int(1)
            call dscal(b_n, unrac2, vip(12), b_incx)
        end if
        if (memo .eq. 1) then
            vip(15) = rpvp
            vip(16) = qp
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, ksip, b_incx, vip(17), b_incy)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, epspp, b_incx, vip(23), b_incy)
        end if
    end if
!
!     Critere de radialite
    if (option(1:9) .ne. 'RIGI_MECA') then
        if (carcri(10) .gt. 0.d0) then
!           CALCUL DE X1, X2
            cp = cinf*(un+(k-un)*exp(-w*pp))
            c2p = c2inf*(un+(k-un)*exp(-w*pp))
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, alfam, b_incx, xm, b_incy)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, alfa, b_incx, xp, b_incy)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            call dscal(b_n, cm/1.5d0, xm, b_incx)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            call dscal(b_n, cp/1.5d0, xp, b_incx)
            if (nbvar .eq. 2) then
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, c2m/1.5d0, alfa2m, b_incx, xm, &
                           b_incy)
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, c2p/1.5d0, alfa2, b_incx, xp, &
                           b_incy)
            end if
            call radial(ndimsi, sigm, sigp, vim(2), vip(2), &
                        1, xm, xp, radi)
!
            if (radi .gt. carcri(10)) then
                iret = 2
            end if
        end if
    end if
!
999 continue
end subroutine
