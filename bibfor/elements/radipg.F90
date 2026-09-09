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
subroutine radipg(sig1, sig2, npg, nbsig, radia, &
                  cosang, ind, relaComp, imate, nvi, &
                  vari1, vari2)
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/nmcham.h"
#include "asterfort/norsig.h"
#include "asterfort/radial.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer(kind=8) :: npg, nbsig, ind, nvi, imate
    real(kind=8) :: sig1(*), sig2(*), radia(*), cosang(*)
    real(kind=8) :: vari1(*), vari2(*)
    character(len=16) :: relaComp
!
!     BUT:
!       CALCUL DE L'INDICATEUR LOCAL DE PERTE DE RADIALITE RADIA
!       I = 1- ABS(SIG1:DSIGMA)/(NORME(SIG1)*NORME(DSIGMA)
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   SIG1     : CONTRAINTES INSTANT -
! IN   SIG2     : CONTRAINTES INSTANT +
! IN   NPG      : NOMBRE DE POINT DE GAUSS
! IN   NBSIG    : NOMBRE DE CMP DE CONTR
! IN   IND      : 0 : CALCUL DE RADI_V, 1 : CALCUL DE ERR_RADI
! IN   COMPOR   : COMPORTEMENT
! IN   IMATE    : ADRESSE MATERIAU CODE
! IN   NVI      : NOMBRE DE VARIABLES INTERNES
! IN   VARI1    : VARIABLES INTERNES INSTANT -
! IN   VARI2    : VARIABLES INTERNES INSTANT +
!
!      SORTIE :
!-------------
! OUT  RADIA    : INDICATEUR DE PERTE DE RADIALITE
! OUT  COSANG   : COSINUS DE L'ANGLE
!
! ......................................................................
!
    integer(kind=8) :: mxcmel
    parameter(mxcmel=162)
!
    integer(kind=8) :: i, k, igau, icine, nbvar, memo, visc, iradi, idelta
!
    real(kind=8) :: dsigma(mxcmel), zero, deux, s1dsig, norm, dnorm, matel(20)
    real(kind=8) :: zernor, tensm(6), tensp(6), indm, indp, xm(6), xp(6)
    real(kind=8) :: coef, cinf, c2inf, mat(50)
    character(len=16) :: relaComp2
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    zero = 0.0d0
    deux = 2.0d0
    zernor = 10.0d0*r8prem()
!
    if (ind .eq. 0) then
!
! ----    CALCUL DE DSIGMA = SIG2 - SIG1 :
!         ----------------------------------
        k = 0
        do igau = 1, npg
            do i = 1, nbsig
                k = k+1
                dsigma(k) = sig2(k)-sig1(k)
!
            end do
        end do
!
! ----    CALCUL DE L'INDICATEUR LOCAL DE PERTE DE RADIALITE
! ----    AUX POINTS D'INTEGRATION :
!         ------------------------
        do igau = 1, npg
!
! ----       CALCUL DU PRODUIT SIG1:(SIG2-SIG1) :
!            ----------------------------------------
            s1dsig = zero
            do i = 1, 3
                s1dsig = s1dsig+sig1(i+(igau-1)*nbsig)*dsigma(i+(igau-1)*nbsig)
            end do
!
            do i = 4, nbsig
                s1dsig = s1dsig+deux*sig1(i+(igau-1)*nbsig)*dsigma(i+(igau-1)*nbsig)
            end do
!
! ----       CALCUL DU SECOND INVARIANT DES TENSEURS DES CONTRAINTES :
!            -------------------------------------------------------
            norm = norsig(sig1(1+(igau-1)*nbsig), nbsig)
            dnorm = norsig(dsigma(1+(igau-1)*nbsig), nbsig)
!
! ----       DANS LE CAS OU NORME(SIG1) = 0  OU NORME(DSIGMA) = 0 :
! ----       ON MET L'INDICATEUR A 0 :
!            -----------------------
            if (norm .le. zernor .or. dnorm .le. zernor) then
                radia(igau) = zero
                cosang(igau) = zero
            else if (dnorm .le. 1.0d4*r8prem()*norm) then
                radia(igau) = zero
                cosang(igau) = zero
            else
                radia(igau) = 1.d0-abs(s1dsig)/norm/dnorm
                cosang(igau) = s1dsig/norm/dnorm
            end if
        end do
!
    else if (ind .eq. 1) then
!
        do igau = 1, npg
!
            iradi = 0
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, sig1(1+(igau-1)*nbsig), b_incx, tensm, b_incy)
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, sig2(1+(igau-1)*nbsig), b_incx, tensp, b_incy)
            b_n = to_blas_int(nbsig-3)
            b_incx = to_blas_int(1)
            call dscal(b_n, sqrt(2.d0), tensm(4), b_incx)
            b_n = to_blas_int(nbsig-3)
            b_incx = to_blas_int(1)
            call dscal(b_n, sqrt(2.d0), tensp(4), b_incx)
!
!           ISOTROPE : LA NORMALE NE DEPEND QUE DE SIG
            if ((relaComp .eq. 'VMIS_ISOT_TRAC') .or. (relaComp .eq. 'VMIS_ISOT_LINE') .or. &
                (relaComp .eq. 'VMIS_ISOT_PUIS')) then
                indm = vari1((igau-1)*nvi+2)
                indp = vari2((igau-1)*nvi+2)
                icine = 0
                iradi = 1
!
!           CINEMATIQUE : LA NORMALE DEPEND DE SIG ET X
            elseif ((relaComp .eq. 'VMIS_ECMI_TRAC') .or. ( &
                    relaComp .eq. 'VMIS_ECMI_LINE')) then
                b_n = to_blas_int(nbsig)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vari1((igau-1)*nvi+3), b_incx, xm, b_incy)
                b_n = to_blas_int(nbsig)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vari2((igau-1)*nvi+3), b_incx, xp, b_incy)
                indm = vari1((igau-1)*nvi+2)
                indp = vari2((igau-1)*nvi+2)
                icine = 1
                iradi = 1
                b_n = to_blas_int(nbsig-3)
                b_incx = to_blas_int(1)
                call dscal(b_n, sqrt(2.d0), xm(4), b_incx)
                b_n = to_blas_int(nbsig-3)
                b_incx = to_blas_int(1)
                call dscal(b_n, sqrt(2.d0), xp(4), b_incx)
!
            else if ((relaComp .eq. 'VMIS_CINE_LINE')) then
                b_n = to_blas_int(nbsig)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vari1((igau-1)*nvi+1), b_incx, xm, b_incy)
                b_n = to_blas_int(nbsig)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vari2((igau-1)*nvi+1), b_incx, xp, b_incy)
                indm = vari1((igau-1)*nvi+7)
                indp = vari2((igau-1)*nvi+7)
                icine = 1
                iradi = 1
                b_n = to_blas_int(nbsig-3)
                b_incx = to_blas_int(1)
                call dscal(b_n, sqrt(2.d0), xm(4), b_incx)
                b_n = to_blas_int(nbsig-3)
                b_incx = to_blas_int(1)
                call dscal(b_n, sqrt(2.d0), xp(4), b_incx)
!
            elseif ((relaComp .eq. 'VMIS_CIN1_CHAB') .or. ( &
                    relaComp .eq. 'VISC_CIN1_CHAB') .or. ( &
                    relaComp .eq. 'VMIS_CIN2_CHAB') .or. ( &
                    relaComp .eq. 'VMIS_CIN2_MEMO') .or. ( &
                    relaComp .eq. 'VISC_CIN2_CHAB') .or. ( &
                    relaComp .eq. 'VISC_CIN2_MEMO')) then
                relaComp2 = relaComp
                call nmcham('RIGI', igau, 1, imate, relaComp2, &
                            matel, mat, nbvar, memo, visc, &
                            idelta, coef)
!              approximation : on supose C constant
                cinf = mat(4)/1.5d0
                indm = vari1((igau-1)*nvi+2)
                indp = vari2((igau-1)*nvi+2)
                b_n = to_blas_int(nbsig)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vari1((igau-1)*nvi+3), b_incx, xm, b_incy)
                b_n = to_blas_int(nbsig)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vari2((igau-1)*nvi+3), b_incx, xp, b_incy)
                b_n = to_blas_int(nbsig)
                b_incx = to_blas_int(1)
                call dscal(b_n, cinf, xm, b_incx)
                b_n = to_blas_int(nbsig)
                b_incx = to_blas_int(1)
                call dscal(b_n, cinf, xp, b_incx)
                if (nbvar .eq. 2) then
                    c2inf = mat(9)/1.5d0
                    b_n = to_blas_int(nbsig)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, c2inf, vari1((igau-1)*nvi+9), b_incx, xm, &
                               b_incy)
                    b_n = to_blas_int(nbsig)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, c2inf, vari2((igau-1)*nvi+9), b_incx, xp, &
                               b_incy)
                end if
                icine = 1
                iradi = 1
                b_n = to_blas_int(nbsig-3)
                b_incx = to_blas_int(1)
                call dscal(b_n, sqrt(2.d0), xm(4), b_incx)
                b_n = to_blas_int(nbsig-3)
                b_incx = to_blas_int(1)
                call dscal(b_n, sqrt(2.d0), xp(4), b_incx)
!
!
            end if
!
!           CALCUL EFFECTUE UNIQUEMENT SI LE COMPORTEMENT LE PERMET
            if (iradi .eq. 1) then
                call radial(nbsig, tensm, tensp, indm, indp, &
                            icine, xm, xp, radia(igau))
                cosang(igau) = sqrt(abs(1.d0-radia(igau)*radia(igau)))
            else
                radia(igau) = 0.d0
                cosang(igau) = 0.d0
            end if
!
        end do
!
    end if
end subroutine
