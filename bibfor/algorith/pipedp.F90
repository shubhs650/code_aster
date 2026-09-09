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
subroutine pipedp(BEHInteg, ndim, typmod, &
                  epsm, sigm, vim, epsp, &
                  epsd, a0, a1)
!
    use Behaviour_type
    implicit none
!
#include "asterc/matfpe.h"
#include "asterf_types.h"
#include "asterfort/betfpp.h"
#include "asterfort/betmat.h"
#include "asterfort/zerop2.h"
#include "blas/ddot.h"
#include "blas/dnrm2.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=8) :: typmod(2)
    integer(kind=8) :: ndim
    real(kind=8) :: epsp(6), epsd(6)
    real(kind=8) :: epsm(6), vim(2), sigm(6), a0, a1
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE - PRED_ELAS)
!
! LOI DE COMPORTEMENT BETON_DOUBLE_DP
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  TYPMOD : TYPE DE MODELISATION
! IN  MATE   : MATERIAU CODE
! IN  VIM    : VARIABLES INTERNES EN T-
! IN  EPSM   : DEFORMATIONS EN T-
! IN  SIGM   : CONTRAINTES EN T-
! IN  EPSP   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES FIXES
! IN  EPSD   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES PILOTEES
! IN  KPG    : NUMERO DU POINT DE GAUSS
! IN  KPG    : NUMERO DU SOUS-POINT DE GAUSS
! IN  ELGEOM : TABLEAUX DES ELEMENTS GEOMETRIQUES SPECIFIQUES AUX
!               LOIS DE COMPORTEMENT (DIMENSION MAXIMALE FIXEE EN DUR)
! OUT A0     : LINEARISATION DU CRITERE : FEL = A0 + A1*ETA
! OUT A1     : IDEM A0
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: un = 1.d0, deux = 2.d0, trois = 3.d0, neuf = 9.d0
    real(kind=8), parameter :: d23 = .66666666666666d0, d13 = .33333333333333d0
    real(kind=8), parameter :: raci2 = sqrt(2.d0)
    integer(kind=8) :: ndimsi, k, nrac1, nrac2
    aster_logical :: trac, comp, notrac, nocomp
    real(kind=8) :: trsigp, trsigd, sigelp(6), sigeld(6)
    real(kind=8) :: eps1(6), eps2(6), pp(6), dd(6)
    real(kind=8) :: d1, d2, g1, g2, g3, g4
    real(kind=8) :: p0, p1, p2, q0, q1, q2, eta, eta1, eta2
    real(kind=8) :: rac1(2), rac2(2)
    real(kind=8) :: e, nu, lambda, deuxmu
    real(kind=8) :: fc, ft, beta
    real(kind=8) :: a, b, c, d
    integer(kind=8) :: ndt, ndi, nr, nvi
    integer(kind=8), parameter :: nmat = 90
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: pc, pt, kuc, kut, ke, tbid, rbid, fcp, ftp
    character(len=8) :: mod
    character(len=3) :: matcst
    blas_int :: b_incx, b_incy, b_n
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
!
! --------------------------------------------------------------------------------------------------
!
    call matfpe(-1)

! - OPTION ET MODELISATION
    ndimsi = 2*ndim

    tbid = 0.d0
    mod = typmod(1)
    call betmat(BEHInteg%materPara%schemePara%fami, &
                BEHInteg%materPara%schemePara%kpg, &
                BEHInteg%materPara%schemePara%ksp, &
                mod, BEHInteg%materPara%jvMaterCode, &
                nmat, tbid, tbid, materd, materf, &
                matcst, ndt, ndi, nr, nvi)
!
    e = materd(1, 1)
    nu = materd(2, 1)
    beta = materd(3, 2)
!
    pc = vim(1)
    pt = vim(2)
    call betfpp(BEHInteg, materf, nmat, pc, pt, &
                3, fc, ft, rbid, rbid, &
                kuc, kut, ke)
!
    fcp = materf(1, 2)
    ftp = materf(2, 2)
    lambda = e*nu/(1.d0+nu)/(1.d0-2.d0*nu)
    deuxmu = e/(1.d0+nu)
    a = raci2*(beta-un)/(deux*beta-un)
    b = raci2/trois*beta/(deux*beta-un)
    c = raci2
    d = deux*raci2/trois
!
    nocomp = .false.
    if (pc .gt. kuc) nocomp = .true.
    notrac = .false.
    if (pt .gt. kut) notrac = .true.
!
! ======================================================================
!                CALCUL DES DEFORMATIONS POUR LINEARISATION
! ======================================================================
!
!     COEFFICIENTS DE LA FORME QUADRATIQUE DU CRITERE
    do k = 1, ndimsi
        sigelp(k) = sigm(k)+(lambda*kron(k)+deuxmu)*(epsp(k)-epsm(k))
        sigeld(k) = (lambda*kron(k)+deuxmu)*epsd(k)
    end do
!
    trsigp = sigelp(1)+sigelp(2)+sigelp(3)
    trsigd = sigeld(1)+sigeld(2)+sigeld(3)
!
    do k = 1, ndimsi
        pp(k) = sigelp(k)-d13*trsigp*kron(k)
        dd(k) = sigeld(k)-d13*trsigd*kron(k)
    end do
!
!     CRITERE DE TRACTION
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p0 = d13*ddot(b_n, pp, b_incx, pp, b_incy)-(d*ft)**2+d23*c*d*ft*trsigp-c**2/neuf*trsigp**2
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p1 = d13*ddot(b_n, pp, b_incx, dd, b_incy)+d13*c*d*ft*trsigd-c**2/neuf*trsigp*trsigd
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p2 = d13*ddot(b_n, dd, b_incx, dd, b_incy)-c**2/neuf*trsigd**2
!
    p0 = p0/(d*ftp)**2
    p1 = p1/(d*ftp)**2
    p2 = p2/(d*ftp)**2
!
!    CRITERE DE COMPRESSION
!
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    q0 = d13*ddot(b_n, pp, b_incx, pp, b_incy)-(b*fc)**2+d23*a*b*fc*trsigp-a**2/neuf*trsigp**2
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    q1 = d13*ddot(b_n, pp, b_incx, dd, b_incy)+d13*a*b*fc*trsigd-a**2/neuf*trsigp*trsigd
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    q2 = d13*ddot(b_n, dd, b_incx, dd, b_incy)-a**2/neuf*trsigd**2
!
    q0 = q0/(b*fcp)**2
    q1 = q1/(b*fcp)**2
    q2 = q2/(b*fcp)**2
!
!    RECHERCHE DES INTERSECTIONS ELLIPSE / DROITE (TRACTION)
    call zerop2(2*p1/p2, p0/p2, rac1, nrac1)
!
!     RECHERCHE DES INTERSECTIONS ELLIPSE / DROITE (COMPRESSION)
    call zerop2(2*q1/q2, q0/q2, rac2, nrac2)
!
!     CALCUL DES COEFFICIENTS LORSQUE LES DEUX CRITERES
!     SONT COMPLETEMENT ECROUIS
!     ---------------------------------------------------------
    if (notrac .and. nocomp) then
        a0 = 0.d0
        a1 = 0.d0
    end if
!
!     CALCUL DES COEFFICIENTS LORSQUE LE CRITERE DE TRACTION
!     EST COMPLETEMENT ECROUI
!     ---------------------------------------------------------
    if (notrac) then
        if (nrac2 .ne. 0) then
            do k = 1, ndimsi
                eps1(k) = epsp(k)+rac2(1)*epsd(k)-epsm(k)
                eps2(k) = epsp(k)+rac2(2)*epsd(k)-epsm(k)
            end do
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            d1 = dnrm2(b_n, eps1, b_incx)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            d2 = dnrm2(b_n, eps2, b_incx)
            if (d1 .le. d2) then
                eta = rac2(1)
            else
                eta = rac2(2)
            end if
        else
            eta = -q1/q2
        end if
        a0 = -q2*eta**2+q0
        a1 = 2*(eta*q2+q1)
    end if
!
!     CALCUL DES COEFFICIENTS LORSQUE LE CRITERE DE COMPRESSION
!     EST COMPLETEMENT ECROUI
!     ---------------------------------------------------------
    if (nocomp) then
        if (nrac1 .ne. 0) then
            do k = 1, ndimsi
                eps1(k) = epsp(k)+rac1(1)*epsd(k)-epsm(k)
                eps2(k) = epsp(k)+rac1(2)*epsd(k)-epsm(k)
            end do
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            d1 = dnrm2(b_n, eps1, b_incx)
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            d2 = dnrm2(b_n, eps2, b_incx)
            if (d1 .le. d2) then
                eta = rac1(1)
            else
                eta = rac1(2)
            end if
        else
            eta = -p1/p2
        end if
        a0 = -p2*eta**2+p0
        a1 = 2*(eta*p2+p1)
    end if
!
!     CALCUL DES COEFFICIENTS LORSQU'AUCUN DES DEUX CRITERES
!     N'EST COMPLETEMENT ECROUI
!     ---------------------------------------------------------
!
    if (.not. nocomp .and. .not. notrac) then
!        LA DROITE COUPE LE CRITERE DE TRACTION
        if (nrac1 .ne. 0) then
            g1 = q0+deux*rac1(1)*q1+rac1(1)**2*q2
            g2 = q0+deux*rac1(2)*q1+rac1(2)**2*q2
            trac = .true.
            if ((g1 .lt. 0) .and. (g2 .lt. 0)) then
                do k = 1, ndimsi
                    eps1(k) = epsp(k)+rac1(1)*epsd(k)-epsm(k)
                    eps2(k) = epsp(k)+rac1(2)*epsd(k)-epsm(k)
                end do
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                d1 = dnrm2(b_n, eps1, b_incx)
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                d2 = dnrm2(b_n, eps2, b_incx)
                if (d1 .le. d2) then
                    eta1 = rac1(1)
                else
                    eta1 = rac1(2)
                end if
            else if ((g1 .lt. 0) .and. (g2 .gt. 0)) then
                eta1 = rac1(1)
            else if ((g1 .gt. 0) .and. (g2 .lt. 0)) then
                eta1 = rac1(2)
            else
                eta1 = -p1/p2
                trac = .false.
            end if
        else
            eta1 = -p1/p2
            trac = .false.
        end if
!
!        LA DROITE COUPE LE CRITERE DE COMPRESSION
        if (nrac2 .ne. 0) then
            g3 = p0+deux*rac2(1)*p1+rac2(1)**2*p2
            g4 = p0+deux*rac2(2)*p1+rac2(2)**2*p2
            comp = .true.
            if ((g3 .lt. 0) .and. (g4 .lt. 0)) then
                do k = 1, ndimsi
                    eps1(k) = epsp(k)+rac2(1)*epsd(k)-epsm(k)
                    eps2(k) = epsp(k)+rac2(2)*epsd(k)-epsm(k)
                end do
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                d1 = dnrm2(b_n, eps1, b_incx)
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                d2 = dnrm2(b_n, eps2, b_incx)
                if (d1 .le. d2) then
                    eta2 = rac2(1)
                else
                    eta2 = rac2(2)
                end if
            else if ((g3 .lt. 0) .and. (g4 .gt. 0)) then
                eta2 = rac2(1)
            else if ((g3 .gt. 0) .and. (g4 .lt. 0)) then
                eta2 = rac2(2)
            else
                eta2 = -q1/q2
                comp = .false.
            end if
        else
            eta2 = -q1/q2
            comp = .false.
        end if
!
        if (trac) then
            if (comp) then
                do k = 1, ndimsi
                    eps1(k) = epsp(k)+eta1*epsd(k)-epsm(k)
                    eps2(k) = epsp(k)+eta2*epsd(k)-epsm(k)
                end do
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                d1 = dnrm2(b_n, eps1, b_incx)
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                d2 = dnrm2(b_n, eps2, b_incx)
                if (d1 .le. d2) then
                    eta = eta1
                    a0 = -p2*eta**2+p0
                    a1 = 2*(eta*p2+p1)
                else
                    eta = eta2
                    a0 = -q2*eta**2+q0
                    a1 = 2*(eta*q2+q1)
                end if
            else
                eta = eta1
                a0 = -p2*eta**2+p0
                a1 = 2*(eta*p2+p1)
            end if
        else
            if (comp) then
                eta = eta2
                a0 = -q2*eta**2+q0
                a1 = 2*(eta*q2+q1)
            else
                do k = 1, ndimsi
                    eps1(k) = epsp(k)+eta1*epsd(k)-epsm(k)
                    eps2(k) = epsp(k)+eta2*epsd(k)-epsm(k)
                end do
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                d1 = dnrm2(b_n, eps1, b_incx)
                b_n = to_blas_int(ndimsi)
                b_incx = to_blas_int(1)
                d2 = dnrm2(b_n, eps2, b_incx)
                if (d1 .le. d2) then
                    eta = eta1
                    a0 = -p2*eta**2+p0
                    a1 = 2*(eta*p2+p1)
                else
                    eta = eta2
                    a0 = -q2*eta**2+q0
                    a1 = 2*(eta*q2+q1)
                end if
            end if
        end if
    end if
!
    call matfpe(1)
!
end subroutine
