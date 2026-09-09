! --------------------------------------------------------------------
! Copyright (C) 2007 NECS - BRUNO ZUBER   WWW.NECS.FR
! Copyright (C) 2007 - 2026 - EDF - www.code-aster.org
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
! aslint: disable=W1504,W1006
!
subroutine nmfi3d(BEHInteg, typmod, &
                  nno, nddl, npg, lgpg, wref, &
                  vff, dfde, option, geom, &
                  deplm, ddepl, sigm, sigp, fint, &
                  ktan, vim, vip, carcri, compor, &
                  matsym, coopg, tm, tp, lMatr, &
                  lVect, lSigm, codret)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmfici.h"
#include "asterfort/r8inir.h"
#include "blas/ddot.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8) :: nno, nddl, npg, lgpg, codret
    real(kind=8) :: wref(npg), vff(nno, npg), dfde(2, nno, npg)
    real(kind=8) :: geom(nddl), deplm(nddl), ddepl(nddl), tm, tp
    real(kind=8) :: fint(nddl), ktan(*), coopg(4, npg)
    real(kind=8) :: sigm(3, npg), sigp(3, npg), vim(lgpg, npg), vip(lgpg, npg)
    aster_logical :: matsym
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    aster_logical, intent(in) :: lMatr, lVect, lSigm
!
! --------------------------------------------------------------------------------------------------
!
!  OPTIONS DE MECANIQUE NON LINEAIRE POUR LES JOINTS 3D (TE0206)
!
! --------------------------------------------------------------------------------------------------
!
! IN  NNO    NOMBRE DE NOEUDS DE LA FACE (*2 POUR TOUT L'ELEMENT)
! IN  NDDL   NOMBRE DE DEGRES DE LIBERTE EN DEPL TOTAL (3 PAR NOEUDS)
! IN  NPG    NOMBRE DE POINTS DE GAUSS
! IN  LGPG   NOMBRE DE VARIABLES INTERNES
! IN  WREF   POIDS DE REFERENCE DES POINTS DE GAUSS
! IN  VFF    VALEUR DES FONCTIONS DE FORME (DE LA FACE)
! IN  DFDE   DERIVEE DES FONCTIONS DE FORME (DE LA FACE)
! IN  OPTION OPTION DE CALCUL
! IN  GEOM   COORDONNEES DES NOEUDS
! IN  DEPLM  DEPLACEMENTS NODAUX AU DEBUT DU PAS DE TEMPS
! IN  DDEPL  INCREMENT DES DEPLACEMENTS NODAUX
! IN  SIGM   CONTR LOCALES AUX POINTS DE GAUSS - (SIGN, SITX, SITY)
! OUT SIGP   CONTR LOCALES AUX POINTS DE GAUSS + (SIGN, SITX, SITY)
! OUT FINT   FORCES NODALES
! OUT KTAN   MATRICE TANGENTE (STOCKEE EN TENANT COMPTE DE LA SYMETRIE)
! IN  VIM    VARIABLES INTERNES AU DEBUT DU PAS DE TEMPS
! OUT VIP    VARIABLES INTERNES A LA FIN DU PAS DE TEMPS
! IN  CRIT   VALEURS DE L'UTILISATEUR POUR LES CRITERES DE CONVERGENCE
! IN  COMPOR NOM DE LA LOI DE COMPORTEMENT
! IN  MATSYM INFORMATION SUR LA MATRICE TANGENTE : SYMETRIQUE OU PAS
! IN  COOPG  COORDONNEES GEOMETRIQUES DES PG + POIDS
! OUT CODRET CODE RETOUR DE L'INTEGRATION
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: ksp = 1, ndim = 3
    integer(kind=8) :: cod(9), ni, mj, kk, p, q, kpg, n
    real(kind=8) :: b(3, 60), sigmo(6), sigma(6)
    real(kind=8) :: sum(3), dsu(3), dsidep(6, 6), poids
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    sum = 0.d0
    dsu = 0.d0
    cod = 0
    if (lVect) then
        fint = 0.d0
    end if
    if (lMatr) then
        if (matsym) then
            call r8inir((nddl*(nddl+1))/2, 0.d0, ktan, 1)
        else
            call r8inir(nddl*nddl, 0.d0, ktan, 1)
        end if
    end if

! - Loop on Gauss points
    do kpg = 1, npg
!
! CALCUL DE LA MATRICE B DONNANT LES SAUT PAR ELEMENTS A PARTIR DES
! DEPLACEMENTS AUX NOEUDS , AINSI QUE LE POIDS DES PG :
!
        call nmfici(nno, nddl, wref(kpg), vff(1, kpg), dfde(1, 1, kpg), &
                    geom, poids, b)
!
! CALCUL DU SAUT DE DEPLACEMENT - : SUM, ET DE L'INCREMENT : DSU
! AU POINT DE GAUSS KPG
!
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(3)
        b_incy = to_blas_int(1)
        sum(1) = ddot(b_n, b(1, 1), b_incx, deplm, b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(3)
        b_incy = to_blas_int(1)
        sum(2) = ddot(b_n, b(2, 1), b_incx, deplm, b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(3)
        b_incy = to_blas_int(1)
        sum(3) = ddot(b_n, b(3, 1), b_incx, deplm, b_incy)
        if (lVect) then
            b_n = to_blas_int(nddl)
            b_incx = to_blas_int(3)
            b_incy = to_blas_int(1)
            dsu(1) = ddot(b_n, b(1, 1), b_incx, ddepl, b_incy)
            b_n = to_blas_int(nddl)
            b_incx = to_blas_int(3)
            b_incy = to_blas_int(1)
            dsu(2) = ddot(b_n, b(2, 1), b_incx, ddepl, b_incy)
            b_n = to_blas_int(nddl)
            b_incx = to_blas_int(3)
            b_incy = to_blas_int(1)
            dsu(3) = ddot(b_n, b(3, 1), b_incx, ddepl, b_incy)
        end if

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)
        BEHInteg%behavESVA%behavESVAGeom%coorElga(kpg, 1:3) = coopg(1:3, kpg)

! ----- Integrator
        sigmo = 0.d0
        do n = 1, 3
            sigmo(n) = sigm(n, kpg)
        end do
        sigma = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    tm, tp, &
                    compor, carcri, multComp, &
                    3, sum, dsu, &
                    3, sigmo, &
                    vim(1, kpg), &
                    sigma, vip(1, kpg), &
                    36, dsidep, cod(kpg))
        if (cod(kpg) .eq. 1) goto 900

! ----- Stresses
        if (lSigm) then
            do n = 1, 3
                sigp(n, kpg) = sigma(n)
            end do
        end if

! ----- Internal forces
        if (lVect) then
!       Il faudrait séparer les deux => petit travail de réflexion
            ASSERT(lSigm)
            do ni = 1, nddl
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                fint(ni) = fint(ni)+poids*ddot(b_n, b(1, ni), b_incx, sigma, b_incy)
            end do
        end if

! ----- Rigidity matrix
        if (lMatr) then
            if (matsym) then
                kk = 0
                do ni = 1, nddl
                    do mj = 1, ni
                        kk = kk+1
                        do p = 1, 3
                            do q = 1, 3
                                ktan(kk) = ktan(kk)+poids*b(p, ni)*dsidep(p, q)*b(q, mj)
                            end do
                        end do
                    end do
                end do
            else
                kk = 0
                do ni = 1, nddl
                    do mj = 1, nddl
                        kk = kk+1
                        do p = 1, 3
                            do q = 1, 3
                                ktan(kk) = ktan(kk)+poids*b(p, ni)*dsidep(p, q)*b(q, mj)
                            end do
                        end do
                    end do
                end do
            end if
        end if
    end do
!
900 continue
    call codere(cod, npg, codret)
!
end subroutine
