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
!
subroutine nmfici(nno, nddl, wref, vff, dfde, &
                  geom, poids, b)
!
!
    implicit none
#include "asterfort/subaco.h"
#include "asterfort/sumetr.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer(kind=8) :: nno, nddl
    real(kind=8) :: wref, vff(nno), dfde(2, nno), geom(3, 2*nno)
    real(kind=8) :: poids, b(3, 3, nddl/3)
!-----------------------------------------------------------------------
!  MATRICE CINEMATIQUE POUR LES JOINTS 3D (EN UN POINT DE GAUSS DONNE)
!  ATTENTION : ON UTILISE LES FONCTIONS DE FORME DES FACE DU JOINT 3D,
!              PAR EXEMPLE FONCTION DE FORME DU QUAD8 POUR LES HEXA20.
!-----------------------------------------------------------------------
! IN  NNO  NOMBRE DE NOEUDS DE LA FACE (2*NNO POUR TOUT L'ELEM LINEAIRE)
! IN  NDDL   NOMBRE DE DEGRES DE LIBERTE EN DEPL TOTAL (3 PAR NOEUDS)
! IN  WREF   POIDS DE REFERENCE DU POINT DE GAUSS
! IN  VFF    VALEUR DES FONCTIONS DE FORME (DE LA FACE)
!            VVF A POUR DIM NNO, Y COMPRIS EN QUADRATIQUE
! IN  DFDE   DERIVEE DES FONCTIONS DE FORME (DE LA FACE)
! IN  GEOM   COORDONNEES DES NOEUDS
! OUT POIDS  POIDS REEL DU POINT DE GAUSS (AVEC DISTORSION)
! OUT B      MATRICE DE PASSAGE UNODAL -> SAUT DE U LOCAL
!-----------------------------------------------------------------------
    integer(kind=8) :: n
    real(kind=8) :: cova(3, 3), metr(2, 2), jac, r(3, 3), noa1
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!
!
!    CALCUL DE LA BASE COVARIANTE LOCALE ET JACOBIEN
    call subaco(nno, dfde, geom(:, 1:nno), cova)
    call sumetr(cova, metr, jac)
    poids = wref*jac
!
!    CALCUL DE LA BASE ORTHONORMEE LOCALE
    noa1 = sqrt(cova(1, 1)**2+cova(2, 1)**2+cova(3, 1)**2)
    r(1, 1) = cova(1, 3)
    r(1, 2) = cova(2, 3)
    r(1, 3) = cova(3, 3)
    r(2, 1) = cova(1, 1)/noa1
    r(2, 2) = cova(2, 1)/noa1
    r(2, 3) = cova(3, 1)/noa1
    r(3, 1) = r(1, 2)*r(2, 3)-r(1, 3)*r(2, 2)
    r(3, 2) = r(1, 3)*r(2, 1)-r(1, 1)*r(2, 3)
    r(3, 3) = r(1, 1)*r(2, 2)-r(1, 2)*r(2, 1)
!
!
!    CONSTRUCTION DE LA MATRICE B
!
    do n = 1, nno
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, r, b_incx, b(1, 1, n), b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, -vff(n), b(1, 1, n), b_incx)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, r, b_incx, b(1, 1, n+nno), b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, vff(n), b(1, 1, n+nno), b_incx)
    end do
!
end subroutine
