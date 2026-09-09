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

subroutine ngmatr(nddl, neps, npg, w, b, dsidep, matsym, matuu, matns)

    implicit none

#include "asterf_types.h"
#include "blas/dgemm.h"
#include "asterfort/assert.h"

    integer(kind=8), intent(in) :: nddl, neps, npg
    real(kind=8), intent(in) :: w(neps, npg), b(neps, npg, nddl)
    real(kind=8), intent(inout):: dsidep(neps, neps, npg)
    aster_logical, intent(in):: matsym
    real(kind=8), intent(out), optional:: matuu((nddl*(nddl+1))/2)
    real(kind=8), intent(out), target, optional :: matns(nddl, nddl)
!
! --------------------------------------------------------------------------------------------------
!
!     Construction de la matrice Bt.D.B (stockage par lignes successives)
!
! --------------------------------------------------------------------------------------------------
! in  nddl    : nombre de degres de liberte
! in  neps    : nombre de composantes de deformation et contrainte
! in  npg     : nombre de points de gauss
! in  w       : poids des points de gauss
! in  b       : matrice cinematique : deformation = b.ddl
! inout dsisdep : matrices tangentes locales D aux points de Gauss (modifiée dans la routine)
! in  matsym  : construit-on une matrice symétrique ou non
! out matuu   : matrice de rigidite symetrique
! out matns   : matrice de rigidite non symetrique
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: kpg, i, j, nepg
    real(kind=8) :: ktgb(0:neps*npg*nddl-1)
    real(kind=8), pointer, dimension(:, :) :: ktan_t => null()
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
! --------------------------------------------------------------------------------------------------

! - Memory management
    if (matsym) then
        ASSERT(present(matuu))
        allocate (ktan_t(nddl, nddl))
    else
        ASSERT(present(matns))
        ktan_t => matns
    end if

! - Dimensions
    nepg = neps*npg

! - PRISE EN CHARGE DU POIDS DU POINT DE GAUSS  WG.DSIDEP
    do i = 1, neps
        dsidep(:, i, :) = dsidep(:, i, :)*w
    end do

! - CALCUL DES PRODUITS INTERMEDIAIRES (WG.DSIDEP).B POUR CHAQUE G
    do kpg = 1, npg
        b_ldc = to_blas_int(nepg)
        b_ldb = to_blas_int(nepg)
        b_lda = to_blas_int(neps)
        b_m = to_blas_int(neps)
        b_n = to_blas_int(nddl)
        b_k = to_blas_int(neps)
        call dgemm('N', 'N', b_m, b_n, b_k, &
                   1.d0, dsidep(1, 1, kpg), b_lda, b(1, kpg, 1), b_ldb, &
                   0.d0, ktgb((kpg-1)*neps), b_ldc)
    end do

! - CALCUL DU PRODUIT FINAL SOMME(G) BT. ((WG.DSIDEP).B)  TRANSPOSE
    b_ldc = to_blas_int(nddl)
    b_ldb = to_blas_int(nepg)
    b_lda = to_blas_int(nepg)
    b_m = to_blas_int(nddl)
    b_n = to_blas_int(nddl)
    b_k = to_blas_int(nepg)
    call dgemm('T', 'N', b_m, b_n, b_k, &
               1.d0, ktgb, b_lda, b, b_ldb, &
               0.d0, ktan_t, b_ldc)

! - Stockage de la matrice symetrique
    if (matsym) then
        forall (i=1:nddl, j=1:nddl, i .ge. j) matuu((i*(i-1))/2+j) = ktan_t(j, i)
        deallocate (ktan_t)
    end if

end subroutine
