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
module HHO_statcond_module
!
    use HHO_type
    use HHO_size_module
    use HHO_utils_module
    use HHO_matrix_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/utmess.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
#include "blas/dpotrf.h"
#include "blas/dpotrs.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Static condensation
!
! Routine to compute static condensation or decondensation
!
! --------------------------------------------------------------------------------------------------
    public :: hhoCondStatic
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCondStatic(cbs, lhs_elem, rhs_elem, &
                             lhs_cond, rhs_cond, lhs_decond, rhs_decond)
!
        implicit none
!
        integer(kind=8), intent(in) :: cbs
        type(HHO_matrix), intent(in) :: lhs_elem
        type(HHO_matrix), intent(inout) :: lhs_cond, lhs_decond
        real(kind=8), dimension(MSIZE_TDOFS_VEC), intent(in)  :: rhs_elem
        real(kind=8), dimension(MSIZE_TDOFS_VEC), intent(inout)  :: rhs_cond, rhs_decond
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the static condensation - symmetric matrice only
!   In hhoCell      : the current HHO Cell
!   In hhoDta       : information on HHO methods
!   In lhs          : lhs with cell and faces terms (is symmetric)
!   In rhs          : rhs with cell and faces terms
!   In l_lhs_sym    : lhs is symmetric ?
!   Out lhs_local   : lhs after static condensation (is symmetric)
!   Out rhs_local   : rhs after static condensation
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: faces_dofs, total_dofs, i
        real(kind=8) :: rhs_T(MSIZE_CELL_VEC), rhs_F(MSIZE_FDOFS_VEC)
        type(HHO_matrix) :: K_TT, K_FT, K_TF, K_FF
        blas_int :: info, b_n, b_lda, b_nrhs, b_ldb, b_ldc, b_m, b_k
        blas_int, parameter :: b_one = to_blas_int(1)
!
! ---- Number of dofs
        total_dofs = lhs_elem%nrows
        faces_dofs = total_dofs-cbs
!
        call K_TT%initialize(cbs, cbs)
        call K_TF%initialize(cbs, faces_dofs)
        call K_FT%initialize(faces_dofs, cbs)
        call K_FF%initialize(faces_dofs, faces_dofs)
!
        K_TT%m(1:cbs, 1:cbs) = lhs_elem%m(faces_dofs+1:total_dofs, faces_dofs+1:total_dofs)
        K_TF%m(1:cbs, 1:faces_dofs) = lhs_elem%m(faces_dofs+1:total_dofs, 1:faces_dofs)
        K_FF%m(1:faces_dofs, 1:faces_dofs) = lhs_elem%m(1:faces_dofs, 1:faces_dofs)
        K_FT%m(1:faces_dofs, 1:cbs) = lhs_elem%m(1:faces_dofs, faces_dofs+1:total_dofs)
        rhs_T(1:cbs) = rhs_elem(faces_dofs+1:total_dofs)
        rhs_F(1:faces_dofs) = rhs_elem(1:faces_dofs)
!
! ---- factorize K_TT
        info = 0
        b_n = to_blas_int(K_TT%nrows)
        b_lda = to_blas_int(K_TT%max_nrows)
        call dpotrf('U', b_n, K_TT%m, b_lda, info)
!
! ---- Sucess ?
        if (info .ne. 0) then
            call utmess('E', 'HHO1_4')
        end if
!
! ---- Solve K_TF = K_TT^-1 * K_TF
        info = 0
        b_n = to_blas_int(K_TT%nrows)
        b_nrhs = to_blas_int(K_TF%ncols)
        b_lda = to_blas_int(K_TT%max_nrows)
        b_ldb = to_blas_int(K_TF%max_nrows)
        call dpotrs('U', b_n, b_nrhs, K_TT%m, b_lda, K_TF%m, b_ldb, info)
!
! ---- Sucess ?
        if (info .ne. 0) then
            call utmess('F', 'HHO1_4')
        end if
!
! ---- Solve rhs_T = K_TT^-1 * rhs
        info = 0
        b_n = to_blas_int(K_TT%nrows)
        b_lda = to_blas_int(K_TT%max_nrows)
        b_ldb = to_blas_int(MSIZE_CELL_VEC)
        call dpotrs('U', b_n, b_one, K_TT%m, b_lda, rhs_T, b_ldb, info)
!
! ---- Sucess ?
        if (info .ne. 0) then
            call utmess('F', 'HHO1_4')
        end if
!
! ---- Compute K_FF = K_FF - K_FT * K_TF
!
        b_ldc = to_blas_int(K_FF%max_nrows)
        b_ldb = to_blas_int(K_TF%max_nrows)
        b_lda = to_blas_int(K_FT%max_nrows)
        b_m = to_blas_int(K_FT%nrows)
        b_n = to_blas_int(K_TF%ncols)
        b_k = to_blas_int(K_TF%nrows)
        call dgemm('N', 'N', b_m, b_n, b_k, -1.d0, K_FT%m, b_lda, &
                    & K_TF%m, b_ldb, 1.d0, K_FF%m, b_ldc)
!
! ---- Compute rhs_F = rhs_F - K_FT * rhs_T
        b_lda = to_blas_int(K_FT%max_nrows)
        b_m = to_blas_int(K_FT%nrows)
        b_n = to_blas_int(K_FT%ncols)
        call dgemv('N', b_m, b_n, -1.d0, K_FT%m, b_lda, rhs_T, b_one, 1.d0, rhs_F, b_one)
!
! ---- local matrix and vector
        call lhs_cond%initialize(total_dofs, total_dofs, 0.d0)
        lhs_cond%m(1:faces_dofs, 1:faces_dofs) = K_FF%m(1:faces_dofs, 1:faces_dofs)
!
        do i = 1, cbs
            lhs_cond%m(faces_dofs+i, faces_dofs+i) = 1.d0
        end do
!
        call lhs_decond%initialize(total_dofs, total_dofs, 0.d0)
        lhs_decond%m(faces_dofs+1:total_dofs, 1:faces_dofs) = -K_TF%m(1:cbs, 1:faces_dofs)
!
        rhs_cond = 0.d0
        rhs_cond(1:faces_dofs) = rhs_F(1:faces_dofs)
!
        rhs_decond = 0.d0
        rhs_decond(faces_dofs+1:total_dofs) = rhs_T(1:cbs)
!
        call K_TT%free()
        call K_TF%free()
        call K_FT%free()
        call K_FF%free()
!
    end subroutine
!
end module
