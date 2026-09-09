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
!
module HHO_matrix_module
!
    use FE_algebra_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/readMatrix.h"
#include "asterfort/writeMatrix.h"
#include "blas/daxpy.h"
#include "blas/dgemv.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - dynamic matrix
!
! --------------------------------------------------------------------------------------------------
!

    type HHO_matrix
        integer(kind=8) :: nrows = 0, ncols = 0
        integer(kind=8) :: max_nrows = 0, max_ncols = 0
        aster_logical :: is_allocated = ASTER_FALSE
! ----- array
        real(kind=8), dimension(:, :), pointer :: m
!
! ----- member function
    contains
        procedure, pass :: initialize => hhoMatriceInit
        procedure, pass :: free => hhoMatriceFree
        procedure, pass :: write => hhoMatriceWrite
        procedure, pass :: read => hhoMatriceRead
        procedure, pass :: setValue => hhoMatriceSetValue
        procedure, pass :: print => hhoMatricePrint
        procedure, pass :: copySymU => hhoMatriceCopySymU
        procedure, pass :: copy => hhoMatriceCopy
        procedure, pass :: add => hhoMatriceAdd
        procedure, pass :: addBlock => hhoMatriceAddSub
        procedure, pass :: addBlock2 => hhoMatriceAddSub2
        procedure, pass :: prune => hhoMatricePrune
        procedure, pass :: dot => hhoMatriceDot
        procedure, pass :: scale => hhoMatriceScale
!
    end type HHO_matrix
!
    public   :: HHO_matrix
    private  :: hhoMatriceInit, hhoMatriceFree, hhoMatriceWrite, hhoMatriceSetValue
    private  :: hhoMatriceRead, hhoMatricePrint, hhoMatriceCopySymU, hhoMatriceCopy
    private  :: hhoMatriceAdd, hhoMatricePrune, hhoMatriceDot
    private  :: hhoMatriceAddSub2, hhoMatriceAddSub, hhoMatriceScale
!
contains
!---------------------------------------------------------------------------------------------------
! -- member function for HHO_matrix type
!---------------------------------------------------------------------------------------------------
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceInit(this, n_rows, n_cols, val)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        integer(kind=8), intent(in) :: n_rows, n_cols
        real(kind=8), intent(in), optional :: val
!
        ASSERT(.not. this%is_allocated)
!
        this%nrows = n_rows
        this%ncols = n_cols
        this%max_nrows = n_rows
        this%max_ncols = n_cols
!
        ASSERT(n_rows > 0 .and. n_cols > 0)
!
        allocate (this%m(n_rows, n_cols))
        this%is_allocated = ASTER_TRUE
!
        if (present(val)) then
            call this%setValue(val)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceFree(this)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
!
        this%nrows = 0
        this%ncols = 0

        if (this%is_allocated) then
            deallocate (this%m)
            this%is_allocated = ASTER_FALSE
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceWrite(this, name, l_sym)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
        character(len=*), intent(in) :: name
        aster_logical, intent(in) :: l_sym
!
        ASSERT(this%is_allocated)
        call writeMatrix(name, this%nrows, this%ncols, l_sym, this%m)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceRead(this, name, l_sym)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
        character(len=*), intent(in) :: name
        aster_logical, intent(in) :: l_sym
!
        ASSERT(this%is_allocated)
        call readMatrix(name, this%nrows, this%ncols, l_sym, this%m)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceSetValue(this, val)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        real(kind=8), intent(in) :: val
!
        ASSERT(this%is_allocated)
        this%m = val
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatricePrint(this)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i
!
!
        write (6, *) "matrix of", this%nrows, "rows x ", this%ncols, "cols"
        do i = 1, this%nrows
            write (6, '(50ES13.6)') this%m(i, 1:this%ncols)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceCopySymU(this, row_offset_, col_offset_)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
        integer(kind=8), intent(in), optional :: row_offset_, col_offset_
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: row_offset, col_offset, il, size, ig, jg
!
        row_offset = 0
        if (present(row_offset_)) row_offset = row_offset_
        col_offset = 0
        if (present(col_offset_)) col_offset = col_offset_
!
        ASSERT(this%ncols-col_offset == this%nrows-row_offset)
        size = this%nrows-row_offset
        do il = 1, size-1
            ig = row_offset+il
            jg = col_offset+il
            this%m(ig+1:this%nrows, jg) = this%m(ig, jg+1:this%ncols)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceCopy(this, mat, row_offset_, col_offset_)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        type(HHO_matrix), intent(in) :: mat
        integer(kind=8), intent(in), optional :: row_offset_, col_offset_
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: row_offset, col_offset
!
        row_offset = 0
        if (present(row_offset_)) row_offset = row_offset_
        col_offset = 0
        if (present(col_offset_)) col_offset = col_offset_
!
        ASSERT(this%nrows >= row_offset+mat%nrows)
        ASSERT(this%ncols >= col_offset+mat%ncols)
!
        this%m(row_offset+1:row_offset+1+mat%nrows, col_offset+1:col_offset+1+mat%ncols) = mat%m
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceAdd(this, mat, alpha_)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        type(HHO_matrix), intent(in) :: mat
        real(kind=8), intent(in), optional :: alpha_
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: j
        real(kind=8) :: alpha
!
        alpha = 1.d0
        if (present(alpha_)) alpha = alpha_
!
        ASSERT(this%nrows == mat%nrows)
        ASSERT(this%ncols == mat%ncols)
!
        do j = 1, this%ncols
            call daxpy_1(this%nrows, alpha, mat%m(:, j), this%m(:, j))
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceAddSub(this, mat, row_offset, col_offset, alpha_)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        type(HHO_matrix), intent(in) :: mat
        integer(kind=8), intent(in) :: row_offset, col_offset
        real(kind=8), intent(in), optional :: alpha_
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: j, j2
        real(kind=8) :: alpha
!
        alpha = 1.d0
        if (present(alpha_)) alpha = alpha_
!
        ASSERT(this%nrows >= mat%nrows+row_offset)
        ASSERT(this%ncols >= mat%ncols+col_offset)
!
        do j = 1, mat%ncols
            j2 = col_offset+j
            call daxpy_1(mat%nrows, alpha, mat%m(:, j), this%m(row_offset+1:, j2))
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceAddSub2(this, mat, max_row, nb_row, nb_col, row_offset, col_offset, alpha_)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        integer(kind=8), intent(in) :: nb_row, nb_col, max_row
        real(kind=8), intent(in) :: mat(max_row, nb_col)
        integer(kind=8), intent(in) :: row_offset, col_offset
        real(kind=8), intent(in), optional :: alpha_
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: j, j2
        real(kind=8) :: alpha
!
        alpha = 1.d0
        if (present(alpha_)) alpha = alpha_
!
        ASSERT(this%nrows >= nb_row+row_offset)
        ASSERT(this%ncols >= nb_col+col_offset)
!
        do j = 1, nb_col
            j2 = col_offset+j
            call daxpy_1(nb_row, alpha, mat(:, j), this%m(row_offset+1:, j2))
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceDot(this, x, y, alpha_, beta_)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(in), optional :: alpha_, beta_
        real(kind=8), intent(inout) :: y(*)
!
! --------------------------------------------------------------------------------------------------
!
!   y = beta*y+ alpha*A*x
! --------------------------------------------------------------------------------------------------
!
        blas_int :: b_n, b_m, b_lda
        blas_int, parameter :: one = to_blas_int(1)
        real(kind=8) :: alpha, beta
!
        alpha = 1.d0
        beta = 0.d0
!
        if (present(alpha_)) then
            alpha = alpha_
        end if
        if (present(beta_)) then
            beta = alpha_
        end if
!
        b_lda = to_blas_int(this%max_nrows)
        b_m = to_blas_int(this%nrows)
        b_n = to_blas_int(this%ncols)
        call dgemv('N', b_m, b_n, alpha, this%m, &
                   b_lda, x, one, beta, y, one)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatricePrune(this, threshold)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        real(kind=8), intent(in) :: threshold
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, j
!
        do j = 1, this%ncols
            do i = 1, this%nrows
                if (abs(this%m(i, j)) < threshold) then
                    this%m(i, j) = 0.d0
                end if
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceScale(this, val)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        real(kind=8), intent(in) :: val
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, j
!
        do j = 1, this%ncols
            do i = 1, this%nrows
                this%m(i, j) = val*this%m(i, j)
            end do
        end do
!
    end subroutine
!
end module
