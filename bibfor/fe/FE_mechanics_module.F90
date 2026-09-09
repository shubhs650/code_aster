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
module FE_mechanics_module
!
    use FE_basis_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "FE_module.h"
#include "asterfort/assert.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE - Finite Element
!
! Module to compute mechanical terms
!
! --------------------------------------------------------------------------------------------------
!
    public :: FEMatFB, FEMatB, matG2F, matG2E, FEMatBB, matG2Epsi
!    private  ::
!
contains
!
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEMatFB(FEBasis, point, BGSEval, f, matFB)
!
        implicit none
!
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: point(3)
        real(kind=8), intent(in), dimension(3, MAX_BS_CG) :: BGSEval
        real(kind=8), intent(out) :: matFB(6, MAX_BS_CG, 3)
        real(kind=8), intent(in) :: f(3, 3)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the matrix [F].[B]
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
!
        integer(kind=8) :: i, i_dim
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
        real(kind=8) :: funcEF(MAX_BS_CG), r
!
        matFB = 0.d0
!
        select case (FEBasis%ndim)
        case (2)
            do i = 1, FEBasis%size
                do i_dim = 1, 2
                    matFB(1, i, i_dim) = f(i_dim, 1)*BGSEval(1, i)
                    matFB(2, i, i_dim) = f(i_dim, 2)*BGSEval(2, i)
                    matFB(4, i, i_dim) = (f(i_dim, 1)*BGSEval(2, i)+f(i_dim, 2)*BGSEval(1, i))/rac2
                end do
            end do
        case (3)
            do i = 1, FEBasis%size
                do i_dim = 1, 3
                    matFB(1, i, i_dim) = f(i_dim, 1)*BGSEval(1, i)
                    matFB(2, i, i_dim) = f(i_dim, 2)*BGSEval(2, i)
                    matFB(3, i, i_dim) = f(i_dim, 3)*BGSEval(3, i)
                    matFB(4, i, i_dim) = (f(i_dim, 1)*BGSEval(2, i)+f(i_dim, 2)*BGSEval(1, i))/rac2
                    matFB(5, i, i_dim) = (f(i_dim, 1)*BGSEval(3, i)+f(i_dim, 3)*BGSEval(1, i))/rac2
                    matFB(6, i, i_dim) = (f(i_dim, 2)*BGSEval(3, i)+f(i_dim, 3)*BGSEval(2, i))/rac2
                end do
            end do
        case default
            ASSERT(ASTER_FALSE)
        end select
!
        if (FEBasis%l_axis) then
            funcEF = FEBasis%func(point)
            r = 0.d0
            do i = 1, FEBasis%size
                r = r+funcEF(i)*FEBasis%coorno(1, i)
            end do
            do i = 1, FEBasis%size
                matFB(3, i, 1) = f(3, 3)*funcEF(i)/r
            end do
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEMatBB(FEBasis, BGSEval, matBB)
!
        implicit none
!
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in), dimension(3, MAX_BS_CG) :: BGSEval
        real(kind=8), intent(out) :: matBB(6, MAX_BS_CG, MAX_BS_CG)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the matrix [B]^T.[B]
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
!
        integer(kind=8) :: i, j
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
!
        matBB = 0.d0
!
        select case (FEBasis%ndim)
        case (2)
            do i = 1, FEBasis%size
                do j = 1, i
                    matBB(1, i, j) = BGSEval(1, i)*BGSEval(1, j)
                    matBB(2, i, j) = BGSEval(2, i)*BGSEval(2, j)
                    matBB(4, i, j) = (BGSEval(1, i)*BGSEval(2, j)+BGSEval(2, i)*BGSEval(1, j))/rac2
!
                    matBB(1, j, i) = matBB(1, i, j)
                    matBB(2, j, i) = matBB(2, i, j)
                    matBB(4, j, i) = matBB(4, i, j)
                end do
            end do
        case (3)
            do i = 1, FEBasis%size
                do j = 1, i
                    matBB(1, i, j) = BGSEval(1, i)*BGSEval(1, j)
                    matBB(2, i, j) = BGSEval(2, i)*BGSEval(2, j)
                    matBB(3, i, j) = BGSEval(3, i)*BGSEval(3, j)
                    matBB(4, i, j) = (BGSEval(1, i)*BGSEval(2, j)+BGSEval(2, i)*BGSEval(1, j))/rac2
                    matBB(5, i, j) = (BGSEval(1, i)*BGSEval(3, j)+BGSEval(3, i)*BGSEval(1, j))/rac2
                    matBB(6, i, j) = (BGSEval(2, i)*BGSEval(3, j)+BGSEval(3, i)*BGSEval(2, j))/rac2
!
                    matBB(1, j, i) = matBB(1, i, j)
                    matBB(2, j, i) = matBB(2, i, j)
                    matBB(3, j, i) = matBB(3, i, j)
                    matBB(4, j, i) = matBB(4, i, j)
                    matBB(5, j, i) = matBB(5, i, j)
                    matBB(6, j, i) = matBB(6, i, j)
                end do
            end do
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEMatB(FEBasis, point, BGSEval, matB)
!
        implicit none
!
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: point(3)
        real(kind=8), intent(in), dimension(3, MAX_BS_CG) :: BGSEval
        real(kind=8), intent(out) :: matB(6, MAX_BS_CG, 3)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity vector
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, BGSEval v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
!
        integer(kind=8) :: i
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
        real(kind=8) :: funcEF(MAX_BS_CG), r
!
        matB = 0.d0
!
        select case (FEBasis%ndim)
        case (2)
            do i = 1, FEBasis%size
                matB(1, i, 1) = BGSEval(1, i)
                matB(2, i, 2) = BGSEval(2, i)
                matB(4, i, 1) = BGSEval(2, i)/rac2
                matB(4, i, 2) = BGSEval(1, i)/rac2
            end do
        case (3)
            do i = 1, FEBasis%size
                matB(1, i, 1) = BGSEval(1, i)
                matB(2, i, 2) = BGSEval(2, i)
                matB(3, i, 3) = BGSEval(3, i)
                matB(4, i, 1) = BGSEval(2, i)/rac2
                matB(5, i, 1) = BGSEval(3, i)/rac2
                matB(4, i, 2) = BGSEval(1, i)/rac2
                matB(6, i, 2) = BGSEval(3, i)/rac2
                matB(5, i, 3) = BGSEval(1, i)/rac2
                matB(6, i, 3) = BGSEval(2, i)/rac2
            end do
        case default
            ASSERT(ASTER_FALSE)
        end select
!
        if (FEBasis%l_axis) then
            funcEF = FEBasis%func(point)
            r = 0.d0
            do i = 1, FEBasis%size
                r = r+funcEF(i)*FEBasis%coorno(1, i)
            end do
            do i = 1, FEBasis%size
                matB(3, i, 1) = funcEF(i)/r
            end do
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    function matG2F(grad) result(f)
!
        implicit none
!
        real(kind=8), intent(in) :: grad(3, 3)
        real(kind=8) :: f(3, 3)
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate F = I + G
!
! --------------------------------------------------------------------------------------------------
!
!
        f = grad
        f(1, 1) = f(1, 1)+1.d0
        f(2, 2) = f(2, 2)+1.d0
        f(3, 3) = f(3, 3)+1.d0
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function matG2E(grad) result(e)
!
        implicit none
!
        real(kind=8), intent(in) :: grad(3, 3)
        real(kind=8) :: e(6)
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate E = 1/2(C-I)
!
! --------------------------------------------------------------------------------------------------
!
!
        integer(kind=8) :: i, j, k
        real(kind=8) :: c(3, 3), tmp
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
!
        do i = 1, 3
            do j = 1, i
                tmp = grad(i, j)+grad(j, i)
!
                do k = 1, 3
                    tmp = tmp+grad(k, i)*grad(k, j)
                end do
!
                c(i, j) = 0.5d0*tmp
!
            end do
        end do
!
        e(1) = c(1, 1)
        e(2) = c(2, 2)
        e(3) = c(3, 3)
        e(4) = c(2, 1)*rac2
        e(5) = c(3, 1)*rac2
        e(6) = c(3, 2)*rac2
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function matG2Epsi(grad) result(e)
!
        implicit none
!
        real(kind=8), intent(in) :: grad(3, 3)
        real(kind=8) :: e(6)
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate Epsi = 1/2(G^T+G)
!
! --------------------------------------------------------------------------------------------------
!
!
        real(kind=8), parameter :: rac2_2 = sqrt(2.d0)/2.d0
!
        e(1) = grad(1, 1)
        e(2) = grad(2, 2)
        e(3) = grad(3, 3)
        e(4) = (grad(1, 2)+grad(2, 1))*rac2_2
        e(5) = (grad(1, 3)+grad(3, 1))*rac2_2
        e(6) = (grad(2, 3)+grad(3, 2))*rac2_2
!
    end function
!
end module
