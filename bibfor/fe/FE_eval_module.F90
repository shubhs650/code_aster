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
module FE_eval_module
!
    use FE_basis_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "FE_module.h"
#include "blas/ddot.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE - Finite Element
!
! Module to eval function
!
! --------------------------------------------------------------------------------------------------
!
    public :: FEEvalFuncRScal, FEEvalGradVec, FEEvalGradMat, FEEvalGradSymMat
    public :: FEEvalFuncRVec, FEEvalFuncFVec, FEEvalFuncFScal
!    private  ::
!
contains
!
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalFuncRScal(FEBasis, val_nodes, point) result(func)
!
        implicit none
!
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: val_nodes(*)
        real(kind=8), intent(in) :: point(3)
        real(kind=8) :: func
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate scalar values from value at nodes
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        real(kind=8) :: funcEF(MAX_BS_CG)
        blas_int :: b_incx, b_incy, b_n
!
        funcEF = FEBasis%func(point)
        b_n = to_blas_int(FEBasis%size)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        func = ddot(b_n, val_nodes, b_incx, funcEF, b_incy)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalFuncFScal(nomfunc, nbpara, nompara, valpara) result(func)
!
        implicit none
!
        character(len=8), intent(in) :: nomfunc
        integer(kind=8), intent(in) :: nbpara
        character(len=8), intent(in) :: nompara(*)
        real(kind=8), intent(inout) :: valpara(*)
        real(kind=8) :: func
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate scalar values from value at nodes
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iret
!
        call fointe('FM', nomfunc, nbpara, nompara, valpara, func, iret)
        ASSERT(iret == 0)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalFuncRVec(FEBasis, val_nodes, point) result(func)
!
        implicit none
!
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: val_nodes(*)
        real(kind=8), intent(in) :: point(3)
        real(kind=8) :: func(3)
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate vector values from value at nodes
!   In FEBasis     : Basis function
!   In ValuesQP    : Values of vector function f at the quadrature points
!   Out rhs        : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        real(kind=8) :: funcEF(MAX_BS_CG)
        integer(kind=8) :: idim
        blas_int :: b_incx, b_incy, b_n
!
        func = 0.d0
        funcEF = FEBasis%func(point)
        b_n = to_blas_int(FEBasis%size)
        b_incx = to_blas_int(FEBasis%ndim)
        b_incy = to_blas_int(1)
        do idim = 1, FEBasis%ndim
            func(idim) = ddot(b_n, val_nodes(idim), b_incx, funcEF, b_incy)
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalFuncFVec(nomfunc, nbpara, nompara, valpara, ndim) result(func)
!
        implicit none
!
        character(len=8), intent(in) :: nomfunc(*)
        integer(kind=8), intent(in) :: nbpara, ndim
        character(len=8), intent(in) :: nompara(*)
        real(kind=8), intent(inout) :: valpara(*)
        real(kind=8) :: func(3)
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate scalar values from value at nodes
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iret, idim
!
        func = 0.d0
        do idim = 1, ndim
            call fointe('FM', nomfunc(idim), nbpara, nompara, valpara, func(idim), iret)
            ASSERT(iret == 0)
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalGradVec(FEBasis, val_nodes, point, BGSEval) result(grad)
!
        implicit none
!
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: val_nodes(*)
        real(kind=8), intent(in) :: point(3)
        real(kind=8) :: grad(3)
        real(kind=8), intent(in), optional :: BGSEval(3, MAX_BS_CG)
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate scalar values from value at nodes
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: i
        real(kind=8) :: gradEF(3, MAX_BS_CG)
!
        grad = 0.d0
        if (present(BGSEval)) then
            do i = 1, FEBasis%size
                grad = grad+val_nodes(i)*BGSEval(1:3, i)
            end do
        else
            gradEF = FEBasis%grad(point)
            do i = 1, FEBasis%size
                grad = grad+val_nodes(i)*gradEF(1:3, i)
            end do
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalGradSymMat(FEBasis, val_nodes, point, BGSEval) result(grads)
!
        implicit none
!
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: val_nodes(FEBasis%ndim, *)
        real(kind=8), intent(in) :: point(3)
        real(kind=8) :: grads(6)
        real(kind=8), intent(in), optional :: BGSEval(3, MAX_BS_CG)
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate scalar values from value at nodes
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        real(kind=8) :: grad(3, 3)
        real(kind=8), parameter :: rac2_2 = sqrt(2.d0)/2.d0
!
        if (present(BGSEval)) then
            grad = FEEvalGradMat(FEBasis, val_nodes, point, BGSEval)
        else
            grad = FEEvalGradMat(FEBasis, val_nodes, point)
        end if
!
        grads = 0.d0
        select case (FEBasis%ndim)
        case (2)
            grads(1) = grad(1, 1)
            grads(2) = grad(2, 2)
            grads(3) = grad(3, 3)
            grads(4) = (grad(2, 1)+grad(1, 2))*rac2_2
        case (3)
            grads(1) = grad(1, 1)
            grads(2) = grad(2, 2)
            grads(3) = grad(3, 3)
            grads(4) = (grad(1, 2)+grad(2, 1))*rac2_2
            grads(5) = (grad(1, 3)+grad(3, 1))*rac2_2
            grads(6) = (grad(2, 3)+grad(3, 2))*rac2_2
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end function
!
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalGradMat(FEBasis, val_nodes, point, BGSEval) result(grad)
!
        implicit none
!
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: val_nodes(*)
        real(kind=8), intent(in) :: point(3)
        real(kind=8) :: grad(3, 3)
        real(kind=8), intent(in), optional :: BGSEval(3, MAX_BS_CG)
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate scalar values from value at nodes
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: i, n, ind
        real(kind=8) :: gradEF(3, MAX_BS_CG), funcEF(MAX_BS_CG), ur, r
!
        grad = 0.d0
        if (present(BGSEval)) then
            ind = 0
            select case (FEBasis%ndim)
            case (2)
                do n = 1, FEBasis%size
                    do i = 1, 2
                        grad(i, 1:2) = grad(i, 1:2)+BGSEval(1:2, n)*val_nodes(ind+i)
                    end do
                    ind = ind+2
                end do
            case (3)
                do n = 1, FEBasis%size
                    do i = 1, 3
                        grad(i, 1:3) = grad(i, 1:3)+BGSEval(1:3, n)*val_nodes(ind+i)
                    end do
                    ind = ind+3
                end do
            case default
                ASSERT(ASTER_FALSE)
            end select
        else
            gradEF = FEBasis%grad(point)
            ind = 0
            select case (FEBasis%ndim)
            case (2)
                do n = 1, FEBasis%size
                    do i = 1, 2
                        grad(i, 1:2) = grad(i, 1:2)+gradEF(1:2, n)*val_nodes(ind+i)
                    end do
                    ind = ind+2
                end do
            case (3)
                do n = 1, FEBasis%size
                    do i = 1, 3
                        grad(i, 1:3) = grad(i, 1:3)+gradEF(1:3, n)*val_nodes(ind+i)
                    end do
                    ind = ind+3
                end do
            case default
                ASSERT(ASTER_FALSE)
            end select
        end if
!
        if (FEBasis%l_axis) then
            funcEF = FEBasis%func(point)
            r = 0.d0
            ur = 0.d0
            do n = 1, FEBasis%size
                r = r+funcEF(n)*FEBasis%coorno(1, n)
                ur = ur+funcEF(n)*val_nodes(2*(n-1)+1)
            end do
            grad(3, 3) = ur/r
        end if
!
    end function
!
end module
