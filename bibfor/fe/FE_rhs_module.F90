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
module FE_rhs_module
!
    use FE_basis_module
    use FE_quadrature_module
    use FE_algebra_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "FE_module.h"
#include "blas/daxpy.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE - Finite Element
!
! Module to compute the rhs member (f,v)
!
! --------------------------------------------------------------------------------------------------
!
    public :: FEMakeRhsScal, FEMakeRhsVec
!    private  ::
!
contains
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FeMakeRhsScal(FEQuad, FEBasis, ValuesQP, rhs)
!
        implicit none
!
        type(FE_Quadrature), intent(in) :: FEQuad
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: ValuesQP(MAX_QP)
        real(kind=8), intent(out) :: rhs(MAX_BS_CG)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, v)_F
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: ipg
        real(kind=8), dimension(MAX_BS_CG) :: BSEval
        real(kind=8) :: coeff
!
        rhs = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%func(FEQuad%points_param(1:3, ipg))
!
! ---- rhs = rhs + weight * BSEval
            coeff = FEQuad%weights(ipg)*ValuesQP(ipg)
            call daxpy_1(FEBasis%size, coeff, BSEval, rhs)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FeMakeRhsVec(FEQuad, FEBasis, ValuesQP, rhs)
!
        implicit none
!
        type(FE_Quadrature), intent(in) :: FEQuad
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: ValuesQP(3, MAX_QP)
        real(kind=8), intent(out) :: rhs(MAX_BV_CG)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, v)_F
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: ipg, idim
        real(kind=8), dimension(MAX_BS_CG) :: BSEval
        real(kind=8) :: coeff(3)
!
        rhs = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%func(FEQuad%points_param(1:3, ipg))
            coeff = FEQuad%weights(ipg)*ValuesQP(1:3, ipg)
!
! ---- rhs = rhs + weight * BSEval
            do idim = 1, FEBasis%ndim
                call daxpy_1xm(FEBasis%size, coeff(idim), BSEval, rhs(idim), FEBasis%ndim)
            end do
        end do
!
    end subroutine
!
end module
