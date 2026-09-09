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
module FE_mass_module
!
    use FE_basis_module
    use FE_quadrature_module
    use HHO_utils_module, only: hhoCopySymPartMat
!
    implicit none
!
    private
#include "asterf_types.h"
#include "FE_module.h"
#include "blas/dsyr.h"
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
    public :: FEMassMatScal
!    private  ::
!
contains
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEMassMatScal(FEQuad, FEBasis, mass, ValuesQP)
!
        implicit none
!
        type(FE_Quadrature), intent(in) :: FEQuad
        type(FE_Basis), intent(in) :: FEBasis
        real(kind=8), intent(out) :: mass(MAX_BS_CG, MAX_BS_CG)
        real(kind=8), intent(in), optional :: ValuesQP(MAX_QP)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the mass matrix
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
        blas_int :: b_incx, b_lda, b_n
!
        mass = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%func(FEQuad%points_param(1:3, ipg))
!
! ---- mass = mass + weight * BSEval^T * BSEVAL
            if (present(ValuesQP)) then
                coeff = FEQuad%weights(ipg)*ValuesQP(ipg)
!
            else
                coeff = FEQuad%weights(ipg)
            end if
            b_n = to_blas_int(FEBasis%size)
            b_incx = to_blas_int(1)
            b_lda = to_blas_int(MAX_BS_CG)
            call dsyr('U', b_n, coeff, BSEval, b_incx, &
                      mass, b_lda)
        end do
!
! ----- Copy the lower part
!
        call hhoCopySymPartMat('U', mass(1:FEBasis%size, 1:FEBasis%size))
!
    end subroutine
!
end module
