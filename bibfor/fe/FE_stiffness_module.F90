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
module FE_stiffness_module
!
    use FE_basis_module
    use FE_quadrature_module
    use FE_algebra_module
    use HHO_utils_module, only: hhoCopySymPartMat
!
! aslint: disable=C1505
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
! Module to compute stiffness terms
!
! --------------------------------------------------------------------------------------------------
!
    public :: FEStiffResiScal, FEStiffJacoScal, FEStiffResiScalAdd, FEStiffJacoScalAdd
    public :: FEMassStiffJacoScalAdd
    public :: FEStiffResiVectSymAdd, FEStiffJacoVectSymAdd, FEStiffGeomVectSymAdd
!    private  ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffResiScalAdd(FEBasis, BGSEval, weight, ValuesQP, vec)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(3, MAX_BS_CG)  :: BGSEval
        real(kind=8), intent(in)            :: weight
        real(kind=8), intent(inout)         :: vec(MAX_BS_CG)
        real(kind=8), intent(in)            :: ValuesQP(3)
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Compute the rigidity vector
!   In FEQuad      : Quadrature
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
!
!
        real(kind=8) :: val_w(3)
!
        val_w = weight*ValuesQP
        if (FEBasis%ndim == 3) then
            call dgemv_T_3xn(BGSEval, FEBasis%size, val_w, vec)
        elseif (FEBasis%ndim == 2) then
            call dgemv_T_2xn(BGSEval, FEBasis%size, val_w, vec)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffResiScal(FEQuad, FEBasis, ValuesQP, vec)
!
        implicit none
!
        type(FE_Quadrature), intent(in)     :: FEQuad
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(out)           :: vec(MAX_BS_CG)
        real(kind=8), intent(in)            :: ValuesQP(3, MAX_QP)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity vector
!   In FEQuad      : Quadrature
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: ipg
        real(kind=8), dimension(3, MAX_BS_CG) :: BSEval
!
        vec = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%grad(FEQuad%points_param(1:3, ipg), FEQuad%jacob(1:3, 1:3, ipg))
!
            call FEStiffResiScalAdd(FEBasis, BSEval, FEQuad%weights(ipg), ValuesQP(1:3, ipg), vec)
        end do
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffJacoScalAdd(FEBasis, BGSEval, weight, ValueQP, mat)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(3, MAX_BS_CG) :: BGSEval
        real(kind=8), intent(in)            :: weight
        real(kind=8), intent(inout)         :: mat(MAX_BS_CG, MAX_BS_CG)
        real(kind=8), intent(in)            :: ValueQP(3, 3)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity matrix
!   In FEQuad      : Quadrature
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: j
        real(kind=8) :: Kgradj(3)
!
        do j = 1, FEBasis%size
            if (FEBasis%ndim == 3) then
                call dgemv_3x3(ValueQP, BGSEval(1, j), Kgradj, weight)
                call dgemv_T_3xn(BGSEval, FEBasis%size, Kgradj, mat(1, j))
            elseif (FEBasis%ndim == 2) then
                call dgemv_2x2(ValueQP, BGSEval(:, j), Kgradj, weight)
                call dgemv_T_2xn(BGSEval, FEBasis%size, Kgradj, mat(:, j))
            else
                ASSERT(ASTER_FALSE)
            end if
        end do
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffJacoScal(FEQuad, FEBasis, ValuesQP, mat)
!
        implicit none
!
        type(FE_Quadrature), intent(in)     :: FEQuad
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(out)           :: mat(MAX_BS_CG, MAX_BS_CG)
        real(kind=8), intent(in)            :: ValuesQP(3, 3, MAX_QP)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity matrix
!   In FEQuad      : Quadrature
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: ipg
        real(kind=8), dimension(3, MAX_BS_CG) :: BSEval
!
        mat = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%grad(FEQuad%points_param(1:3, ipg), FEQuad%jacob(1:3, 1:3, ipg))

            call FEStiffJacoScalAdd(FEBasis, BSEval, FEQuad%weights(ipg), ValuesQP(1:3, 1:3, ipg), &
                                    mat)
!
        end do
!
! ----- Copy the lower part
!
        call hhoCopySymPartMat('U', mat(1:FEBasis%size, 1:FEBasis%size))
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffResiVectSymAdd(FEBasis, def, weight, stress, vec)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(6, MAX_BS_CG, 3) :: def
        real(kind=8), intent(in)            :: weight
        real(kind=8), intent(inout)         :: vec(*)
        real(kind=8), intent(in)            :: stress(6)
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Compute the rigidity vector (with symetric stess)
!   In FEQuad      : Quadrature
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: stress_w(6)
!
        stress_w = weight*stress
        select case (FEBasis%ndim)
        case (2)
            call dgemv_T_4xn(def(1, 1, 1), FEBasis%size, stress_w, vec(1), 2)
            call dgemv_T_4xn(def(1, 1, 2), FEBasis%size, stress_w, vec(2), 2)
        case (3)
            call dgemv_T_6xn(def(1, 1, 1), FEBasis%size, stress_w, vec(1), 3)
            call dgemv_T_6xn(def(1, 1, 2), FEBasis%size, stress_w, vec(2), 3)
            call dgemv_T_6xn(def(1, 1, 3), FEBasis%size, stress_w, vec(3), 3)
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffJacoVectSymAdd(FEBasis, def, weight, dsidep, l_matsym, mat)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(6, MAX_BS_CG, 3) :: def
        real(kind=8), intent(in)            :: weight
        aster_logical, intent(in)           :: l_matsym
        real(kind=8), intent(inout)         :: mat(*)
        real(kind=8), intent(in)            :: dsidep(6, 6)
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Compute the rigidity vector (with symetric stess)
!   In FEQuad      : Quadrature
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_n, i_d, j_d, kkd
        real(kind=8) :: sig_w(6)
!
        select case (FEBasis%ndim)
        case (2)
            if (l_matsym) then
                do i_n = 1, FEBasis%size
                    do i_d = 1, 2
                        kkd = (2*(i_n-1)+i_d-1)*(2*(i_n-1)+i_d)/2
                        call dgemv_T_4x4(dsidep, def(1, i_n, i_d), sig_w, weight)
                        call dgemv_T_4xn(def(1, 1, 1), i_n-1, sig_w, mat(kkd+1), 2)
                        call dgemv_T_4xn(def(1, 1, 2), i_n-1, sig_w, mat(kkd+2), 2)

                        ! j_n = i_n
                        kkd = kkd+2*(i_n-1)
                        do j_d = 1, i_d
                            mat(kkd+j_d) = mat(kkd+j_d)+ &
                                           def(1, i_n, j_d)*sig_w(1)+def(2, i_n, j_d)*sig_w(2)+ &
                                           def(3, i_n, j_d)*sig_w(3)+def(4, i_n, j_d)*sig_w(4)
                        end do
                    end do
                end do
            else
                do i_n = 1, FEBasis%size
                    do i_d = 1, 2
                        kkd = 2*FEBasis%size*(2*(i_n-1)+i_d-1)
                        call dgemv_T_4x4(dsidep, def(1, i_n, i_d), sig_w, weight)
                        call dgemv_T_4xn(def(1, 1, 1), FEBasis%size, sig_w, mat(kkd+1), 2)
                        call dgemv_T_4xn(def(1, 1, 2), FEBasis%size, sig_w, mat(kkd+2), 2)
                    end do
                end do
            end if
        case (3)
            if (l_matsym) then
                do i_n = 1, FEBasis%size
                    do i_d = 1, 3
                        kkd = (3*(i_n-1)+i_d-1)*(3*(i_n-1)+i_d)/2
                        call dgemv_T_6x6(dsidep, def(1, i_n, i_d), sig_w, weight)
                        call dgemv_T_6xn(def(1, 1, 1), i_n-1, sig_w, mat(kkd+1), 3)
                        call dgemv_T_6xn(def(1, 1, 2), i_n-1, sig_w, mat(kkd+2), 3)
                        call dgemv_T_6xn(def(1, 1, 3), i_n-1, sig_w, mat(kkd+3), 3)

                        ! j_n = i_n
                        kkd = kkd+3*(i_n-1)
                        do j_d = 1, i_d
                            mat(kkd+j_d) = mat(kkd+j_d)+ &
                                           def(1, i_n, j_d)*sig_w(1)+def(2, i_n, j_d)*sig_w(2)+ &
                                           def(3, i_n, j_d)*sig_w(3)+def(4, i_n, j_d)*sig_w(4)+ &
                                           def(5, i_n, j_d)*sig_w(5)+def(6, i_n, j_d)*sig_w(6)
                        end do
                    end do
                end do
            else
                do i_n = 1, FEBasis%size
                    do i_d = 1, 3
                        kkd = 3*FEBasis%size*(3*(i_n-1)+i_d-1)
                        call dgemv_T_6x6(dsidep, def(1, i_n, i_d), sig_w, weight)
                        call dgemv_T_6xn(def(1, 1, 1), FEBasis%size, sig_w, mat(kkd+1), 3)
                        call dgemv_T_6xn(def(1, 1, 2), FEBasis%size, sig_w, mat(kkd+2), 3)
                        call dgemv_T_6xn(def(1, 1, 3), FEBasis%size, sig_w, mat(kkd+3), 3)
                    end do
                end do
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffGeomVectSymAdd(FEBasis, pff, weight, stress, l_matsym, mat)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(6, MAX_BS_CG, MAX_BS_CG) :: pff
        real(kind=8), intent(in)            :: weight
        aster_logical, intent(in)           :: l_matsym
        real(kind=8), intent(inout)         :: mat(*)
        real(kind=8), intent(in)            :: stress(6)
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Compute the rigidity vector (with symetric stess)
!   In FEQuad      : Quadrature
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_n, kkd, j_n, i_d
        real(kind=8) :: stress_w(6), tmp(MAX_BS_CG)
!
        stress_w = weight*stress
!
        select case (FEBasis%ndim)
        case (2)
            if (l_matsym) then
                do i_n = 1, FEBasis%size

                    tmp = 0.d0
                    call dgemv_T_4xn(pff(1, 1, i_n), i_n, stress_w, tmp, 1)
                    do i_d = 1, 2
                        kkd = (2*(i_n-1)+i_d-1)*(2*(i_n-1)+i_d)/2+i_d
                        do j_n = 1, i_n
                            mat(kkd) = mat(kkd)+tmp(j_n)
                            kkd = kkd+2
                        end do
                    end do
                end do
            else
                do i_n = 1, FEBasis%size
                    tmp = 0.d0
                    call dgemv_T_4xn(pff(1, 1, i_n), FEBasis%size, stress_w, tmp, 1)
                    do i_d = 1, 2
                        kkd = 2*FEBasis%size*(2*(i_n-1)+i_d-1)+i_d
                        do j_n = 1, FEBasis%size
                            mat(kkd) = mat(kkd)+tmp(j_n)
                            kkd = kkd+2
                        end do
                    end do
                end do
            end if
        case (3)
            if (l_matsym) then
                do i_n = 1, FEBasis%size
                    tmp = 0.d0
                    call dgemv_T_6xn(pff(1, 1, i_n), i_n, stress_w, tmp, 1)
                    do i_d = 1, 3
                        kkd = (3*(i_n-1)+i_d-1)*(3*(i_n-1)+i_d)/2+i_d
                        do j_n = 1, i_n
                            mat(kkd) = mat(kkd)+tmp(j_n)
                            kkd = kkd+3
                        end do
                    end do
                end do
            else
                do i_n = 1, FEBasis%size
                    tmp = 0.d0
                    call dgemv_T_6xn(pff(1, 1, i_n), FEBasis%size, stress_w, tmp, 1)
                    do i_d = 1, 3
                        kkd = 3*FEBasis%size*(3*(i_n-1)+i_d-1)+i_d
                        do j_n = 1, FEBasis%size
                            mat(kkd) = mat(kkd)+tmp(j_n)
                            kkd = kkd+3
                        end do
                    end do
                end do
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
    subroutine FEMassStiffJacoScalAdd(BSEval, BGSEval, weight, ValueQP, mat)
!
        implicit none
!
        real(kind=8), intent(in)    :: BSEval(MAX_BS_CG)
        real(kind=8), intent(in)    :: weight, BGSEval(3, MAX_BS_CG)
        real(kind=8), intent(in)    :: ValueQP(3)
        real(kind=8), intent(inout) :: mat(MAX_BS_CG, MAX_BS_CG)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity matrix
!   In FEQuad      : Quadrature
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        real(kind=8) :: gloDt(3, MAX_BS_CG)

        gloDt(1, :) = ValueQP(1)*BSEval
        gloDt(2, :) = ValueQP(2)*BSEval
        gloDt(3, :) = ValueQP(3)*BSEval
        mat = mat+weight*matmul(transpose(BGSEval), gloDt)

!
    end subroutine
!
end module
