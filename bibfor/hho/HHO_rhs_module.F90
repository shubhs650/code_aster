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
module HHO_rhs_module
!
    use HHO_basis_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use FE_algebra_module
!
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/HHO_size_module.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Module to compute the rhs member (f,v)
!
! --------------------------------------------------------------------------------------------------
!
    public :: hhoMakeRhsFaceScal, hhoMakeRhsFaceVec
    public :: hhoMakeRhsCellScal, hhoMakeRhsCellVec
    private :: hhoMakeRhsFaceScalBase, hhoMakeRhsCellScalBase
!
contains
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsFaceScalBase(hhoQuad, hhoBasisFace, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Quadrature), intent(in) :: hhoQuad
        type(HHO_basis_face), intent(inout) :: hhoBasisFace
        real(kind=8), intent(in) :: ValuesQP(MSIZE_QP_FACE)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_FACE_SCAL)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vF)_F
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of scalar function f at the quadrature points
!   In degree       : degree of the basis
!   Out rhs         : (f, vF)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: ipg, size
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: BSFEval
        real(kind=8) :: coeff
!
        rhs = 0.d0
!
! -- init face basis
        size = hhoBasisFace%BSSize(0, degree)
!
! -- Loop on quadrature point
        do ipg = 1, hhoQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), 0, degree, BSFEval)
!
! ---- rhs = rhs + weight * BSFEval
            coeff = hhoQuad%weights(ipg)*ValuesQP(ipg)
            call daxpy_1(size, coeff, BSFEval, rhs)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsFaceScal(hhoFace, hhoQuad, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: ValuesQP(MSIZE_QP_FACE)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_FACE_SCAL)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vF)_F
!   In hhoFace      : the current HHO Face
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of scalar function f at the quadrature points
!   In degree       : degree of the basis
!   Out rhs         : (f, vF)_F term
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_face) :: hhoBasisFace
!
! -- init face basis
        call hhoBasisFace%initialize(hhoFace)
!
        call hhoMakeRhsFaceScalBase(hhoQuad, hhoBasisFace, ValuesQP, degree, rhs)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsFaceVec(hhoFace, hhoQuad, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: ValuesQP(3, MSIZE_QP_FACE)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_FACE_VEC)
!
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vF)_F
!   In hhoFace      : the current HHO Face
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of vectorial function f at the quadrature points
!   In degree  : degree of the basis
!   Out rhs  : term (f, vF)_F (rhs member)
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_face) :: hhoBasisFace
        integer(kind=8) :: size, idir, begin, i
        real(kind=8) :: Values(MSIZE_QP_FACE), rhs_dir(MSIZE_FACE_SCAL)
!
! -- init face basis
        call hhoBasisFace%initialize(hhoFace)
        size = hhoBasisFace%BSSize(0, degree)
!
        rhs = 0.d0
        Values = 0.d0
!
        begin = 1
        do idir = 1, hhoFace%ndim+1
            do i = 1, hhoQuad%nbQuadPoints
                Values(i) = ValuesQP(idir, i)
            end do
            call hhoMakeRhsFaceScalBase(hhoQuad, hhoBasisFace, Values, degree, rhs_dir)
            call dcopy_1(size, rhs_dir, rhs(begin))
            begin = begin+size
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsCellScalBase(hhoQuad, hhoBasisCell, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Quadrature), intent(in) :: hhoQuad
        type(HHO_basis_cell), intent(inout) :: hhoBasisCell
        real(kind=8), intent(in) :: ValuesQP(MSIZE_QP_CELL)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_CELL_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vT)_T
!   In hhoCell      : the current HHO Cell
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of scalar function f at the quadrature points
!   Out rhs         : (f, vT)_T term
!   In  degree      : degree of the basis
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: ipg, size
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: BSCEval
        real(kind=8) :: coeff
!
! -- init face basis
        size = hhoBasisCell%BSSize(0, degree)
!
        rhs = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, hhoQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, degree, BSCEval)
!
! ---- rhs = rhs + weight * BSCEval
            coeff = hhoQuad%weights(ipg)*ValuesQP(ipg)
            call daxpy_1(size, coeff, BSCEval, rhs)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsCellScal(hhoCell, hhoQuad, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: ValuesQP(MSIZE_QP_CELL)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_CELL_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vT)_T
!   In hhoCell      : the current HHO Cell
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of scalar function f at the quadrature points
!   Out rhs         : (f, vT)_T term
!   In  degree      : degree of the basis
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        type(HHO_basis_cell) :: hhoBasisCell
!
! -- init face basis
        call hhoBasisCell%initialize(hhoCell)
!
        call hhoMakeRhsCellScalBase(hhoQuad, hhoBasisCell, ValuesQP, degree, rhs)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsCellVec(hhoCell, hhoQuad, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: ValuesQP(3, MSIZE_QP_CELL)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_CELL_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vT)_T
!   In hhoCell      : the current HHO Cell
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of vectorial function f at the quadrature points
!   Out rhs         : term (f, vT)_T (rhs member)
!   In degree  : degree of the basis
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        integer(kind=8) :: size, idir, begin, i
        real(kind=8) :: Values(MSIZE_QP_CELL), rhs_dir(MSIZE_CELL_SCAL)
!
! -- init face basis
        call hhoBasisCell%initialize(hhoCell)
        size = hhoBasisCell%BSSize(0, degree)
!
        rhs = 0.d0
        Values = 0.d0
!
        begin = 1
        do idir = 1, hhoCell%ndim
            do i = 1, hhoQuad%nbQuadPoints
                Values(i) = ValuesQP(idir, i)
            end do
            call hhoMakeRhsCellScalBase(hhoQuad, hhoBasisCell, Values, degree, rhs_dir)
            call dcopy_1(size, rhs_dir, rhs(begin))
            begin = begin+size
        end do
!
    end subroutine
!
end module
