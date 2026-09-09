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
module FE_Basis_module
!
    use FE_topo_module
!
    implicit none
!
    private
#include "asterc/r8gaem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfvf.h"
#include "asterfort/tecael.h"
#include "asterfort/lteatt.h"
#include "asterfort/utmess.h"
#include "FE_module.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "FE_basis_module.h"
! --------------------------------------------------------------------------------------------------
!
! FE - generic
!
! Module to generate basis function used for fe methods
!
! --------------------------------------------------------------------------------------------------
!
    type FE_Basis
!
        integer(kind=8) :: typeEF = EF_LAGRANGE
        integer(kind=8) :: size
! ----- Dimension topologique
        integer(kind=8) :: ndim = 0
        aster_logical :: l_skin = ASTER_FALSE
        aster_logical :: l_axis = ASTER_FALSE
! ----- Type maille
        character(len=8) :: typema = ''
! ----- Nombre de noeuds
        integer(kind=8) :: nbnodes = 0
! ----- Coordonnees des noeuds
        real(kind=8), dimension(3, MT_NNOMAX3D) :: coorno = 0.d0
!
! ----- member function
    contains
        procedure, pass :: initCell => init_cell
        procedure, pass :: initFace => init_face
        procedure, pass :: func => feBSCEval
        procedure, pass :: grad => feBSCGradEv
    end type
!
! --------------------------------------------------------------------------------------------------
    public :: FE_Basis
    private :: feBSCEval, feBSCGradEv, init_cell, init_face, FE_grad_lagr
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine init_cell(this, FECell)
!
        implicit none
!
        class(FE_Basis), intent(inout) :: this
        type(FE_Cell), intent(in) :: FECell
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   Initialization
! --------------------------------------------------------------------------------------------------
!
        this%typema = FECell%typemas
        this%typeEF = EF_LAGRANGE
        this%ndim = FECell%ndim
        this%nbnodes = FECell%nbnodes
        this%coorno = FECell%coorno
        this%l_skin = ASTER_FALSE
!
        if (lteatt('AXIS', 'OUI')) then
            this%l_axis = ASTER_TRUE
        end if
!
        if (this%typeEF == EF_LAGRANGE) then
            call elrfno(this%typema, this%size)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine init_face(this, FESkin)
!
        implicit none
!
        class(FE_Basis), intent(inout) :: this
        type(FE_Skin), intent(in) :: FESkin
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   Initialization
! --------------------------------------------------------------------------------------------------
!
        this%typema = FESkin%typemas
        this%ndim = FESkin%ndim+1
        this%nbnodes = FESkin%nbnodes
        this%coorno(1:3, 1:MT_NNOMAX2D) = FESkin%coorno
        this%typeEF = EF_LAGRANGE
        this%l_skin = ASTER_TRUE
!
        if (this%typeEF == EF_LAGRANGE) then
            call elrfno(this%typema, this%size)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
    !
!===================================================================================================
!
!===================================================================================================
!
    function FE_grad_lagr(this, point, jacob_) result(BSGrad2)
!
        implicit none
!
        class(FE_Basis), intent(in) :: this
        real(kind=8), intent(in) :: point(3)
        real(kind=8), optional, intent(in) :: jacob_(3, 3)
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   Initialization
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, iadzi, iazk24
        real(kind=8) :: jaco(3, 3), cojac(3, 3), jacob
        real(kind=8), dimension(3, MAX_BS_CG) :: BSGrad, BSGrad2
!
        ASSERT(.not. this%l_skin)
        BSGrad = 0.d0
        call elrfdf(this%typema, point, BSGrad)
!
! ---  Compute the jacobienne
        if (present(jacob_)) then
            jaco = jacob_
        else
            jaco = 0.d0
            if (this%ndim == 3) then
                do i = 1, this%nbnodes
                    jaco(1, 1:3) = jaco(1, 1:3)+this%coorno(1:3, i)*BSGrad(1, i)
                    jaco(2, 1:3) = jaco(2, 1:3)+this%coorno(1:3, i)*BSGrad(2, i)
                    jaco(3, 1:3) = jaco(3, 1:3)+this%coorno(1:3, i)*BSGrad(3, i)
                end do
            else if (this%ndim == 2) then
                do i = 1, this%nbnodes
                    jaco(1, 1:2) = jaco(1, 1:2)+this%coorno(1:2, i)*BSGrad(1, i)
                    jaco(2, 1:2) = jaco(2, 1:2)+this%coorno(1:2, i)*BSGrad(2, i)
                end do
                jaco(3, 3) = 1.d0
            else if (this%ndim == 1) then
                do i = 1, this%nbnodes
                    jaco(1, 1) = jaco(1, 1)+this%coorno(1, i)*BSGrad(1, i)
                end do
                jaco(2, 2) = 1.d0
                jaco(3, 3) = 1.d0
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
        cojac(1, 1) = jaco(2, 2)*jaco(3, 3)-jaco(2, 3)*jaco(3, 2)
        cojac(2, 1) = jaco(3, 1)*jaco(2, 3)-jaco(2, 1)*jaco(3, 3)
        cojac(3, 1) = jaco(2, 1)*jaco(3, 2)-jaco(3, 1)*jaco(2, 2)
        cojac(1, 2) = jaco(1, 3)*jaco(3, 2)-jaco(1, 2)*jaco(3, 3)
        cojac(2, 2) = jaco(1, 1)*jaco(3, 3)-jaco(1, 3)*jaco(3, 1)
        cojac(3, 2) = jaco(1, 2)*jaco(3, 1)-jaco(3, 2)*jaco(1, 1)
        cojac(1, 3) = jaco(1, 2)*jaco(2, 3)-jaco(1, 3)*jaco(2, 2)
        cojac(2, 3) = jaco(2, 1)*jaco(1, 3)-jaco(2, 3)*jaco(1, 1)
        cojac(3, 3) = jaco(1, 1)*jaco(2, 2)-jaco(1, 2)*jaco(2, 1)
!
        if (this%ndim == 3) then
            jacob = jaco(1, 1)*cojac(1, 1)+jaco(1, 2)*cojac(2, 1)+jaco(1, 3)*cojac(3, 1)
        else if (this%ndim == 2) then
            jacob = cojac(3, 3)
        else if (this%ndim == 1) then
            jacob = jaco(1, 1)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (abs(jacob) .le. 1.d0/r8gaem()) then
            call tecael(iadzi, iazk24)
            call utmess('F', 'ALGORITH2_59', sk=zk24(iazk24-1+3) (1:8))
        end if
!
        BSGrad2 = 0.d0
        if (this%ndim == 3) then
            cojac = cojac/jacob
            do i = 1, this%size
                BSGrad2(1, i) = ( &
                                cojac(1, 1)*BSGrad(1, i)+cojac(1, 2)*BSGrad(2, i)+cojac(1, 3)*BS&
                                &Grad(3, i) &
                                )
                BSGrad2(2, i) = ( &
                                cojac(2, 1)*BSGrad(1, i)+cojac(2, 2)*BSGrad(2, i)+cojac(2, 3)*BS&
                                &Grad(3, i) &
                                )
                BSGrad2(3, i) = ( &
                                cojac(3, 1)*BSGrad(1, i)+cojac(3, 2)*BSGrad(2, i)+cojac(3, 3)*BS&
                                &Grad(3, i) &
                                )
            end do
        else if (this%ndim == 2) then
            cojac(1:2, 1:2) = cojac(1:2, 1:2)/jacob
            do i = 1, this%size
                BSGrad2(1, i) = (cojac(1, 1)*BSGrad(1, i)+cojac(1, 2)*BSGrad(2, i))
                BSGrad2(2, i) = (cojac(2, 1)*BSGrad(1, i)+cojac(2, 2)*BSGrad(2, i))
            end do
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function feBSCEval(this, point) result(basisScalEval)
!
        implicit none
!
        class(FE_Basis), intent(in) :: this
        real(kind=8), dimension(3), intent(in) :: point
        real(kind=8), dimension(MAX_BS_CG) :: basisScalEval
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   evaluate fe basis scalar
!   In this                 : FE_Basis_scalar_cell
!   In point                : point where evaluate
!   Out basisScalEval       : evaluation of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        basisScalEval = 0.d0
        if (this%typeEF == EF_LAGRANGE) then
            call elrfvf(this%typema, point, basisScalEval)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!
!===================================================================================================
!
!===================================================================================================
!
    function feBSCGradEv(this, point, jacob_) result(BSGradEval)
!
        implicit none
!
        class(FE_Basis), intent(in) :: this
        real(kind=8), dimension(3), intent(in) :: point
        real(kind=8), dimension(3, 3), optional, intent(in) :: jacob_
        real(kind=8), dimension(3, MAX_BS_CG) :: BSGradEval
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   evaluate fe basis scalar
!   In this                 : FE_Basis_scalar_cell
!   In point                : point where evaluate
!   Out BSGradEval   : evaluation of the gradient of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        if (this%typeEF == EF_LAGRANGE) then
            if (present(jacob_)) then
                BSGradEval = FE_grad_lagr(this, point, jacob_)
            else
                BSGradEval = FE_grad_lagr(this, point)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
end module
