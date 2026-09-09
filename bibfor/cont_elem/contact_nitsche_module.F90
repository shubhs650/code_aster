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
module contact_nitsche_module
!
    use contact_type
    use contact_algebra_module
    use contact_module
!
    implicit none
!
    private
!
#include "asterf_types.h"
#include "asterfort/jevech.h"
#include "asterfort/reereg.h"
#include "asterfort/trace_mat.h"
#include "blas/dgemv.h"
#include "contact_module.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Nitsche methods
!
! --------------------------------------------------------------------------------------------------
!
!
    type ContactNitsche
! Young modulus
        real(kind=8) :: E = 0.d0
! Poisson ratio
        real(kind=8) :: nu = 0.d0
!
! stress at nodes
        real(kind=8), dimension(3, 3, 9) :: stress_nodes = 0.d0
!
    end type
!
    public :: ContactNitsche
    public :: dofsMapping, remappingVect, remappingMatr, nbDofsNitsche
    public :: evalStressNodes, evalStress, getMaterialProperties
    public :: dStress_n_du, dStress_nn_du
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    function dofsMapping(geom)
!
        implicit none
!
        integer(kind=8) :: dofsMapping(54)
        type(ContactGeom), intent(in) :: geom
!
! --------------------------------------------------------------------------------------------------
!
!   dofs mapping between face to volume
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node, i_elem, nb_dofs_slav, nb_dofs_volu
!
        dofsMapping = 0
!
        do i_node = 1, geom%nb_node_slav
            do i_elem = 1, geom%elem_dime
                dofsMapping((i_node-1)*geom%elem_dime+i_elem) = (geom%mapVolu2Surf(i_node)-1 &
                                                                 )*geom%elem_dime+i_elem
            end do
        end do
!
        nb_dofs_slav = geom%nb_node_slav*geom%elem_dime
        nb_dofs_volu = geom%nb_node_volu*geom%elem_dime
        do i_node = 1, geom%nb_node_mast
            do i_elem = 1, geom%elem_dime
                dofsMapping(nb_dofs_slav+(i_node-1)*geom%elem_dime+i_elem) = nb_dofs_volu+( &
                                                                             i_node-1 &
                                                                             )*geom%elem_dime+i_e&
                                                                             &lem
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine remappingVect(geom, dofsmap, vect, new_vect, coeff)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        integer(kind=8), intent(in) :: dofsMap(54)
        real(kind=8), intent(in) :: vect(MAX_LAGA_DOFS), coeff
        real(kind=8), intent(inout) :: new_vect(MAX_NITS_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Remapping vector (+= operation)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_dof, nb_dofs
!
        nb_dofs = (geom%nb_node_slav+geom%nb_node_mast)*geom%elem_dime
!
        do i_dof = 1, nb_dofs
            new_vect(dofsMap(i_dof)) = new_vect(dofsMap(i_dof))+coeff*vect(i_dof)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine remappingMatr(geom, dofsmap, matr, new_matr, coeff)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        integer(kind=8), intent(in) :: dofsMap(54)
        real(kind=8), intent(in) :: matr(MAX_LAGA_DOFS, MAX_LAGA_DOFS), coeff
        real(kind=8), intent(inout) :: new_matr(MAX_NITS_DOFS, MAX_NITS_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Remapping matrix (+= operation)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_dof, nb_dofs, j_dof, i_glo, j_glo
!
        nb_dofs = (geom%nb_node_slav+geom%nb_node_mast)*geom%elem_dime
!
        do j_dof = 1, nb_dofs
            j_glo = dofsMap(j_dof)
            do i_dof = 1, nb_dofs
                i_glo = dofsMap(i_dof)
                new_matr(i_glo, j_glo) = new_matr(i_glo, j_glo)+coeff*matr(i_dof, j_dof)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine gradDisp(geom, dshape_func_vo, grad, eps)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dshape_func_vo(3, MT_NNOMAX3D)
        real(kind=8), intent(out) :: grad(3, 3), eps(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of raytracing operator
!   (De^m(u)[v], Dgap(u)[v]) = (tau^m, -n^s )^{-1}(v^s - v^m  + g*Dn^s)
!
! --------------------------------------------------------------------------------------------------
!
!
        integer(kind=8) :: n, i, j
        real(kind=8) :: dfunc_dx(3, MT_NNOMAX3D)
!
        call gradFuncDispVolu(geom, dshape_func_vo, dfunc_dx)
!
        grad = 0.d0
!
        do n = 1, geom%nb_node_volu
            do i = 1, geom%elem_dime
                do j = 1, geom%elem_dime
                    grad(i, j) = grad(i, j)+dfunc_dx(j, n)*geom%depl_volu_curr(i, n)
                end do
            end do
        end do
!
        eps = (grad+transpose(grad))/2.d0
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine lameCoeff(E, nu, mu, lambda)
!
        implicit none
!
        real(kind=8), intent(in) :: E, nu
        real(kind=8), intent(out) :: mu, lambda
!
! --------------------------------------------------------------------------------------------------
!
!   Lame coefficient
!
! --------------------------------------------------------------------------------------------------
!
        lambda = E*nu/(1.d0+nu)/(1.d0-2.d0*nu)
        mu = 0.5d0*E/(1.d0+nu)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine evalStressNodes(geom, nits)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        type(ContactNitsche), intent(inout) :: nits
!
! --------------------------------------------------------------------------------------------------
!
!   Remapping matrix (+= operation)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: coor_qp_vo(3), mu, lambda
        real(kind=8) :: dshape_func_vo(3, MT_NNOMAX3D), grad(3, 3), eps(3, 3)
        integer(kind=8) :: i_node, iret
!
        call lameCoeff(nits%E, nits%nu, mu, lambda)
!
        do i_node = 1, geom%nb_node_slav
!
! ----- Projection of node on volumic slave cell (volumic parametric space)
!
            coor_qp_vo = 0.d0
            call reereg('S', geom%elem_volu_code, geom%nb_node_volu, geom%coor_volu_curr, &
                        geom%coor_slav_curr(1:3, i_node), geom%elem_dime, coor_qp_vo, iret, &
                        ndim_coor_=3)
!
! ----- Eval shape function and gradient
!
            call shapeFuncDispVolu(geom%elem_volu_code, coor_qp_vo, dshape_=dshape_func_vo)
            call gradDisp(geom, dshape_func_vo, grad, eps)
!
! ----- Eval stress
!
            nits%stress_nodes(1:3, 1:3, i_node) = 2.d0*mu*eps+lambda*trace_mat(3, eps)*Iden3()
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function evalStress(nits, nb_node_slav, shape_func_sl)
!
        implicit none
!
        type(ContactNitsche), intent(in) :: nits
        integer(kind=8), intent(in) :: nb_node_slav
        real(kind=8), intent(in) :: shape_func_sl(9)
        real(kind=8) :: evalStress(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Eval stress tensor
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, j
!
        evalStress = 0.d0
!
        do j = 1, 3
            do i = 1, 3
                evalStress(i, j) = evalPoly( &
                                   nb_node_slav, shape_func_sl, nits%stress_nodes(i, j, :))
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getMaterialProperties(nits)
!
        implicit none
!
        type(ContactNitsche), intent(inout) :: nits
!
! --------------------------------------------------------------------------------------------------
!
!   Read material properties
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: jcont
!
        call jevech('PCONFR', 'L', jcont)
!
        nits%E = zr(jcont+45)
        nits%nu = zr(jcont+46)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine nbDofsNitsche(geom, total_dofs, face_dofs, slav_dofs)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        integer(kind=8), intent(out) :: total_dofs, face_dofs, slav_dofs
!
! --------------------------------------------------------------------------------------------------
!
!   Compute number of dofs for Nitsche
!
! --------------------------------------------------------------------------------------------------
!
        total_dofs = geom%nb_dofs
        face_dofs = (geom%nb_node_slav+geom%nb_node_mast)*geom%elem_dime
        slav_dofs = geom%nb_node_volu*geom%elem_dime
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function dStress_n_du(geom, nits, norm_slav_init)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        type(ContactNitsche), intent(in) :: nits
        real(kind=8), intent(in) :: norm_slav_init(3)
        real(kind=8) :: dStress_n_du(MAX_NITS_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of normal stress
!   D stress_n^s(u^s)[v^s] = Aep(u^s) : grad(v^s)  * N^s
!
! --------------------------------------------------------------------------------------------------
!
        dStress_n_du = 0.d0
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dStress_nn_du(geom, stress_n, dStress_n, norm_slav, dNs)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: norm_slav(3), stress_n(3)
        real(kind=8), intent(in) :: dStress_n(MAX_NITS_DOFS, 3), dNs(MAX_LAGA_DOFS, 3)
        real(kind=8) :: dStress_nn_du(MAX_NITS_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of normal, normal stress
!   D stress_nn^s(u^s)[v^s] = Aep(u^s) : grad(v^s)  * N^s . n^s + stress(u^s) * N^s * Dn^s(u^s)[v^s]
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: dNs_sn(MAX_LAGA_DOFS)
        integer(kind=8) :: total_dofs, face_dofs, slav_dofs, i_dof
        integer(kind=8) :: dofsMap(54), slav_face_dofs
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
!
        dStress_nn_du = 0.d0
!
        call nbDofsNitsche(geom, total_dofs, face_dofs, slav_dofs)
        slav_face_dofs = geom%nb_node_slav*geom%elem_dime
!
! --- Term Aep(u^s) : grad(v^s)  * N^s . n^s
!
        do i_dof = 1, slav_dofs
            dStress_nn_du(i_dof) = dot_product(dStress_n(i_dof, 1:3), norm_slav)
        end do
!
! --- Term: stress(u^s) * N^s * D n^s(u^s)[v^s]
!
        dNs_sn = 0.d0
        b_lda = to_blas_int(MAX_LAGA_DOFS)
        b_m = to_blas_int(slav_face_dofs)
        b_n = to_blas_int(geom%elem_dime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, 1.d0, dNs, &
                   b_lda, stress_n, b_incx, 1.d0, dNs_sn, &
                   b_incy)
!
! --- Remapping
!
        dofsMap = dofsMapping(geom)
! do i_dof = 1, slav_face_dofs
!     dStress_nn_du(dofsMap(i_dof)) = dStress_nn_du(dofsMap(i_dof)) + dNs_sn(i_dof)
! end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
end module
