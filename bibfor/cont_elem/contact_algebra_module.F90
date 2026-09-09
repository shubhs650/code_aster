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
module contact_algebra_module
!
    use contact_type
    use contact_module
!
    implicit none
!
    private
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/mmtang.h"
#include "asterfort/matinv.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
#include "blas/dger.h"
#include "MeshTypes_type.h"
#include "contact_module.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Linear algebra and basic derivatives
!
! --------------------------------------------------------------------------------------------------
!
    public :: Heaviside, jump_norm, jump_tang, jump
    public :: dGap_du, d2Gap_du2, dZetaM_du, dPi_du
    public :: dTs_du_ns, dTs_du_nm, dTm_du_nm, dTs_du, dTm_du
    public :: dNs_du, d2Ns_du2_ns, d2Ns_du2_nm
    public :: projBs, projRm, projTn, otimes, Iden3, dvT_du
    public :: dProjBs_dx, dprojBs_ds, dprojBs_dn, dv_du
    public :: metricTensor, invMetricTensor, cmpTang, secondFundForm
    public :: dprojCoulomb_du, dprojCoulomb_dlc, dprojCoulomb_dlf
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function projRm(x)
!
        implicit none
!
        real(kind=8), intent(in) :: x
!
! --------------------------------------------------------------------------------------------------
!
!   Project x on R^- = min(x, 0)
!
! --------------------------------------------------------------------------------------------------
!
        if (x <= TOLE_BORNE) then
            projRm = x
        else
            projRm = 0.d0
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function Heaviside(x)
!
        implicit none
!
        real(kind=8), intent(in) :: x
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate Heaviside function
!
! --------------------------------------------------------------------------------------------------
!
        if (x >= -(TOLE_BORNE)) then
            Heaviside = 1.d0
        else
            Heaviside = 0.d0
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function otimes(v1, v2)
!
        implicit none
!
        real(kind=8), intent(in) :: v1(3), v2(3)
        real(kind=8) :: otimes(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate otimes product: (v1 \otimes v2 )_{i,j} = v1_{i} * v2_{j}
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, j
!
        do j = 1, 3
            do i = 1, 3
                otimes(i, j) = v1(i)*v2(j)
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function Iden3()
!
        implicit none
!
        real(kind=8) :: Iden3(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate identity matrix
!
! --------------------------------------------------------------------------------------------------
!
        Iden3 = 0.d0
        Iden3(1, 1) = 1.d0
        Iden3(2, 2) = 1.d0
        Iden3(3, 3) = 1.d0
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function projTn(normal)
!
        implicit none
!
        real(kind=8), intent(in) :: normal(3)
        real(kind=8) :: projTn(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Projection operator onto the tangent plane corresponding to normal vector n
!   Tn = Id - n \otimes n
!
! --------------------------------------------------------------------------------------------------
!
        projTn = Iden3()-otimes(normal, normal)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dProjTn_du(geom, tau_slav, dTs_ns, norm_slav)
!
        implicit none
!
        real(kind=8), intent(in) :: norm_slav(3)
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: tau_slav(3, 2)
        real(kind=8), intent(in) :: dTs_ns(MAX_LAGA_DOFS, 2)
        real(kind=8) :: dProjTn_du(MAX_LAGA_DOFS, 3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   First derivative of the projection operator onto the tangent plane
!   dTn_du =  -dNs_du \otimes n - n \otimes dNs_du
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: dns_du_(MAX_LAGA_DOFS, 3)
        integer(kind=8) :: i_dof
!
        dns_du_ = dNs_du(geom, tau_slav, dTs_ns)
!
        do i_dof = 1, MAX_LAGA_DOFS
            dProjTn_du(i_dof, 1:3, 1:3) = -otimes(dns_du_(i_dof, 1:3), norm_slav) &
                                          -otimes(norm_slav, dns_du_(i_dof, 1:3))
        end do
!
    end function
!
!===================================================================================================
!===================================================================================================
!
    function metricTensor(tau)
!
        implicit none
!
        real(kind=8), intent(in) :: tau(3, 2)
        real(kind=8) :: metricTensor(2, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Metric tensor: m[i,j] = tau_i . tau_j
!
! --------------------------------------------------------------------------------------------------
!
        metricTensor = matmul(transpose(tau), tau)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function invMetricTensor(geom, metricTens)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: metricTens(2, 2)
        real(kind=8) :: invMetricTensor(2, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Metric tensor: invm[i,j] = m^-1[i,j]
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: det
!
        invMetricTensor = 0.d0
!
        if (geom%elem_dime == 3) then
            det = metricTens(1, 1)*metricTens(2, 2)-metricTens(1, 2)*metricTens(2, 1)
            invMetricTensor(1, 1) = metricTens(2, 2)/det
            invMetricTensor(1, 2) = -metricTens(1, 2)/det
            invMetricTensor(2, 1) = -metricTens(2, 1)/det
            invMetricTensor(2, 2) = metricTens(1, 1)/det
        else
            invMetricTensor(1, 1) = 1.d0/metricTens(1, 1)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function secondFundForm(elem_dime, nb_node, elem_coor, ddff, norm)
!
        implicit none
!
        integer(kind=8), intent(in) :: elem_dime, nb_node
        real(kind=8), intent(in) :: elem_coor(3, MT_NNOMAX2D), norm(3), ddff(3, MT_NNOMAX2D)
        real(kind=8) :: secondFundForm(2, 2)
!
! --------------------------------------------------------------------------------------------------
!   Eval second fondamental form
!   H[i,j] = d^2 x / de_i de_j * norm
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: idim, ino
        real(kind=8) :: ddgeo1(3), ddgeo2(3), ddgeo3(3)
!
        secondFundForm = 0.d0
        ddgeo1 = 0.d0
        ddgeo2 = 0.d0
        ddgeo3 = 0.d0
!
! - CALCUL DE DDGEOMM
!
        do idim = 1, elem_dime
            do ino = 1, nb_node
                ddgeo1(idim) = ddgeo1(idim)+ddff(1, ino)*elem_coor(idim, ino)
                ddgeo2(idim) = ddgeo2(idim)+ddff(2, ino)*elem_coor(idim, ino)
                ddgeo3(idim) = ddgeo3(idim)+ddff(3, ino)*elem_coor(idim, ino)
            end do
        end do
!
! --- Evaluate H
!
        secondFundForm(1, 1) = ddgeo1(1)*norm(1)+ddgeo1(2)*norm(2)+ddgeo1(3)*norm(3)
        secondFundForm(1, 2) = ddgeo3(1)*norm(1)+ddgeo3(2)*norm(2)+ddgeo3(3)*norm(3)
        secondFundForm(2, 1) = secondFundForm(1, 2)
        secondFundForm(2, 2) = ddgeo2(1)*norm(1)+ddgeo2(2)*norm(2)+ddgeo2(3)*norm(3)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function cmpTang(invMetricTens, tau, vec)
!
        implicit none
!
        real(kind=8), intent(in) :: invMetricTens(2, 2)
        real(kind=8), intent(in) :: tau(3, 2), vec(3)
        real(kind=8) :: cmpTang(2)
!
! --------------------------------------------------------------------------------------------------
!
!   Compute coefficient in tangential basis (tau_1, tau_2)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: rhs(2)
!
        rhs = matmul(transpose(tau), vec)
!
        cmpTang = matmul(invMetricTens, rhs)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function projBs(param, x, s, norm_slav)
!
        implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3), norm_slav(3)
        real(kind=8), intent(in) :: s
        real(kind=8) :: projBs(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Project x on the sphere center to zero of radius s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3, 3), xT(3)
!
        projBs = 0.d0
        if (param%type_fric .ne. FRIC_TYPE_NONE) then
            if (abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
!
                if (norm_xT <= s) then
                    projBs = xT
                else
                    projBs = s*xT/norm_xT
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojBs_dx(param, x, s, norm_slav)
!
        implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3)
        real(kind=8), intent(in) :: s, norm_slav(3)
        real(kind=8) :: dprojBs_dx(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative along x of the projection on the sphere center to zero of radius s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3, 3), xT_uni(3), xT(3)
!
        dprojBs_dx = 0.d0
        if (param%type_fric .ne. FRIC_TYPE_NONE) then
            if (abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
!
                if (norm_xT <= s) then
                    dprojBs_dx = Tn
                else
                    xT_uni = xT/norm_xT
                    dprojBs_dx = s*norm_xT*(Tn-otimes(xT_uni, xT_uni))
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojBs_ds(param, x, s, norm_slav)
!
        implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3), norm_slav(3)
        real(kind=8), intent(in) :: s
        real(kind=8) :: dprojBs_ds(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative along s of the projection on the sphere center to zero of radius s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3, 3), xT(3)
!
        dprojBs_ds = 0
        if (param%type_fric .ne. FRIC_TYPE_NONE) then
            if (abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
!
                if (norm_xT > s) then
                    dprojBs_ds = xT/norm_xT
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojBs_dn(param, x, s, norm_slav)
!
        implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3)
        real(kind=8), intent(in) :: s, norm_slav(3)
        real(kind=8) :: dprojBs_dn(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative along the normal of the projection on the sphere center to zero of radius s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3, 3), xT_uni(3), x_n, xT(3)
!
        dprojBs_dn = 0
        if (param%type_fric .ne. FRIC_TYPE_NONE) then
            if (abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
                x_n = dot_product(x, norm_slav)
!
                if (norm_xT <= s) then
                    dprojBs_dn = -x_n*Tn-otimes(norm_slav, xT)
                else
                    xT_uni = xT/norm_xT
                    dprojBs_dn = -s/norm_xT*( &
                                 x_n*(Tn-otimes(xT_uni, xT_uni))+otimes(norm_slav, xT))
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojCoulomb_du(geom, param, gamma_f, x, s, ds_du, dv_du_, norm_slav, &
                             tau_slav, dts_ns)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3), ds_du(MAX_LAGA_DOFS), gamma_f, tau_slav(3, 2)
        real(kind=8), intent(in) :: s, norm_slav(3), dv_du_(MAX_LAGA_DOFS, 3)
        real(kind=8), intent(in) ::  dts_ns(MAX_LAGA_DOFS, 2)
        real(kind=8) :: dprojCoulomb_du(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative of the projection of the friction augmented Lagrange multiplier
!   on the Coulomb cone with respect to the displacement
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_dof
        real(kind=8) :: norm_xT, Tn(3, 3), xT_uni(3), xT(3), coeff, TnoT(3, 3)
        real(kind=8) :: dTn(MAX_LAGA_DOFS, 3, 3)
        blas_int :: b_MAX_LAGA_DOFS, b_3, b_1
!
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_3 = 3
        b_1 = 1
!
        dprojCoulomb_du = 0
        if (param%type_fric .ne. FRIC_TYPE_NONE) then
            if (abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                dTn = dProjTn_du(geom, tau_slav, dts_ns, norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
                !
                if (norm_xT <= s) then
                    ! -gamma_f * Tn * dv_du + DTn[du] x
                    do i_dof = 1, MAX_LAGA_DOFS
                        dprojCoulomb_du(i_dof, 1:3) = -gamma_f*matmul(Tn, dv_du_(i_dof, 1:3))+ &
                                                      matmul(dTn(i_dof, :, :), x)
                    end do
                else
                    ! ds_du * (xT/|xT|) + s/|xT| * (Id - xT/|xT| \otimes xT/|xT|)*
                    !                                                    (-gamma Tn dv_du + dTn x)
                    xT_uni = xT/norm_xT
                    !     call dger(b_MAX_LAGA_DOFS, b_3, 1.d0, ds_du, b_1, &
                    !               xT_uni, b_1, dprojCoulomb_du, b_MAX_LAGA_DOFS)
                    do i_dof = 1, MAX_LAGA_DOFS
                        dprojCoulomb_du(i_dof, 1:3) = ds_du(i_dof)*xT_uni(1:3)
                    end do
                    TnoT = Iden3()-otimes(xT_uni, xT_uni)
                    coeff = s/norm_xT
                    do i_dof = 1, b_MAX_LAGA_DOFS
                        dprojCoulomb_du(i_dof, 1:3) = dprojCoulomb_du(i_dof, 1:3)+ &
                                                      coeff*(matmul(TnoT, -gamma_f* &
                                                                    matmul(Tn, dv_du_(i_dof, 1:3)) &
                                                                    +matmul(dTn(i_dof, :, :), x)))
                    end do
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojCoulomb_dlc(param, x, s, ds_dl, norm_slav)
!
        implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3), ds_dl(MAX_LAGA_DOFS)
        real(kind=8), intent(in) :: s, norm_slav(3)
        real(kind=8) :: dprojCoulomb_dlc(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative of the projection of the friction augmented Lagrange multiplier
!   on the Coulomb cone with respect to the contact Lagrange multiplier
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3, 3), xT_uni(3), xT(3)
        blas_int :: b_MAX_LAGA_DOFS, b_3, b_1
!
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_3 = 3
        b_1 = 1
!
        dprojCoulomb_dlc = 0
        if (param%type_fric .ne. FRIC_TYPE_NONE) then
            if (abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
                !
                if (norm_xT <= s) then
                    dprojCoulomb_dlc = 0
                else
                    ! ds_dl * x/|x|
                    xT_uni = xT/norm_xT
                    call dger(b_MAX_LAGA_DOFS, b_3, 1.d0, ds_dl, b_1, &
                              xT_uni, b_1, dprojCoulomb_dlc, b_MAX_LAGA_DOFS)
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojCoulomb_dlf(param, x, s, mu_g, ds_dl, norm_slav, gamma_speed)
!
        implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3), mu_g(MAX_LAGA_DOFS, 3), gamma_speed(3)
        real(kind=8), intent(in) :: s, norm_slav(3), ds_dl(MAX_LAGA_DOFS)
        real(kind=8) :: dprojCoulomb_dlf(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative of the projection of the friction augmented Lagrange multiplier
!   on the Couomb cone with respect to the friction Lagrange multiplier
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_dof
        real(kind=8) :: norm_xT, Tn(3, 3), xT_uni(3), xT(3), coeff, TnoT(3, 3)
        blas_int :: b_MAX_LAGA_DOFS, b_3, b_1
!
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_3 = 3
        b_1 = 1
!
        dprojCoulomb_dlf = 0.d0
        if (param%type_fric .ne. FRIC_TYPE_NONE) then
            if (abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
                !
                if (norm_xT <= s) then
                    ! Tn * dmu_g
                    coeff = 1.d0
                    do i_dof = 1, MAX_LAGA_DOFS
                        dprojCoulomb_dlf(i_dof, 1:3) = matmul(Tn, mu_g(i_dof, 1:3))
                    end do
                else
                    !  ds_dlag*xT/|xT| + s/|xT| * (Id - xT/|xT| \otimes xT/|xT|)*Tn*
                    !                                                      (dmu_g-gamma_f*speed)
                    xT_uni = xT/norm_xT
                    coeff = s/norm_xT
                    TnoT = Iden3()-otimes(xT_uni, xT_uni)
                    do i_dof = 1, MAX_LAGA_DOFS
                        dprojCoulomb_dlf(i_dof, 1:3) = ds_dl(i_dof)*xT_uni+ &
                                                       coeff*matmul(TnoT, &
                                                                    matmul(Tn, mu_g(i_dof, :)- &
                                                                           gamma_speed))
                    end do
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dNs_du(geom, tau_slav, dTs_ns)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: tau_slav(3, 2)
        real(kind=8), intent(in) :: dTs_ns(MAX_LAGA_DOFS, 2)
        real(kind=8) :: dNs_du(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of slave normal
!   D n^s(u^s)[v^s] = -tau^s * invMetricTens^s * D tau^s(u^s)[v^s]  * n^s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: lhs_inv(3, 2), metricTens(2, 2), invMetricTens(2, 2)
        blas_int :: b_dime_m1, b_MAX_LAGA_DOFS, b_3, b_nb_dofs, b_dime
!
        b_dime_m1 = to_blas_int(geom%elem_dime-1)
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_nb_dofs = to_blas_int(geom%nb_dofs)
        b_dime = to_blas_int(geom%elem_dime)
        b_3 = to_blas_int(3)
!
        metricTens = metricTensor(tau_slav)
        invMetricTens = invMetricTensor(geom, metricTens)
!
! ----- Term: -tau_slav * invMetricTens_slav
!
        lhs_inv = -matmul(tau_slav, invMetricTens)
!
! --- Compute solution: Dns = lhs^-1 * dTs_ns
!
        dNs_du = 0.d0
        call dgemm("N", "T", b_nb_dofs, b_dime, b_dime_m1, &
                   1.d0, dTs_ns, b_MAX_LAGA_DOFS, lhs_inv, b_3, &
                   0.d0, dNs_du, b_MAX_LAGA_DOFS)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function d2Ns_du2_ns(geom, tau_slav, dTs_ns)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: tau_slav(3, 2)
        real(kind=8), intent(in) :: dTs_ns(MAX_LAGA_DOFS, 2)
        real(kind=8) :: d2Ns_du2_ns(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate second derivative of slave normal
!   D^2 n^s(u^s)[v^s, dv^s] . n^s =
!      - (D tau^s(u^s)[v^s]  * n^s)^T * invMetricTens^s * (D tau^s(u^s)[v^s]  * n^s)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: invA_dT_ns(MAX_LAGA_DOFS, 2), metricTens(2, 2), invMetricTens(2, 2)
        blas_int :: b_dime_m1, b_MAX_LAGA_DOFS, b_2, b_nb_dofs
!
        b_dime_m1 = to_blas_int(geom%elem_dime-1)
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_nb_dofs = to_blas_int(geom%nb_dofs)

!
        metricTens = metricTensor(tau_slav)
        invMetricTens = invMetricTensor(geom, metricTens)
!
! ----- Term: -invMetricTens_slav * (D tau^s(u^s)[v^s]  * n^s)
!
        invA_dT_ns = 0.d0
        b_2 = to_blas_int(2)
        call dgemm("N", "T", b_nb_dofs, b_dime_m1, b_dime_m1, &
                   -1.d0, dTs_ns, b_MAX_LAGA_DOFS, invMetricTens, b_2, &
                   0.d0, invA_dT_ns, b_MAX_LAGA_DOFS)
!
! --- Compute D2 ns . n^s
!
        d2Ns_du2_ns = 0.d0
        call dgemm("N", "T", b_nb_dofs, b_nb_dofs, b_dime_m1, &
                   1.d0, dTs_ns, b_MAX_LAGA_DOFS, invA_dT_ns, b_MAX_LAGA_DOFS, &
                   0.d0, d2Ns_du2_ns, b_MAX_LAGA_DOFS)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function d2Ns_du2_nm(geom, tau_slav, norm_mast, dNs, dTs, &
                         dTs_ns, dTs_nm)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: tau_slav(3, 2), norm_mast(3)
        real(kind=8), intent(in) :: dNs(MAX_LAGA_DOFS, 3), dTs(MAX_LAGA_DOFS, 3, 2)
        real(kind=8), intent(in) :: dTs_nm(MAX_LAGA_DOFS, 2), dTs_ns(MAX_LAGA_DOFS, 2)
        real(kind=8) :: d2Ns_du2_nm(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate second derivative of slave normal
!   D^2 {n^s}(u^s)[v^s, du^s] . n^m  =
!       -( D tau^s(u^s)[v^s](u^s)[v^s] n^s)^T * invMetricTens^s * ( D tau^s(u^s)[v^s] n^m)
!       + (tau^s)^T n^m ( D invMetricTens^s(u^s){\vdu} D tau^s(u^s)[v^s] n^s
!                        -  invMetricTens^s D tau^s(u^s)[v^s] D n^s(u^s)[du^s] )
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: invA_dT_ns(MAX_LAGA_DOFS, 2), metricTens(2, 2), invMetricTens(2, 2)
        real(kind=8) :: ts_nm(2), invA_ts_nm(2), dinvMetric(MAX_LAGA_DOFS, 3, 2)
        integer(kind=8) :: i_tau
        blas_int :: b_dime_m1, b_MAX_LAGA_DOFS, b_2, b_nb_dofs, b_dime
!
        b_dime_m1 = to_blas_int(geom%elem_dime-1)
        b_dime = to_blas_int(geom%elem_dime)
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_nb_dofs = to_blas_int(geom%nb_dofs)
!
        d2Ns_du2_nm = 0.d0
        metricTens = metricTensor(tau_slav)
        invMetricTens = invMetricTensor(geom, metricTens)
        ts_nm = matmul(transpose(tau_slav), norm_mast)
        invA_ts_nm = matmul(invMetricTens, ts_nm)
!
! --- Term: D invMetricTens^s(u^s){\vdu} -> not implemented
!
        dinvMetric = 0.d0
!
! --- Term: -invMetricTens_slav * (D tau^s(u^s)[v^s]  * n^s)
!
        invA_dT_ns = 0.d0
        b_2 = to_blas_int(2)
        call dgemm("N", "T", b_nb_dofs, b_dime_m1, b_dime_m1, &
                   -1.d0, dTs_ns, b_MAX_LAGA_DOFS, invMetricTens, b_2, &
                   0.d0, invA_dT_ns, b_MAX_LAGA_DOFS)
!
! --- Term: -( D tau^s(u^s)[v^s](u^s)[v^s] n^s)^T * invMetricTens^s * ( D tau^s(u^s)[v^s] n^m)
!
        call dgemm("N", "T", b_nb_dofs, b_nb_dofs, b_dime_m1, &
                   1.d0, dTs_nm, b_MAX_LAGA_DOFS, invA_dT_ns, b_MAX_LAGA_DOFS, &
                   0.d0, d2Ns_du2_nm, b_MAX_LAGA_DOFS)
!
! --- Term: -((tau^s)^T n^m  invMetricTens^s ) . D tau^s(u^s)[v^s] D n^s(u^s)[du^s]
!
        do i_tau = 1, geom%elem_dime-1
            call dgemm("N", "T", b_nb_dofs, b_nb_dofs, b_dime, &
                       -invA_ts_nm(i_tau), dTs(:, :, i_tau), b_MAX_LAGA_DOFS, dNs, &
                       b_MAX_LAGA_DOFS, 1.d0, d2Ns_du2_nm, b_MAX_LAGA_DOFS)
        end do
!
! --- Term: (tau^s)^T n^m ( D invMetricTens^s(u^s){\vdu} D tau^s(u^s)[v^s] n^s )
!
        do i_tau = 1, geom%elem_dime-1
            call dgemm("N", "T", b_nb_dofs, b_nb_dofs, b_dime_m1, &
                       -invA_ts_nm(i_tau), dinvMetric(:, :, i_tau), b_MAX_LAGA_DOFS, dTs_ns, &
                       b_MAX_LAGA_DOFS, 1.d0, d2Ns_du2_nm, b_MAX_LAGA_DOFS)
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dT_du(elem_dime, nb_node, dfunc)
!
        implicit none
!
        integer(kind=8), intent(in) :: elem_dime, nb_node
        real(kind=8), intent(in) :: dfunc(2, MT_NNOMAX2D)
        real(kind=8) :: dT_du(MT_NNOMAX3D, 3, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of  tangent
!
!   D t^a = dfunc_dzeta^a * e_j
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node, i_dim, index, i_tau
!
        dT_du = 0.d0
!
        index = 0
        do i_node = 1, nb_node
            do i_dim = 1, elem_dime
                index = index+1
                do i_tau = 1, elem_dime-1
                    dT_du(index, i_dim, i_tau) = dfunc(i_tau, i_node)
                end do
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dTs_du(geom, dfunc_slav)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dfunc_slav(2, MT_NNOMAX2D)
        real(kind=8) :: dTs_du(MAX_LAGA_DOFS, 3, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of slave tangent
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: dT(MT_NNOMAX3D, 3, 2)
        integer(kind=8) :: i_node, i_dim, index, index2
!
        dTs_du = 0.d0
!
        dT = dT_du(geom%elem_dime, geom%nb_node_slav, dfunc_slav)
!
! --- Slave side
!
        index = 0
        index2 = 0
        do i_node = 1, geom%nb_node_slav
            do i_dim = 1, geom%elem_dime
                index = index+1
                index2 = index2+1
                dTs_du(index, 1:geom%elem_dime, 1:geom%elem_dime-1) = dT( &
                                                                      index2, 1:geom%elem_dime, &
                                                                      1:geom%elem_dime-1 &
                                                                      )
            end do
!
            index = index+geom%indi_lagc(i_node)
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dTm_du(geom, dfunc_mast)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dfunc_mast(2, MT_NNOMAX2D)
        real(kind=8) :: dTm_du(MAX_LAGA_DOFS, 3, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of master tangent
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: dT(MT_NNOMAX3D, 3, 2)
        integer(kind=8) :: i_node, i_dim, index, index2
!
        dTm_du = 0.d0
!
        dT = dT_du(geom%elem_dime, geom%nb_node_mast, dfunc_mast)
!
! --- Master side
!
        index = geom%nb_dofs-geom%nb_node_mast*geom%elem_dime
        index2 = 0
        do i_node = 1, geom%nb_node_mast
            do i_dim = 1, geom%elem_dime
                index = index+1
                index2 = index2+1
                dTm_du(index, 1:geom%elem_dime, 1:geom%elem_dime-1) = dT( &
                                                                      index2, 1:geom%elem_dime, &
                                                                      1:geom%elem_dime-1 &
                                                                      )
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dT_du_norm(elem_dime, nb_node, dfunc, norm)
!
        implicit none
!
        integer(kind=8), intent(in) :: elem_dime, nb_node
        real(kind=8), intent(in) :: dfunc(2, MT_NNOMAX2D), norm(3)
        real(kind=8) :: dT_du_norm(MT_NNOMAX3D, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of tangent dot with normal
!
!   D t^a . n = dfunc_dzeta^a * e_j . norm
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node, i_dim, index, i_tau
!
        dT_du_norm = 0.d0
!
        index = 0
        do i_node = 1, nb_node
            do i_dim = 1, elem_dime
                index = index+1
                do i_tau = 1, elem_dime-1
                    dT_du_norm(index, i_tau) = dfunc(i_tau, i_node)*norm(i_dim)
                end do
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dTs_du_ns(geom, dfunc_slav, norm_slav)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dfunc_slav(2, MT_NNOMAX2D), norm_slav(3)
        real(kind=8) :: dTs_du_ns(MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of slave tangent
!   D tau^s (u)[v] . n^s
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: dT_n(MT_NNOMAX3D, 2)
        integer(kind=8) :: i_node, i_dim, index, index2
!
        dTs_du_ns = 0.d0
!
        dT_n = dT_du_norm(geom%elem_dime, geom%nb_node_slav, dfunc_slav, norm_slav)
!
! --- Slave side
!
        index = 0
        index2 = 0
        do i_node = 1, geom%nb_node_slav
            do i_dim = 1, geom%elem_dime
                index = index+1
                index2 = index2+1
                dTs_du_ns(index, 1:geom%elem_dime-1) = dT_n(index2, 1:geom%elem_dime-1)
            end do
!
            index = index+geom%indi_lagc(i_node)
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dTs_du_nm(geom, dfunc_slav, norm_mast)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dfunc_slav(2, MT_NNOMAX2D), norm_mast(3)
        real(kind=8) :: dTs_du_nm(MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of slave tangent
!   D tau^s (u)[v] . n^m
! --------------------------------------------------------------------------------------------------
!
        dTs_du_nm = dTs_du_ns(geom, dfunc_slav, norm_mast)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dTm_du_nm(geom, dfunc_mast, norm_mast)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dfunc_mast(2, MT_NNOMAX2D), norm_mast(3)
        real(kind=8) :: dTm_du_nm(MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of master tangent
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: dT_n(MT_NNOMAX3D, 2)
        integer(kind=8) :: i_node, i_dim, index, index2
!
        dTm_du_nm = 0.d0
!
        dT_n = dT_du_norm(geom%elem_dime, geom%nb_node_mast, dfunc_mast, norm_mast)
!
! --- Master side
!
        index = geom%nb_dofs-geom%nb_node_mast*geom%elem_dime
        index2 = 0
        do i_node = 1, geom%nb_node_mast
            do i_dim = 1, geom%elem_dime
                index = index+1
                index2 = index2+1
                dTm_du_nm(index, 1:geom%elem_dime-1) = dT_n(index2, 1:geom%elem_dime-1)
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dPi_du(geom, norm_slav, tau_mast, gap, jump_v, &
                      dNs, dGap, dZetaM)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: norm_slav(3), tau_mast(3, 2), gap
        real(kind=8), intent(in) :: jump_v(MAX_LAGA_DOFS, 3), dNs(MAX_LAGA_DOFS, 3)
        real(kind=8), intent(out), optional :: dGap(MAX_LAGA_DOFS), dZetaM(MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of raytracing operator
!   (De^m(u)[v], Dgap(u)[v]) = (tau^m, -n^s )^{-1}(v^s - v^m  + g*Dn^s)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: det, rhs(MAX_LAGA_DOFS, 3), lhs(3, 3), lhs_inv(3, 3)
        real(kind=8) :: sol(3, MAX_LAGA_DOFS)
        blas_int :: b_MAX_LAGA_DOFS, b_3, b_nb_dofs
!
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_3 = to_blas_int(3)
        b_nb_dofs = to_blas_int(geom%nb_dofs)
!
        dZetaM = 0.d0
        dGap = 0.d0
        sol = 0.d0
!
! --- lhs : (tau^m, -n^s )
!
        if (geom%elem_dime == 2) then
            lhs(1:3, 1) = tau_mast(1:3, 1)
            lhs(1:3, 2) = -norm_slav
            lhs(1:3, 3) = [0.d0, 0.d0, 1.d0]
        else
            lhs(1:3, 1:2) = tau_mast
            lhs(1:3, 3) = -norm_slav
        end if
!
        call matinv("S", 3, lhs, lhs_inv, det)
!
! --- rhs : (v^s - v^m  + g*Dn^s)
!
        rhs = jump_v+gap*dNs
!
! --- Compute solution: (De^m(u)[v], Dgap(u)[v]) = lhs^-1 * rhs
!
        call dgemm("N", "T", b_3, b_nb_dofs, b_3, &
                   1.d0, lhs_inv, b_3, rhs, b_MAX_LAGA_DOFS, &
                   0.d0, sol, b_3)
!
! --- Extract terms
!
        if (present(dZetaM)) then
            dZetaM(1:geom%nb_dofs, 1:geom%elem_dime-1) = transpose( &
                                                         sol(1:geom%elem_dime-1, 1:geom%nb_dofs))
        end if
!
        if (present(dGap)) then
            dGap(1:geom%nb_dofs) = sol(geom%elem_dime, 1:geom%nb_dofs)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function dZetaM_du(geom, norm_slav, tau_mast, gap, jump_v, &
                       dNs)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: norm_slav(3), tau_mast(3, 2), gap
        real(kind=8), intent(in) :: jump_v(MAX_LAGA_DOFS, 3), dNs(MAX_LAGA_DOFS, 3)
        real(kind=8) :: dZetaM_du(MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   First derivative of projected master convertive coordinates
!   dZetaM_du = invMetricTensor^m * tau^m^T*(v^s - v^m + Dgn * n^s + gn*Dn^s)
!
! --------------------------------------------------------------------------------------------------
!
        call dPi_du(geom, norm_slav, tau_mast, gap, jump_v, &
                    dNs, dZetaM=dZetaM_du)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dGap_du(geom, norm_slav, tau_mast, gap, jump_v, &
                     dNs)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: norm_slav(3), tau_mast(3, 2), gap
        real(kind=8), intent(in) :: jump_v(MAX_LAGA_DOFS, 3), dNs(MAX_LAGA_DOFS, 3)
        real(kind=8) :: dGap_du(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of gap for raytracing
!   D gap(u)[v] = -(v^s - v^m + gap*Dn^s[v^s]).n^m/(n^m.n^s)
!
! --------------------------------------------------------------------------------------------------
!
        call dPi_du(geom, norm_slav, tau_mast, gap, jump_v, &
                    dNs, dGap=dGap_du)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dvT_du2(geom, speed, norm_slav, proj_tole, dGap_du_, &
                     invMetricTens, tau_slav, dZetaM_du_, dts_ns, coor_qp_sl, dshape_func_ma)
!
        implicit none
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: tau_slav(3, 2), norm_slav(3), speed(3), coor_qp_sl(2)
        real(kind=8), intent(in) :: dGap_du_(MAX_LAGA_DOFS), dts_ns(MAX_LAGA_DOFS, 2)
        real(kind=8), intent(in) :: proj_tole, invMetricTens(2, 2), dshape_func_ma(2, MT_NNOMAX2D)
        real(kind=8), intent(in) :: dZetaM_du_(MAX_LAGA_DOFS, 2)
        real(kind=8) :: dvT_du2(MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of tangential speed
!   D vT(u)[v] ~  A^{-1} *\tau_s * ((D Tn(u)[v]) v + Tn (D v(u)[v]))
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_dof
        real(kind=8) :: dv_du_(MAX_LAGA_DOFS, 3), Tn(3, 3), dTn(MAX_LAGA_DOFS, 3, 3)
        real(kind=8) :: dvT_g(MAX_LAGA_DOFS, 3)
!
! ------ First derivative of tangential speed
!
        dv_du_ = dv_du(geom, proj_tole, dGap_du_, dZetaM_du_, coor_qp_sl, dshape_func_ma)
!
! ------ Tangential components in global coordinates
!
        Tn = projTn(norm_slav)
        dTn = dProjTn_du(geom, tau_slav, dts_ns, norm_slav)
!
        do i_dof = 1, MAX_LAGA_DOFS
            dvT_g(i_dof, 1:3) = matmul(dTn(i_dof, :, :), speed)+matmul(Tn, dv_du_(i_dof, :))
        end do
!
! ------ Switch to local coordinates
!
        do i_dof = 1, MAX_LAGA_DOFS
            dvT_du2(i_dof, 1:2) = cmpTang(invMetricTens, tau_slav, dvT_g(i_dof, 1:3))
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
    function dvT_du(geom, speed, norm_slav, proj_tole, dGap_du_, &
                    invMetricTens, tau_slav, dZetaM_du_, dts_ns, coor_qp_sl, dshape_func_ma)
!
        implicit none
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: tau_slav(3, 2), norm_slav(3), speed(3), coor_qp_sl(2)
        real(kind=8), intent(in) :: dGap_du_(MAX_LAGA_DOFS), dts_ns(MAX_LAGA_DOFS, 2)
        real(kind=8), intent(in) :: proj_tole, invMetricTens(2, 2), dshape_func_ma(2, MT_NNOMAX2D)
        real(kind=8), intent(in) :: dZetaM_du_(MAX_LAGA_DOFS, 2)
        real(kind=8) :: dvT_du(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of tangential relative speed in global coordinates
!   D vT(u)[v] = ((D Tn(u)[v]) v + Tn (D v(u)[v]))
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_dof
        real(kind=8) :: dv_du_(MAX_LAGA_DOFS, 3), Tn(3, 3), dTn(MAX_LAGA_DOFS, 3, 3)
!
! ------ First derivative of tangential speed
!
        dv_du_ = dv_du(geom, proj_tole, dGap_du_, dZetaM_du_, coor_qp_sl, dshape_func_ma)
!
! ------ Tangential components in global coordinates
!
        Tn = projTn(norm_slav)
        dTn = dProjTn_du(geom, tau_slav, dts_ns, norm_slav)
!
        do i_dof = 1, MAX_LAGA_DOFS
            dvT_du(i_dof, 1:3) = matmul(dTn(i_dof, :, :), speed)+matmul(Tn, dv_du_(i_dof, :))
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dv_du(geom, proj_tole, dGap_du_, dZetaM_du_, coor_qp_sl, dshape_func_ma)
!
        implicit none
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dGap_du_(MAX_LAGA_DOFS), proj_tole, coor_qp_sl(2)
        real(kind=8), intent(in) :: dZetaM_du_(MAX_LAGA_DOFS, 2), dshape_func_ma(2, MT_NNOMAX2D)
        real(kind=8) :: dv_du(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of relative speed in global coordinates
!   D v(u)[v] = (tau^m_0 dZetaM_du - D g[v] n^s_0)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_slav_prev(3)
        real(kind=8) :: dGap_du_ns(MAX_LAGA_DOFS, 3)
        real(kind=8) :: tau_mast_prev(3, 2), tau_dZeta(3)
        integer(kind=8) :: i_dof
!
! ------ Compute outward previous slave normal
        call apnorm(geom%nb_node_slav, geom%elem_slav_code, geom%elem_dime, geom%coor_slav_prev, &
                    coor_qp_sl(1), coor_qp_sl(2), norm_slav_prev)
! ------ Compute outward previous tangent base
        tau_mast_prev = 0.d0
        call mmtang(geom%elem_dime, geom%nb_node_mast, geom%coor_mast_prev, dshape_func_ma, &
                    tau_mast_prev(:, 1), tau_mast_prev(:, 2))
!
        dGap_du_ns = 0.d0
        do i_dof = 1, MAX_LAGA_DOFS
            dGap_du_ns(i_dof, 1:3) = dGap_du_(i_dof)*norm_slav_prev
        end do
!
        dv_du = 0.d0
        do i_dof = 1, MAX_LAGA_DOFS
            tau_dZeta = matmul(tau_mast_prev, dZetaM_du_(i_dof, 1:2))
            dv_du(i_dof, 1:3) = (tau_dZeta-dGap_du_ns(i_dof, 1:3))
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function d2Gap_du2(geom, norm_slav, norm_mast, gap, Hm, dNs, dGap, dZetaM, dTm_nm, d2Ns_nm)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: norm_slav(3), norm_mast(3), Hm(2, 2), gap
        real(kind=8), intent(in) :: dNs(MAX_LAGA_DOFS, 3), dGap(MAX_LAGA_DOFS)
        real(kind=8), intent(in) :: d2Ns_nm(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
        real(kind=8), intent(in) :: dZetaM(MAX_LAGA_DOFS, 2), dTm_nm(MAX_LAGA_DOFS, 2)
        real(kind=8) :: d2Gap_du2(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate second derivative of gap
!   D^2 gap(u)[v,w] =
!   - ( D gap(u)[v]* D n^s[w^s] + D gap(u)[w]* D n^s[v^s] + gap(u) D^2 n^s[v^s, w^s] ).n^m/(n^m.n^s)
!   + ( D tau^m[v] .n^m * De[w] + D tau^m[w] .n^m * De[v] + H^m * De[v]*De[w] )/(n^m.n^s)
!   whit H^m = d tau^m / de .n^m
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm(3), dNs_n(MAX_LAGA_DOFS), dZetaM_H(MAX_LAGA_DOFS, 2), inv_ns_nm
        blas_int :: b_dime_m1, b_MAX_LAGA_DOFS, b_nb_dofs, b_dime, b_incx, b_incy, b_2
!
        b_dime_m1 = to_blas_int(geom%elem_dime-1)
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_nb_dofs = to_blas_int(geom%nb_dofs)
        b_dime = to_blas_int(geom%elem_dime)
!
        d2Gap_du2 = 0.d0
!
! --- Term: n^m/(n^m.n^s)
!
        inv_ns_nm = 1.d0/dot_product(norm_mast, norm_slav)
        norm = norm_mast*inv_ns_nm
!
! --- Term: D(n^s[v^s]).n^m/(n^m.n^s)
!
        dNs_n = 0.d0
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_nb_dofs, b_dime, 1.d0, dNs, &
                   b_MAX_LAGA_DOFS, norm, b_incx, 1.d0, dNs_n, &
                   b_incy)
!
! --- Term: -(D gap(u)[v]* D n^s[w^s] + D gap(u)[w]* D n^s[v^s]).n^m/(n^m.n^s)
!
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dger(b_nb_dofs, b_nb_dofs, -1.d0, dGap, b_incx, &
                  dNs_n, b_incy, d2Gap_du2, b_MAX_LAGA_DOFS)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dger(b_nb_dofs, b_nb_dofs, -1.d0, dNs_n, b_incx, &
                  dGap, b_incy, d2Gap_du2, b_MAX_LAGA_DOFS)
!
! --- Term: -gap(u) D^2 n^s[v^s, w^s].n^m/(n^m.n^s)
!
        d2Gap_du2(1:geom%nb_dofs, &
                  1:geom%nb_dofs) = d2Gap_du2(1:geom%nb_dofs, 1:geom%nb_dofs) &
                          -(gap*inv_ns_nm)*d2Ns_nm(1:geom%nb_dofs, 1:geom%nb_dofs)
!
! --- Term: ( D tau^m[v] .n^m * De[w] + D tau^m[w] .n^m * De[v] )/(n^m.n^s)
!
        call dgemm("N", "T", b_nb_dofs, b_nb_dofs, b_dime_m1, &
                   inv_ns_nm, dTm_nm, b_MAX_LAGA_DOFS, dZetaM, b_MAX_LAGA_DOFS, &
                   1.d0, d2Gap_du2, b_MAX_LAGA_DOFS)
        call dgemm("N", "T", b_nb_dofs, b_nb_dofs, b_dime_m1, &
                   inv_ns_nm, dZetaM, b_MAX_LAGA_DOFS, dTm_nm, b_MAX_LAGA_DOFS, &
                   1.d0, d2Gap_du2, b_MAX_LAGA_DOFS)
!
! --- Term: H^m * De[v]
!
        dZetaM_H = 0.d0
        b_2 = to_blas_int(2)
        call dgemm("N", "N", b_nb_dofs, b_dime_m1, b_dime_m1, &
                   1.d0, dZetaM, b_MAX_LAGA_DOFS, Hm, b_2, &
                   0.d0, dZetaM_H, b_MAX_LAGA_DOFS)
!
! --- Term: (H^m * De[v])*De[w] / (n^m.n^s)
!
        call dgemm("N", "T", b_nb_dofs, b_nb_dofs, b_dime_m1, &
                   inv_ns_nm, dZetaM, b_MAX_LAGA_DOFS, dZetaM_H, b_MAX_LAGA_DOFS, &
                   1.d0, d2Gap_du2, b_MAX_LAGA_DOFS)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function jump(geom, func_slav, func_mast)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_slav(9), func_mast(9)
        real(kind=8) :: jump(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate  jump of test function:
!   jump = (v^s - v^m)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node, i_dim, index
!
        jump = 0.d0
        index = 0
!
! --- Slave side
!
        do i_node = 1, geom%nb_node_slav
            do i_dim = 1, geom%elem_dime
                index = index+1
                jump(index, i_dim) = func_slav(i_node)
            end do
!
            index = index+geom%indi_lagc(i_node)
        end do
!
! --- Master side
!
        do i_node = 1, geom%nb_node_mast
            do i_dim = 1, geom%elem_dime
                index = index+1
                jump(index, i_dim) = -func_mast(i_node)
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function jump_tang(geom, func_slav, func_mast, norm_slav)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_slav(9), func_mast(9)
        real(kind=8), intent(in) :: norm_slav(3)
        real(kind=8) :: jump_tang(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate tangential jum of test function:
!   jump_tang = Tn * jump
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node, i_dim, index, j_dim
        real(kind=8) :: Tn(3, 3)
!
        jump_tang = 0.d0
        index = 0
        Tn = projTn(norm_slav)
!
! --- Slave side
!
        do i_node = 1, geom%nb_node_slav
            do i_dim = 1, geom%elem_dime
                index = index+1
                do j_dim = 1, geom%elem_dime
                    jump_tang(index, j_dim) = func_slav(i_node)*Tn(j_dim, i_dim)
                end do
            end do
!
            index = index+geom%indi_lagc(i_node)
        end do
!
! --- Master side
!
        do i_node = 1, geom%nb_node_mast
            do i_dim = 1, geom%elem_dime
                index = index+1
                do j_dim = 1, geom%elem_dime
                    jump_tang(index, j_dim) = -func_mast(i_node)*Tn(j_dim, i_dim)
                end do
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function jump_norm(geom, jump_v, norm_)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: jump_v(MAX_LAGA_DOFS, 3)
        real(kind=8), intent(in) :: norm_(3)
        real(kind=8) :: jump_norm(MAX_LAGA_DOFS)
        blas_int :: b_incx, b_incy, b_MAX_LAGA_DOFS, b_nb_dofs, b_dime
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate normal jump of test function:
!   jump_norm = jump.norm
!
! --------------------------------------------------------------------------------------------------
!
        b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
        b_nb_dofs = to_blas_int(geom%nb_dofs)
        b_dime = to_blas_int(geom%elem_dime)
!
        jump_norm = 0.d0
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv("N", b_nb_dofs, b_dime, 1.d0, jump_v, &
                   b_MAX_LAGA_DOFS, norm_, b_incx, 0.d0, jump_norm, &
                   b_incy)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
end module
