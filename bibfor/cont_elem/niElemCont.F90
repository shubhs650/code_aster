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
subroutine niElemCont(parameters, geom, nits, coor_qp_sl, hF, &
                      stress_nn, gap, gamma_c, projRmVal, l_cont_qp, &
                      stress_t, vT, gamma_f, projBsVal, l_fric_qp, &
                      dGap, d2Gap, jump_t, dStress_nn)
!
    use contact_module
    use contact_type
    use contact_algebra_module
    use contact_nitsche_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/apnorm.h"
#include "contact_module.h"
#include "MeshTypes_type.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    type(ContactNitsche), intent(in) :: nits
    real(kind=8), intent(in) :: coor_qp_sl(2), hF
    real(kind=8), intent(out) :: stress_nn, gap, gamma_c, projRmVal
    real(kind=8), intent(out) :: stress_t(2), vT(2), gamma_f, projBsVal(2)
    aster_logical, intent(out) :: l_cont_qp, l_fric_qp
    real(kind=8), intent(out), optional :: dGap(MAX_LAGA_DOFS)
    real(kind=8), intent(out), optional :: d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
    real(kind=8), intent(out), optional :: jump_t(MAX_LAGA_DOFS, 3)
    real(kind=8), intent(out), optional :: dStress_nn(MAX_NITS_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute elementary contact quantities
!
! ------------------------------------
! In  elem_dime        : dimension of elements
! In  l_axis           : .true. for axisymmetric element
! In  nb_lagr          : total number of Lagrangian dof on contact element
! In  indi_lagc        : node where Lagrangian dof is present (1) or not (0)
! In  stress_nn           : value of contact lagrangian
! In  l_norm_smooth    : indicator for normals smoothing
! In  nb_node_slav     : number of nodes of for slave side from contact element
! In  elem_slav_code   : code element for slave side from contact element
! In  nb_node_mast     : number of nodes of for master side from contact element
! In  elem_mast_code   : code element for master side from contact element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: shape_func_sl(MT_NNOMAX2D), dshape_func_sl(2, MT_NNOMAX2D)
    real(kind=8) ::  ddshape_func_sl(3, MT_NNOMAX2D)
    real(kind=8) :: shape_func_ma(MT_NNOMAX2D), dshape_func_ma(2, MT_NNOMAX2D)
    real(kind=8) ::  ddshape_func_ma(3, MT_NNOMAX2D)
    real(kind=8) :: shape_func_vo(MT_NNOMAX3D), dshape_func_vo(3, MT_NNOMAX3D)
    real(kind=8) :: norm_slav(3), norm_mast(3), tau_slav(3, 2), tau_mast(3, 2)
    real(kind=8) :: H, coor_qp_ma(2), sn_gap, st_v(3), coor_qp_vo(3)
    real(kind=8) :: metricTens(2, 2), invMetricTens(2, 2)
    real(kind=8) :: metricTens_mast(2, 2), invMetricTens_mast(2, 2), H_mast(2, 2)
    real(kind=8) :: thres_qp, speed(3), projBsVal3(3)
    real(kind=8) :: dNs(MAX_LAGA_DOFS, 3), dGap_(MAX_LAGA_DOFS)
    real(kind=8) :: jump_v(MAX_LAGA_DOFS, 3), dZetaM(MAX_LAGA_DOFS, 2)
    real(kind=8) :: dTs_ns(MAX_LAGA_DOFS, 2), dTm_nm(MAX_LAGA_DOFS, 2)
    real(kind=8) :: dTs(MAX_LAGA_DOFS, 3, 2)
    real(kind=8) :: d2Ns_nm(MAX_LAGA_DOFS, MAX_LAGA_DOFS), dTs_nm(MAX_LAGA_DOFS, 2)
    real(kind=8) :: stress(3, 3), stress_n(3), norm_slav_init(3)
    real(kind=8) :: dStress_n(MAX_NITS_DOFS, 3)
!
! ----- Project quadrature point (on master side)
!
    call projQpSl2Ma(geom, coor_qp_sl, parameters%proj_tole, &
                     coor_qp_ma, gap, tau_slav, norm_slav, tau_mast, norm_mast)
!
! ----- Change parametric coordinates (on slave side)
!
    call projQpSl2Vo(geom, coor_qp_sl, coor_qp_vo)
!
! ----- Evaluate shape function for displacement (slave and master)
!
    call shapeFuncDisp(geom%elem_dime, geom%nb_node_slav, geom%elem_slav_code, coor_qp_sl, &
                       shape_func_sl, dshape_func_sl, ddshape_func_sl)
    call shapeFuncDisp(geom%elem_dime, geom%nb_node_mast, geom%elem_mast_code, coor_qp_ma, &
                       shape_func_ma, dshape_func_ma, ddshape_func_ma)
    call shapeFuncDispVolu(geom%elem_volu_code, coor_qp_vo, shape_func_vo, dshape_func_vo)
!
! ------ Compute outward slave normal (inital configuration)
!
    call apnorm(geom%nb_node_slav, geom%elem_slav_code, geom%elem_dime, geom%coor_slav_init, &
                coor_qp_sl(1), coor_qp_sl(2), norm_slav_init)
!
! ----- Evaluate stress_nn and gamma_c at quadrature point
!
    stress = evalStress(nits, geom%nb_node_slav, shape_func_sl)
    stress_n = matmul(stress, norm_slav_init)
    stress_nn = dot_product(stress_n, norm_slav)
    gamma_c = evalPoly(geom%nb_node_slav, shape_func_sl, parameters%coef_cont)/hF
    sn_gap = stress_nn+gamma_c*gap
!
! ----- Contact activate at quadrature point ( H = 0 or 1 )
!
    if (parameters%type_cont == CONT_TYPE_UNIL) then
        if (parameters%cont_init == PAIR_CONT_INTE) then
            H = Heaviside(-sn_gap)
            projRmVal = projRm(sn_gap)
        elseif (parameters%cont_init == FRIC_ALGO_FALSE) then
            H = 0.d0
            projRmVal = 0.d0
        elseif (parameters%cont_init == FRIC_ALGO_TRUE) then
            H = 1.d0
            projRmVal = sn_gap
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        H = 1.d0
        projRmVal = sn_gap
    end if
    l_cont_qp = (H > 0.5d0)
!
! ----- Evaluate contact derivative
!
    metricTens_mast = metricTensor(tau_mast)
    invMetricTens_mast = invMetricTensor(geom, metricTens_mast)
    H_mast = secondFundForm(geom%elem_dime, geom%nb_node_mast, geom%coor_mast_pair, &
                            ddshape_func_ma, norm_mast)
!
    jump_v = jump(geom, shape_func_sl, shape_func_ma)
    dTs_ns = dTs_du_ns(geom, dshape_func_sl, norm_slav)
    dNs = dNs_du(geom, tau_slav, dTs_ns)
    call dPi_du(geom, norm_slav, tau_mast, gap, jump_v, dNs, dGap_, dZetaM)
!
! ----- Evaluate stress_t and gamma_f at quadrature point
!
    vT = 0.d0
    stress_t = 0.d0
    gamma_f = 0.d0
    projBsVal = 0.d0
    l_fric_qp = ASTER_FALSE
!
    if (parameters%l_fric) then
        call thresEval(parameters, l_cont_qp, projRmVal, thres_qp, l_fric_qp)
!
        metricTens = metricTensor(tau_slav)
        invMetricTens = invMetricTensor(geom, metricTens)
!
        stress_t(1) = 0.d0
        stress_t(2) = 0.d0
!
        gamma_f = evalPoly(geom%nb_node_slav, shape_func_sl, parameters%coef_fric)/hF
!
        speed = speedEval(geom, coor_qp_sl, coor_qp_ma, gap)
        vT = cmpTang(invMetricTens, tau_slav, speed)
!
        st_v = (stress_nn*norm_slav+matmul(tau_slav, stress_t))-gamma_f*speed
        projBsVal3 = projBs(parameters, st_v, thres_qp, norm_slav)
        projBsVal = cmpTang(invMetricTens, tau_slav, projBsVal3)
    end if
!
    if (present(dGap)) then
!
! ----- Compute d (gap(u))[v] / du
!
        dGap = dGap_
!
    end if
!
    if (present(d2Gap)) then
!
! ----- Compute d^2 (gap(u))[v, w] / du^2
!
        dTs = dTs_du(geom, dshape_func_sl)
        dTs_nm = dTs_du_nm(geom, dshape_func_sl, norm_mast)
        dTm_nm = dTm_du_nm(geom, dshape_func_ma, norm_mast)
        d2Ns_nm = d2Ns_du2_nm(geom, tau_slav, norm_mast, dNs, dTs, dTs_ns, dTs_nm)
        d2Gap = d2Gap_du2(geom, norm_slav, norm_mast, gap, H_mast, &
                          dNs, dGap_, dZetaM, dTm_nm, d2Ns_nm)

    end if
!
    if (present(jump_t)) then
!
! ----- Compute tangential jump
!
        if (parameters%l_fric) then
            jump_t = jump_tang(geom, shape_func_sl, shape_func_ma, norm_slav)
        else
            jump_t = 0.d0
        end if
    end if
!
    if (present(dStress_nn)) then
!
! ----- Compute d stress_nn / du
!
        dStress_n = dStress_n_du(geom, nits, norm_slav_init)
        dStress_nn = dStress_nn_du(geom, stress_n, dStress_n, norm_slav, dNs)
    end if
!
end subroutine
