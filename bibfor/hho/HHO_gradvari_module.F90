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
! aslint: disable=W1504

module HHO_GV_module
!
    use NonLin_Datastructure_type
    use Behaviour_type
    use Behaviour_module
    use HHO_basis_module
    use HHO_compor_module
    use HHO_Dirichlet_module
    use HHO_eval_module
    use HHO_LargeStrainMeca_module
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_SmallStrainMeca_module
    use HHO_Ther_module
    use HHO_type
    use HHO_utils_module
    use HHO_algebra_module
    use HHO_matrix_module
    use FE_algebra_module
!
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/deflg4.h"
#include "asterfort/desymt46.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/lagmodtonommod.h"
#include "asterfort/lcdetf.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmcomp.h"
#include "asterfort/pk2topk1.h"
#include "asterfort/poslog.h"
#include "asterfort/prelog.h"
#include "asterfort/prodmt.h"
#include "asterfort/rcvalb.h"
#include "asterfort/readVector.h"
#include "asterfort/sigtopk1.h"
#include "asterfort/tecach.h"
#include "blas/dsyr.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - mechanics - GRAD_VARI
!
! Specific routines for mechanics
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: HHO_GV_State
    public :: hhoGradVariLC, hhoCalcOpGv, hhoDataGVinit, hhoExtrField
    private :: check_behavior, gdef_log, hhoAssGVRhs, hhoAssGVLhs, numGVMap
    private :: initialize_gv, hhoCalcStabCoeffGV
    private :: hhoComputeLhsLargeLM, hhoComputeLhsLargeML
    private :: hhoComputeLhsLargeVM, hhoComputeLhsLargeMV
    private :: hhoComputeLhsSmallLM, hhoComputeLhsSmallML
    private :: hhoComputeLhsSmallVM, hhoComputeLhsSmallMV
!
! --------------------------------------------------------------------------------------------------
!
    type HHO_GV_State
!
        aster_logical :: l_debug = ASTER_FALSE
! ----- GRAD_VARI
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: vari_prev = 0.d0
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: vari_curr = 0.d0
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: vari_incr = 0.d0
! ----- LAGR
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: lagv_prev = 0.d0
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: lagv_curr = 0.d0
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: lagv_incr = 0.d0
! ----- Gradient reconstruction and stabilisation
        type(HHO_matrix) :: grad
        type(HHO_matrix) :: stab
! ----- member function
    contains
        procedure, pass :: initialize => initialize_gv
        procedure, pass :: free => free_gv
    end type HHO_GV_State
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradVariLC(hhoCell, hhoDataMk, hhoDataGv, &
                             hhoQuadCellRigi, hhoMecaState, hhoCS, &
                             hhoGVState, lhs, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoDataMk, hhoDataGv
        type(HHO_Quadrature), intent(in) :: hhoQuadCellRigi
        type(HHO_Meca_State), intent(in) :: hhoMecaState
        type(HHO_Compor_State), intent(inout) :: hhoCS
        type(HHO_GV_State), intent(in) :: hhoGVState
        type(HHO_matrix), intent(out) :: lhs
        real(kind=8), intent(out) :: rhs(MSIZE_TDOFS_MIX)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the local contribution for GRAD_VARI
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   Out lhs         : local contribution (lhs)
!   Out rhs         : local contribution (rhs)
!
!   phi = phi^e(eps) + phi^p(p) + c/2(grad vari, grad_vari) + phi^nl(vari, lagv, p)
!   eps_gene = (eps, vari, lagv, grad vari)
!   sig_gene = (sig, dphi^nl/dvari, dphi^nl/dlagv, c grad vari)
!
!   Attention, il faut que la quadrature soit suffisante pour tout les termes comme
!   (c_phi, c_phi), (c_phi, g_phi), (g_phi, g_phi)
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: ksp = 1
        type(HHO_basis_cell) :: hhoBasisCell
        real(kind=8), dimension(MSIZE_CELL_MAT) :: mk_bT, G_prev_coeff, G_curr_coeff
        real(kind=8), dimension(MSIZE_CELL_VEC) :: gv_bT, GV_prev_coeff, GV_curr_coeff
        real(kind=8) :: G_prev(3, 3), G_curr(3, 3)
        real(kind=8) :: GV_prev(3), GV_curr(3)
        real(kind=8) :: F_prev(3, 3), F_curr(3, 3), Pk1(3, 3)
        real(kind=8) :: var_prev, var_curr, lag_prev, lag_curr
        real(kind=8) :: sig_vari, sig_lagv, sig_gv(3), dsv_dv, dsv_dl, dsl_dl, dsgv_dgv(3, 3)
        real(kind=8) :: dsv_dF(3, 3), dsl_dF(3, 3), dsv_dEps(6), dsl_dEps(6)
        real(kind=8) :: Eps_prev(6), Eps_curr(6), Cauchy(6)
        real(kind=8) :: dPK1_dF(3, 3, 3, 3), dPK1_dv(3, 3), dPK1_dl(3, 3)
        real(kind=8) :: dSig_dEps(6, 6), dSig_dv(6), dSig_dl(6)
        real(kind=8) :: jac_prev, jac_curr, coorpg(3), weight, coeff, mk_stab, gv_stab
        real(kind=8) ::  BSCEval(MSIZE_CELL_SCAL)
        type(HHO_matrix) :: mk_AT, gv_AT, gv_TMP, mv_AT, ml_AT, vm_AT, lm_AT
        type(HHO_matrix) :: lhs_mv, lhs_ml, lhs_mm, lhs_ll, lhs_vm, lhs_vv, lhs_vl, lhs_lm, lhs_lv
        type(HHO_matrix) :: mv_lhs_ax, ml_lhs_ax, vm_lhs_ax, lm_lhs_ax
        real(kind=8) :: rhs_vari(MSIZE_TDOFS_SCAL), rhs_lagv(MSIZE_CELL_SCAL)
        real(kind=8) :: rhs_mk(MSIZE_TDOFS_VEC)
        integer(kind=8) :: mapMeca(MSIZE_TDOFS_VEC), mapVari(MSIZE_TDOFS_SCAL)
        integer(kind=8) :: mapLagv(MSIZE_CELL_SCAL)
        integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs, mk_gbs, mk_gbs_sym, mk_gbs_cmp
        integer(kind=8) :: mk_faces_dofs
        integer(kind=8) :: gv_cbs, gv_fbs, gv_total_dofs, gv_gbs, gv_faces_dofs, gv_cell_offset
        integer(kind=8) :: cod(MSIZE_QP_CELL), kpg, mk_gbs_tot, mk_gbs_axis
        aster_logical :: l_lhs, l_rhs, forc_noda
        blas_int :: b_n
        blas_int, parameter :: b_one = to_blas_int(1)
! --------------------------------------------------------------------------------------------------
!
        cod = 0
!
! ----- Type of behavior
        call check_behavior(hhoCS)
!
! ----- Vector and/or matrix
        forc_noda = hhoCS%option == "FORC_NODA"
        l_lhs = L_MATR(hhoCS%option)
        l_rhs = L_VECT(hhoCS%option) .or. forc_noda
!
! ------ number of dofs
!
        call hhoMecaNLDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs, &
                           mk_gbs, mk_gbs_sym, mk_gbs_axis)
        call hhoTherNLDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs, gv_gbs)
        gv_faces_dofs = gv_total_dofs-gv_cbs
        gv_cell_offset = gv_faces_dofs+1
        if (hhoCS%l_largestrain) then
            mk_gbs_tot = mk_gbs
        else
            mk_gbs_tot = mk_gbs_sym
        end if
        mk_faces_dofs = mk_total_dofs-mk_cbs
        mk_gbs_cmp = (mk_gbs-mk_gbs_axis)/(hhoCell%ndim*hhoCell%ndim)
!
! ------ initialization
!
        if (l_lhs) then
            call lhs_mm%initialize(mk_total_dofs, mk_total_dofs, 0.d0)
            call lhs_ll%initialize(gv_cbs, gv_cbs, 0.d0)
            call lhs_mv%initialize(mk_total_dofs, gv_total_dofs, 0.d0)
            call lhs_ml%initialize(mk_total_dofs, gv_cbs, 0.d0)
            call lhs_vm%initialize(gv_total_dofs, mk_total_dofs, 0.d0)
            call lhs_vv%initialize(gv_total_dofs, gv_total_dofs, 0.d0)
            call lhs_vl%initialize(gv_total_dofs, gv_cbs, 0.d0)
            call lhs_lm%initialize(gv_cbs, mk_total_dofs, 0.d0)
            call lhs_lv%initialize(gv_cbs, gv_total_dofs, 0.d0)
            call mk_AT%initialize(mk_gbs_tot, mk_gbs_tot, 0.d0)
            call gv_AT%initialize(gv_gbs, gv_gbs, 0.d0)
            call gv_TMP%initialize(gv_gbs, gv_total_dofs, 0.d0)
            call mv_AT%initialize(mk_gbs_tot, gv_total_dofs, 0.d0)
            call ml_AT%initialize(mk_gbs_tot, gv_cbs, 0.d0)
            call vm_AT%initialize(gv_total_dofs, mk_gbs_tot, 0.d0)
            call lm_AT%initialize(gv_cbs, mk_gbs_tot, 0.d0)
        end if
!
        mk_bT = 0.d0
        gv_bT = 0.d0
        G_prev_coeff = 0.d0
        G_curr_coeff = 0.d0
        GV_prev_coeff = 0.d0
        GV_curr_coeff = 0.d0
        rhs_vari = 0.d0
        rhs_lagv = 0.d0
        rhs_mk = 0.d0

! ----- Prepare external state variables (geometry)
        if (.not. forc_noda) then
            call behaviourPrepESVAGeomHHO(hhoCell, hhoQuadCellRigi, hhoCS%BEHInteg)
        end if
!
! ----- init basis
        call hhoBasisCell%initialize(hhoCell)
!
! ----- compute G_prev = gradrec * depl_prev (sym or not)
!
        call hhoMecaState%grad%dot(hhoMecaState%depl_prev, G_prev_coeff)
!
! ----- compute G_curr = gradrec * depl_curr (sym or not)
!
        call hhoMecaState%grad%dot(hhoMecaState%depl_curr, G_curr_coeff)
!
! ----- compute GV_prev = gradrec * vari_prev
!
        call hhoGVState%grad%dot(hhoGVState%vari_prev, GV_prev_coeff)
!
! ----- compute GV_curr = gradrec * vari_curr
!
        call hhoGVState%grad%dot(hhoGVState%vari_curr, GV_curr_coeff)
!
! ----- Loop on quadrature point
!
        do kpg = 1, hhoQuadCellRigi%nbQuadPoints
            coorpg(1:3) = hhoQuadCellRigi%points(1:3, kpg)
            weight = hhoQuadCellRigi%weights(kpg)

! --------- Eval basis function at the quadrature point
            call hhoBasisCell%BSEval(coorpg(1:3), 0, &
                                     max(hhoDataMk%grad_degree(), hhoDataMk%cell_degree(), &
                                         hhoDataGv%grad_degree(), hhoDataGv%cell_degree()), &
                                     BSCEval)

! --------- Eval gradient at T- and T+
            if (hhoCS%l_largestrain) then
                G_prev = hhoEvalMatCell(hhoCell%ndim, mk_gbs, BSCEval, G_prev_coeff)
                G_curr = hhoEvalMatCell(hhoCell%ndim, mk_gbs, BSCEval, G_curr_coeff)
!
! ------------- Eval gradient of the deformation at T- and T+
                call hhoCalculF(G_prev, F_prev)
                call hhoCalculF(G_curr, F_curr)

! ------------- Check the jacobian jac >= r8prem
! ------------- be carreful with c_plan, I don't know the result
                call lcdetf(hhoCell%ndim, F_prev, jac_prev)
                cod(kpg) = merge(1, 0, jac_prev .le. r8prem())
                if (cod(kpg) .ne. 0) goto 999
!
                call lcdetf(hhoCell%ndim, F_curr, jac_curr)
                cod(kpg) = merge(1, 0, jac_curr .le. r8prem())
                if (cod(kpg) .ne. 0) goto 999
            else
                Eps_prev = hhoEvalSymMatCell(hhoCell%ndim, mk_gbs_sym, BSCEval, G_prev_coeff)
                Eps_curr = hhoEvalSymMatCell(hhoCell%ndim, mk_gbs_sym, BSCEval, G_curr_coeff)
            end if
!
            GV_prev = hhoEvalVecCell(hhoCell%ndim, gv_gbs, BSCEval, GV_prev_coeff)
            GV_curr = hhoEvalVecCell(hhoCell%ndim, gv_gbs, BSCEval, GV_curr_coeff)
!
            var_prev = hhoEvalScalCell(gv_cbs, BSCEval, hhoGVState%vari_prev(gv_cell_offset:))
            var_curr = hhoEvalScalCell(gv_cbs, BSCEval, hhoGVState%vari_curr(gv_cell_offset:))
!
            lag_prev = hhoEvalScalCell(gv_cbs, BSCEval, hhoGVState%lagv_prev)
            lag_curr = hhoEvalScalCell(gv_cbs, BSCEval, hhoGVState%lagv_curr)

! --------- Set main parameters for behaviour (on point)
            call behaviourSetParaPoin(kpg, ksp, hhoCS%BEHInteg)

! --------- Compute behavior
            if (forc_noda) then
                call forc_noda_stress(hhoCS, hhoCell%ndim, kpg, F_curr, Pk1, &
                                      Cauchy, sig_vari, sig_lagv, sig_gv)
            else if (hhoCS%l_largestrain) then
                call gdef_log(hhoCS, hhoCell%ndim, kpg, &
                              hhoMecaState%time_prev, hhoMecaState%time_curr, F_prev, F_curr, &
                              var_prev, var_curr, lag_prev, lag_curr, GV_prev, &
                              GV_curr, Pk1, sig_vari, sig_lagv, sig_gv, &
                              dPK1_dF, dPK1_dv, dPK1_dl, dsv_dF, dsv_dv, &
                              dsv_dl, dsl_dF, dsl_dl, dsgv_dgv, cod(kpg))
            else
                call petit(hhoCS, hhoCell%ndim, kpg, hhoMecaState%time_prev, &
                           hhoMecaState%time_curr, Eps_prev, Eps_curr, var_prev, var_curr, &
                           lag_prev, lag_curr, GV_prev, GV_curr, Cauchy, &
                           sig_vari, sig_lagv, sig_gv, dSig_dEps, dSig_dv, &
                           dSig_dl, dsv_dEps, dsv_dv, dsv_dl, dsl_dEps, &
                           dsl_dl, dsgv_dgv, cod(kpg))
            end if
!
! -------- Test the code of the LDC
!
            if (cod(kpg) .eq. 1) goto 999
!
! ------- Compute rhs
!
            if (l_rhs) then
                if (hhoCS%l_largestrain) then
! ---------- += weight * (PK1, g_phi)
                    call hhoComputeRhsLarge(hhoCell, Pk1, weight, BSCEval, mk_gbs, mk_bT)
                else
! ---------- += weight * (Cauchy, gs_phi)
                    call hhoComputeRhsSmall(hhoCell, Cauchy, weight, BSCEval, mk_gbs_cmp, mk_bT)
                end if
! ---------- += weight * (sig_gv, g_phi)
                call hhoComputeRhsRigiTher(hhoCell, sig_gv, weight, BSCEval, gv_gbs, gv_bT)
! ---------- += weight * (sig_vari, c_phi)
                call hhoComputeRhsMassTher(sig_vari, weight, BSCEval, gv_cbs, &
                                           rhs_vari(gv_cell_offset:))
! ---------- += weight * (sig_lagv, c_phi)
                call hhoComputeRhsMassTher(sig_lagv, weight, BSCEval, gv_cbs, rhs_lagv)
            end if
!
! ------- Compute lhs
!
            if (l_lhs) then
                if (hhoCS%l_largestrain) then
! ---------- += weight * (dPK1_dF : g_phi, g_phi)
                    call hhoComputeLhsLarge(hhoCell, dPK1_dF, weight, BSCEval, mk_gbs, &
                                            mk_AT)
!
! ---------- += weight * (g_phi, dPK1_dv : c_phi) -> lhs_mv
                    call hhoComputeLhsLargeMV(hhoCell, dPK1_dv, weight, BSCEval, gv_cbs, &
                                              mk_gbs, mv_AT)
! ---------- += weight * (g_phi, dPK1_dl : c_phi) -> lhs_ml
                    call hhoComputeLhsLargeML(hhoCell, dPK1_dl, weight, BSCEval, gv_cbs, &
                                              mk_gbs, ml_AT)
! ---------- += weight * (dsv_dF : g_phi, c_phi) -> lhs_vm
                    call hhoComputeLhsLargeVM(hhoCell, dsv_dF, weight, BSCEval, gv_cbs, &
                                              mk_gbs, vm_AT)
! ---------- += weight * (dsl_dF : g_phi, c_phi) -> lhs_lm
                    call hhoComputeLhsLargeLM(hhoCell, dsl_dF, weight, BSCEval, gv_cbs, &
                                              mk_gbs, lm_AT)
                else
! ---------- += weight * (dSig_deps : gs_phi, gs_phi)
                    call hhoComputeLhsSmall(hhoCell, dSig_deps, hhoCS%matsym, weight, &
                                            BSCEval, mk_gbs_sym, mk_gbs_cmp, mk_AT)
!
! ---------- += weight * (gs_phi, dSig_dv : c_phi) -> lhs_mv
                    call hhoComputeLhsSmallMV(hhoCell, dSig_dv, weight, BSCEval, gv_cbs, &
                                              mk_gbs_cmp, mv_AT)
! ---------- += weight * (gs_phi, dSig_dl : c_phi) -> lhs_ml
                    call hhoComputeLhsSmallML(hhoCell, dSig_dl, weight, BSCEval, gv_cbs, &
                                              mk_gbs_cmp, ml_AT)
! ---------- += weight * (dsv_dEps : gs_phi, c_phi) -> lhs_vm
                    call hhoComputeLhsSmallVM(hhoCell, dsv_dEps, weight, BSCEval, gv_cbs, &
                                              mk_gbs_sym, mk_gbs_cmp, vm_AT)
! ---------- += weight * (dsl_dEps : gs_phi, c_phi) -> lhs_lm
                    call hhoComputeLhsSmallLM(hhoCell, dsl_dEps, weight, BSCEval, gv_cbs, &
                                              mk_gbs_sym, mk_gbs_cmp, lm_AT)
                end if
! ---------- += weight * (dgv_dv : g_phi, g_phi)
                call hhoComputeLhsRigiTher(hhoCell, dsgv_dgv, weight, BSCEval, gv_gbs, gv_AT)
! ---------- += weight * (dsv_dv : c_phi, c_phi)
                coeff = weight*dsv_dv
                b_n = to_blas_int(gv_cbs)
                call dsyr('U', b_n, coeff, BSCEval, b_one, &
                          lhs_vv%m(gv_cell_offset:gv_total_dofs, gv_cell_offset:gv_total_dofs), &
                          b_n)
! ---------- += weight * (dsv_dl : c_phi, c_phi)
                coeff = weight*dsv_dl
                call dsyr('U', b_n, coeff, BSCEval, b_one, &
                          lhs_vl%m(gv_cell_offset:gv_total_dofs, 1:gv_cbs), b_n)
! ---------- += weight * (dsl_dl : c_phi, c_phi)
                call hhoComputeLhsMassTher(dsl_dl, weight, BSCEval, gv_cbs, lhs_ll)
            end if
!
        end do
!
999     continue
!
! - SYNTHESE DES CODES RETOURS
!
        call codere(cod, hhoQuadCellRigi%nbQuadPoints, hhoCS%codret)
!
        call numGVMap(hhoCell, hhoDataMk, hhoDataGv, mapMeca, mapVari, mapLagv)
        call hhoCalcStabCoeffMeca(hhoDataMk, &
                                  hhoCS%BEHInteg%materPara%schemePara%fami, &
                                  hhoMecaState%time_curr, &
                                  hhoQuadCellRigi)
        mk_stab = hhoDataMk%coeff_stab()
        gv_stab = hhoCalcStabCoeffGV(hhoCS%BEHInteg%materPara%schemePara%fami)
!
! ------- Compute rhs
!
        if (l_rhs) then
! ----- compute rhs += Gradrec**T * bT
            call hho_dgemv_T(1.d0, hhoMecaState%grad, mk_bT, 1.d0, rhs_mk)
! ----- compute rhs += stab
            call hho_dsymv_U(mk_stab, hhoMecaState%stab, hhoMecaState%depl_curr, 1.d0, rhs_mk)
! ----- compute rhs += Gradrec**T * bT
            call hho_dgemv_T(1.d0, hhoGVState%grad, gv_bT, 1.d0, rhs_vari)
! ----- compute rhs += stab
            call hho_dsymv_U(gv_stab, hhoGVState%stab, hhoGVState%vari_curr, 1.d0, rhs_vari)
! ----- assembly
            call hhoAssGVRhs(hhoCell, hhoDataMk, hhoDataGv, mapMeca, mapVari, mapLagv, &
                             rhs_mk, rhs_vari, rhs_lagv, rhs)
        end if
!
! ------- Compute lhs
!
        if (l_lhs) then
! ----- Add symmetry
!
            call lhs_vv%copySymU(gv_faces_dofs, gv_faces_dofs)
            call lhs_vl%copySymU(gv_faces_dofs, 0)
            call lhs_ll%copySymU()
!
! ----- Add gradient: += gradrec**T * AT * gradrec
! ----- step1: TMP = AT * gradrec
            if (hhoCS%l_largestrain) then
                call hhoAssembleLhsLarge(hhoMecaState%grad, &
                                         mk_AT, lhs_mm)
            else
                call hhoAssembleLhsSmall(hhoCS, hhoMecaState%grad, mk_AT, lhs_mm)
            end if
!
            call hho_dgemm_NN(1.d0, gv_AT, hhoGVState%grad, 0.d0, gv_TMP)
! ----- step2: lhs += gradrec**T * TMP
!
            call hho_dgemm_TN(1.d0, hhoGVState%grad, gv_TMP, 1.d0, lhs_vv)
!
            call hho_dgemm_TN(1.d0, hhoMecaState%grad, mv_AT, 0.d0, lhs_mv)
!
            call hho_dgemm_TN(1.d0, hhoMecaState%grad, ml_AT, 0.d0, lhs_ml)
!
            call hho_dgemm_NN(1.d0, vm_AT, hhoMecaState%grad, 0.d0, lhs_vm)
!
            call hho_dgemm_NN(1.d0, lm_AT, hhoMecaState%grad, 0.d0, lhs_lm)
!
! ----- Add stabilization
! ----- += coeff * stab_mk
            call lhs_mm%add(hhoMecaState%stab, mk_stab)
! ----- += coeff * stab_vv
            call lhs_vv%add(hhoGVState%stab, gv_stab)
!
! ----- the symmetry is checked inside behaviour
            lhs_lv%m = transpose(lhs_vl%m)
!
! ----- assembly
            call hhoAssGVLhs(hhoCell, hhoDataMk, hhoDataGv, mapMeca, mapVari, mapLagv, &
                             lhs_mm, lhs_mv, lhs_ml, lhs_vm, lhs_vv, &
                             lhs_vl, lhs_lm, lhs_lv, lhs_ll, lhs)
        end if
!
        call lhs_mm%free()
        call lhs_ll%free()
        call lhs_mv%free()
        call lhs_ml%free()
        call lhs_vm%free()
        call lhs_vv%free()
        call lhs_vl%free()
        call lhs_lm%free()
        call lhs_lv%free()
        call gv_AT%free()
        call gv_TMP%free()
        call mv_AT%free()
        call ml_AT%free()
        call vm_AT%free()
        call lm_AT%free()
        call mv_lhs_ax%free()
        call ml_lhs_ax%free()
        call vm_lhs_ax%free()
        call lm_lhs_ax%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine check_behavior(hhoCS)
!
        implicit none
!
        type(HHO_Compor_State), intent(in) :: hhoCS
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Check behavior
!   In behavior     : type of behavior
! --------------------------------------------------------------------------------------------------
!
        select case (hhoCS%compor(DEFO))
        case ('GDEF_LOG')
            ASSERT(hhoCS%l_largestrain)
        case ('PETIT')
            ASSERT(.not. hhoCS%l_largestrain)
        case default
            ASSERT(ASTER_FALSE)
        end select
        ASSERT(.not. hhoCS%c_plan)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine gdef_log(hhoCS, ndim, kpg, time_prev, &
                        time_curr, F_prev, F_curr, var_prev, var_curr, &
                        lag_prev, lag_curr, GV_prev, GV_curr, PK1_curr, &
                        sig_vari, sig_lagv, sig_gv, dPK1_dF, dPK1_dv, &
                        dPK1_dl, dsv_dF, dsv_dv, dsv_dl, dsl_dF, &
                        dsl_dl, dsgv_dgv, cod)
!
        implicit none
!
        type(HHO_Compor_State), intent(inout) :: hhoCS
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: time_prev, time_curr
        real(kind=8), intent(in) :: F_prev(3, 3)
        real(kind=8), intent(in) :: F_curr(3, 3)
        real(kind=8), intent(in) :: var_prev, var_curr
        real(kind=8), intent(in) :: lag_prev, lag_curr
        real(kind=8), intent(in) :: GV_prev(3), GV_curr(3)
        real(kind=8), intent(out) :: PK1_curr(3, 3)
        real(kind=8), intent(out) :: sig_vari, sig_lagv, sig_gv(3)
        real(kind=8), intent(out) :: dsv_dv, dsv_dl, dsl_dl, dsgv_dgv(3, 3)
        real(kind=8), intent(out) :: dPK1_dF(3, 3, 3, 3), dPK1_dv(3, 3), dPK1_dl(3, 3)
        real(kind=8), intent(out) :: dsv_dF(3, 3), dsl_dF(3, 3)
        integer(kind=8), intent(out) :: cod
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the behavior laws for GDEF_LOG
!   In ndim         : dimension of the problem
!   In fami         : familly of quadrature points
!   In typmod       : type of modelization
!   In compor       : type of behavior
!   In option       : option of computations
!   In carcri       : local criterion of convergence
!   In lgpg         : size of internal variables for 1 pg
!   In kpg          : i-th quadrature point
!   In time_prev    : previous time T-
!   In time_curr    : current time T+
!   In multcomp     : ?
!   In cplan        : plane stress hypothesis
!   In F_prev       : previous deformation gradient at T-
!   In F_curr       : curr deformation gradient at T+
!   In sig_prev_pg  : cauchy stress at T-  (XX, YY, ZZ, XY, XZ, YZ)
!   In vi_prev_pg   : internal variables at T-
!   Out sig_curr    : cauchy stress at T+  (XX, YY, ZZ, XY, XZ, YZ)
!   Out vi_curr     : internal variables at T+
!   Out Pk1_curr    : piola-kirschooff 1 at T+
!   Out module_tang : tangent modulus dPK1/dF(Fp)
!   Out cod         : info on integration of the LDC
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: gn(3, 3), lamb(3), logl(3), epslPrev(6), epslIncr(6)
        real(kind=8) :: tlogPrev(6), tlogCurr(6)
        real(kind=8) :: dtde(6, 6), PK2_prev(6), PK2_curr(6), sig(6)
        real(kind=8) :: dpk2dc(6, 6), ftot(3, 3, 3, 3)
        real(Kind=8) :: eplcm(3*ndim+2), eplci(3*ndim+2)
        real(Kind=8) :: silcm(3*ndim+2), silcp(3*ndim+2), dsde(3*ndim+2, 3*ndim+2)
        real(kind=8) :: sigPrev(11), viPrev(hhoCS%lgpg), viCurr(hhoCS%lgpg)
        real(kind=8) :: dT_dv(6), dT_dl(6)
        real(kind=8) :: dsv_de(6), dsl_de(6)
        real(kind=8) :: dsl_dv, norm, pe(3, 3, 3, 3)
        integer(kind=8) :: lgpg, neu, neg, ntot
        aster_logical :: lCorr, lMatr, lSigm, lVari
!
        lCorr = L_CORR(hhoCS%option)
        lMatr = L_MATR(hhoCS%option)
        lSigm = L_SIGM(hhoCS%option)
        lVari = L_VARI(hhoCS%option)
!
        neu = 2*ndim
        neg = 2+ndim
        ntot = neu+neg
!
        lgpg = hhoCS%lgpg
        sigPrev(1:ntot) = hhoCS%sig_prev((kpg-1)*ntot+1:kpg*ntot)
        viPrev = hhoCS%vari_prev((kpg-1)*lgpg+1:kpg*lgpg)

! ----- Compute pre-processing Elog
        call prelog(ndim, lgpg, viPrev, gn, lamb, &
                    logl, F_prev, F_curr, epslPrev, epslIncr, &
                    tlogPrev, lCorr, cod)
        if (cod .ne. 0) then
            goto 999
        end if
! Preparation des deformations generalisees de ldc en t- et t+
        eplcm(1:neu) = epslPrev(1:neu)
        eplci(1:neu) = epslIncr(1:neu)
        eplcm(neu+1) = var_prev
        eplci(neu+1) = var_curr-var_prev
        eplcm(neu+2) = lag_prev
        eplci(neu+2) = lag_curr-lag_prev
        eplcm(neu+2+1:neu+2+ndim) = GV_prev(1:ndim)
        eplci(neu+2+1:neu+2+ndim) = GV_curr(1:ndim)-GV_prev(1:ndim)
! Preparation des contraintes generalisees de ldc en t-
        silcm(1:neu) = tlogPrev(1:neu)
        silcm(neu+1:ntot) = sigPrev(neu+1:ntot)

! ----- Compute Stress and module_tangent
        silcp = 0.d0
        viCurr = 0.d0
        dsde = 0.d0
        call nmcomp(hhoCS%BEHInteg, &
                    ndim, hhoCS%option, hhoCS%typmod, &
                    time_prev, time_curr, &
                    hhoCS%compor, hhoCS%carcri, hhoCS%multComp, &
                    ntot, eplcm, eplci, &
                    ntot, silcm, &
                    viPrev, &
                    silcp, viCurr, &
                    ntot*ntot, dsde, cod)
!
! ----- Test the code of the LDC
!
        if (cod .eq. 1) goto 999
!
! ----- Compute post-processing Elog
!
        dtde = 0.d0
        dtde(1:neu, 1:neu) = dsde(1:neu, 1:neu)
        tlogCurr = 0.d0
        tlogCurr(1:neu) = silcp(1:neu)
!
        call poslog(hhoCS%BEHInteg, &
                    lCorr, lMatr, lSigm, lVari, tlogPrev, &
                    tlogCurr, F_prev, lgpg, viCurr, ndim, &
                    F_curr, dtde, sigPrev, hhoCS%c_plan, &
                    time_curr, gn, &
                    lamb, logl, sig, dpk2dc, PK2_prev, &
                    PK2_curr, cod)
!
! Archivage des contraintes mecaniques en t+ (tau tilda) dans les vi

        if (lVari) then
            hhoCS%vari_curr((kpg-1)*lgpg+1:kpg*lgpg) = viCurr
        end if
!
! ----- Test the code of the LDC
!
        if (cod .ne. 0) goto 999
!
        if (.not. lCorr) then
            PK2_curr = PK2_prev
            tlogCurr = tlogPrev
        end if
!
        if (lSigm) then
            hhoCS%sig_curr((kpg-1)*ntot+1:(kpg-1)*ntot+neu) = sig(1:neu)
            hhoCS%sig_curr((kpg-1)*ntot+neu+1:kpg*ntot) = silcp(neu+1:ntot)
        end if
!
! ----- Compute PK1
!
        call pk2topk1(ndim, PK2_curr, F_curr, Pk1_curr)
!
        sig_vari = silcp(neu+1)
        sig_lagv = silcp(neu+2)
        sig_gv = 0.d0
        sig_gv(1:ndim) = silcp(neu+2+1:ntot)
!
        dPK1_dF = 0.d0
        dPK1_dv = 0.d0
        dPK1_dl = 0.d0
        dsv_dF = 0.d0
        dsv_dv = 0.d0
        dsv_dl = 0.d0
        dsl_dF = 0.d0
        dsl_dl = 0.d0
        dsgv_dgv = 0.d0
!
        if (lMatr) then
!
! ----- Unpack lagrangian tangent modulus
!
            call desymt46(dpk2dc, ftot)
!
! ----- Compute nominal tangent modulus
!
            call lagmodtonommod(ftot, PK2_curr, F_curr, dPK1_dF)
!
            dT_dv = 0.d0
            dT_dl = 0.d0
            dsv_de = 0.d0
            dsl_de = 0.d0
!
            dsv_dv = dsde(neu+1, neu+1)
            dsv_dl = dsde(neu+1, neu+2)
            dsl_dv = dsde(neu+2, neu+1)
            dsl_dl = dsde(neu+2, neu+2)
            dT_dv(1:neu) = dsde(1:neu, neu+1)
            dT_dl(1:neu) = dsde(1:neu, neu+2)
            dsv_de(1:neu) = dsde(neu+1, 1:neu)
            dsl_de(1:neu) = dsde(neu+2, 1:neu)
            dsgv_dgv(1:ndim, 1:ndim) = dsde(neu+3:ntot, neu+3:ntot)
!
            call deflg4(gn, lamb, logl, F_curr, pe)
!
!           dP_da = dT_da * dElog_dF
            call prodmt(dT_dv, pe, dPK1_dv)
            call prodmt(dT_dl, pe, dPK1_dl)
!
!           da_dF = da_dElog * dElog_dF
            call prodmt(dsv_de, pe, dsv_dF)
            call prodmt(dsl_de, pe, dsl_dF)
!
! ----- Verify symmetry
            norm = max(1.d0, dsv_dl)
            ASSERT(abs(dsv_dl-dsl_dv) < 1d-10*norm)
        end if
! print *, hhoCS%option
! print *, kpg, lgpg, ntot, viCurr(1)
! print *, eplcm
! print *, eplci
! print*, var_curr, lag_curr, sig_gv
! print*, lag_curr + (var_curr-viCurr(1)) , " vs ", sig_vari
! print*, (var_curr-viCurr(1)) , " vs ", sig_lagv
! print*, GV_curr, " vs ", sig_gv
! print *, sig_vari, sig_lagv, sig_gv
! print *, dsv_dv, dsv_dl, dsl_dl, dsgv_dgv
! print*, kpg, dsv_dv, dsv_dl, dsl_dl
! print*, dPK1_dv
! print*, dPK1_dl
! print*, dsv_dF
! print*, dsl_dF
!
999     continue
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine petit(hhoCS, ndim, kpg, time_prev, &
                     time_curr, Eps_prev, Eps_curr, var_prev, var_curr, &
                     lag_prev, lag_curr, GV_prev, GV_curr, Sig_curr, &
                     sig_vari, sig_lagv, sig_gv, dSig_dEps, dSig_dv, &
                     dSig_dl, dsv_dEps, dsv_dv, dsv_dl, dsl_dEps, &
                     dsl_dl, dsgv_dgv, cod)
!
        implicit none
!
        type(HHO_Compor_State), intent(inout) :: hhoCS
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: time_prev
        real(kind=8), intent(in) :: time_curr
        real(kind=8), intent(in) :: Eps_prev(6)
        real(kind=8), intent(in) :: Eps_curr(6)
        real(kind=8), intent(in) :: var_prev, var_curr
        real(kind=8), intent(in) :: lag_prev, lag_curr
        real(kind=8), intent(in) :: GV_prev(3), GV_curr(3)
        real(kind=8), intent(out) :: Sig_curr(6)
        real(kind=8), intent(out) :: sig_vari, sig_lagv, sig_gv(3)
        real(kind=8), intent(out) :: dsv_dv, dsv_dl, dsl_dl, dsgv_dgv(3, 3)
        real(kind=8), intent(out) :: dSig_dEps(6, 6), dSig_dv(6), dSig_dl(6)
        real(kind=8), intent(out) :: dsv_dEps(6), dsl_dEps(6)
        integer(kind=8), intent(out) :: cod
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the behavior laws for GDEF_LOG
!   In ndim         : dimension of the problem
!   In fami         : familly of quadrature points
!   In typmod       : type of modelization
!   In compor       : type of behavior
!   In option       : option of computations
!   In carcri       : local criterion of convergence
!   In lgpg         : size of internal variables for 1 pg
!   In kpg          : i-th quadrature point
!   In time_prev    : previous time T-
!   In time_curr    : current time T+
!   In multcomp     : ?
!   In cplan        : plane stress hypothesis
!   In F_prev       : previous deformation gradient at T-
!   In F_curr       : curr deformation gradient at T+
!   In sig_prev_pg  : cauchy stress at T-  (XX, YY, ZZ, XY, XZ, YZ)
!   In vi_prev_pg   : internal variables at T-
!   Out sig_curr    : cauchy stress at T+  (XX, YY, ZZ, XY, XZ, YZ)
!   Out vi_curr     : internal variables at T+
!   Out Pk1_curr    : piola-kirschooff 1 at T+
!   Out module_tang : tangent modulus dPK1/dF(Fp)
!   Out cod         : info on integration of the LDC
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: Cauchy_prev(6), Cauchy_curr(6), sig(6)
        real(Kind=8) :: eplcm(3*ndim+2), eplci(3*ndim+2)
        real(Kind=8) :: silcm(3*ndim+2), silcp(3*ndim+2), dsde(3*ndim+2, 3*ndim+2)
        real(kind=8) :: sigPrev(11), viPrev(hhoCS%lgpg), viCurr(hhoCS%lgpg)
        real(kind=8) :: dsl_dv, norm
        integer(kind=8) :: lgpg, neu, neg, ntot
        aster_logical :: lCorr, lMatr, lSigm, lVari
!
        lCorr = L_CORR(hhoCS%option)
        lMatr = L_MATR(hhoCS%option)
        lSigm = L_SIGM(hhoCS%option)
        lVari = L_VARI(hhoCS%option)
!
        neu = 2*ndim
        neg = 2+ndim
        ntot = neu+neg
!
        lgpg = hhoCS%lgpg
        sigPrev(1:ntot) = hhoCS%sig_prev((kpg-1)*ntot+1:kpg*ntot)
        call tranfoMatToSym(ndim, sigPrev(1:neu), Cauchy_prev)
        viPrev = hhoCS%vari_prev((kpg-1)*lgpg+1:kpg*lgpg)
! Preparation des deformations generalisees de ldc en t- et t+
        eplcm(1:neu) = Eps_prev(1:neu)
        eplci(1:neu) = Eps_curr(1:neu)-Eps_prev(1:neu)
        eplcm(neu+1) = var_prev
        eplci(neu+1) = var_curr-var_prev
        eplcm(neu+2) = lag_prev
        eplci(neu+2) = lag_curr-lag_prev
        eplcm(neu+2+1:neu+2+ndim) = GV_prev(1:ndim)
        eplci(neu+2+1:neu+2+ndim) = GV_curr(1:ndim)-GV_prev(1:ndim)
! Preparation des contraintes generalisees de ldc en t-
        silcm(1:neu) = Cauchy_prev(1:neu)
        silcm(neu+1:ntot) = sigPrev(neu+1:ntot)
!
! ----- Compute Stress and module_tangent
!
        silcp = 0.d0
        viCurr = 0.d0
        dsde = 0.d0
        call nmcomp(hhoCS%BEHInteg, &
                    ndim, hhoCS%option, hhoCS%typmod, &
                    time_prev, time_curr, &
                    hhoCS%compor, hhoCS%carcri, hhoCS%multComp, &
                    ntot, eplcm, eplci, &
                    ntot, silcm, &
                    viPrev, &
                    silcp, viCurr, &
                    ntot*ntot, dsde, cod)

! ----- Test the code of the LDC
        if (cod .eq. 1) goto 999

! Archivage des contraintes mecaniques en t+ (tau tilda) dans les vi
        if (lVari) then
            hhoCS%vari_curr((kpg-1)*lgpg+1:kpg*lgpg) = viCurr
        end if
!
! ----- Compute post-processing Elog
!
        Cauchy_curr = 0.d0
        Cauchy_curr(1:neu) = silcp(1:neu)
!
! --------- For new prediction and nmisot.F90
        if (L_PRED(hhoCS%option)) then
            Cauchy_curr = 0.d0
        end if
!
        if (.not. lCorr) then
            Cauchy_curr = Cauchy_prev
        end if
!
        if (lSigm) then
            call tranfoSymToMat(ndim, Cauchy_curr, sig)
            hhoCS%sig_curr((kpg-1)*ntot+1:(kpg-1)*ntot+neu) = sig(1:neu)
            hhoCS%sig_curr((kpg-1)*ntot+neu+1:kpg*ntot) = silcp(neu+1:ntot)
        end if
!
! ----- Compute stress
!
        Sig_curr = 0.d0
        Sig_curr(1:neu) = Cauchy_curr(1:neu)
        sig_vari = silcp(neu+1)
        sig_lagv = silcp(neu+2)
        sig_gv = 0.d0
        sig_gv(1:ndim) = silcp(neu+2+1:ntot)
!
        dSig_dEps = 0.d0
        dSig_dv = 0.d0
        dSig_dl = 0.d0
        dsv_dEps = 0.d0
        dsv_dv = 0.d0
        dsv_dl = 0.d0
        dsl_dEps = 0.d0
        dsl_dl = 0.d0
        dsgv_dgv = 0.d0
!
        if (lMatr) then
!
            dSig_dEps(1:neu, 1:neu) = dsde(1:neu, 1:neu)
            dsv_dv = dsde(neu+1, neu+1)
            dsv_dl = dsde(neu+1, neu+2)
            dsl_dv = dsde(neu+2, neu+1)
            dsl_dl = dsde(neu+2, neu+2)
            dSig_dv(1:neu) = dsde(1:neu, neu+1)
            dSig_dl(1:neu) = dsde(1:neu, neu+2)
            dsv_dEps(1:neu) = dsde(neu+1, 1:neu)
            dsl_dEps(1:neu) = dsde(neu+2, 1:neu)
            dsgv_dgv(1:ndim, 1:ndim) = dsde(neu+3:ntot, neu+3:ntot)
!
! print *, hhoCS%option
! print *, kpg, lgpg, ntot
! print *, eplcm
! print *, eplci
! print *, sig_vari, sig_lagv, sig_gv
! print *, dsv_dv, dsv_dl, dsl_dl, dsgv_dgv
!if (abs(dsl_dl) < 1d-8) dsl_dl = 1.d0
! ----- Verify symmetry
            norm = max(1.d0, dsv_dl)
            ASSERT(abs(dsv_dl-dsl_dv) < 1d-10*norm)
        end if
!
999     continue
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine numGVMap(hhoCell, hhoDataMk, hhoDataGv, mapMeca, mapVari, mapLagv)
!
        implicit none
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoDataMk, hhoDataGv
        integer(kind=8), intent(out) :: mapMeca(MSIZE_TDOFS_VEC)
        integer(kind=8), intent(out) :: mapVari(MSIZE_TDOFS_SCAL)
        integer(kind=8), intent(out) :: mapLagv(MSIZE_CELL_SCAL)
!
!--------------------------------------------------------------------------------------------------
! Numbering map for GRAD_VARI
!
!--------------------------------------------------------------------------------------------------
        integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs, gv_cbs, gv_fbs, gv_total_dofs
        integer(kind=8) :: i_face, i_dof, num_tot, num_gv, num_vari, num_mk
!
!
        call hhoMecaDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs)
        call hhoTherDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs)
!
        mapMeca = 0
        mapVari = 0
        mapLagv = 0
!
        num_tot = 0
        num_gv = 0
        num_mk = 0
        num_vari = 0
!
        do i_face = 1, hhoCell%nbfaces
            do i_dof = 1, mk_fbs
                num_tot = num_tot+1
                num_mk = num_mk+1
                mapMeca(num_mk) = num_tot
            end do
            do i_dof = 1, gv_fbs
                num_tot = num_tot+1
                num_gv = num_gv+1
                mapVari(num_gv) = num_tot
            end do
        end do
!
        do i_dof = 1, mk_cbs
            num_tot = num_tot+1
            num_mk = num_mk+1
            mapMeca(num_mk) = num_tot
        end do
        do i_dof = 1, gv_cbs
            num_tot = num_tot+1
            num_gv = num_gv+1
            mapVari(num_gv) = num_tot
        end do
        do i_dof = 1, gv_cbs
            num_tot = num_tot+1
            mapLagv(i_dof) = num_tot
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoAssGVRhs(hhoCell, hhoDataMk, hhoDataGv, mapMeca, mapVari, mapLagv, &
                           rhs_meca, rhs_vari, rhs_lagv, rhs)
!
        implicit none
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoDataMk, hhoDataGv
        integer(kind=8), intent(in) :: mapMeca(MSIZE_TDOFS_VEC)
        integer(kind=8), intent(in) :: mapVari(MSIZE_TDOFS_SCAL)
        integer(kind=8), intent(in) :: mapLagv(MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: rhs_meca(MSIZE_TDOFS_VEC)
        real(kind=8), intent(in) :: rhs_vari(MSIZE_TDOFS_SCAL)
        real(kind=8), intent(in) :: rhs_lagv(MSIZE_CELL_SCAL)
        real(kind=8), intent(out) :: rhs(MSIZE_TDOFS_MIX)
!--------------------------------------------------------------------------------------------------
! Assembly RHS for GRAD_VARI
!
!--------------------------------------------------------------------------------------------------
        integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs, gv_cbs, gv_fbs, gv_total_dofs
        integer(kind=8) :: i_dof
!
        call hhoMecaDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs)
        call hhoTherDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs)
!
        rhs = 0.d0
!
        do i_dof = 1, mk_total_dofs
            rhs(mapMeca(i_dof)) = rhs_meca(i_dof)
        end do
        do i_dof = 1, gv_total_dofs
            rhs(mapVari(i_dof)) = rhs_vari(i_dof)
        end do
        do i_dof = 1, gv_cbs
            rhs(mapLagv(i_dof)) = rhs_lagv(i_dof)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoAssGVLhs(hhoCell, hhoDataMk, hhoDataGv, mapMeca, mapVari, mapLagv, &
                           lhs_mm, lhs_mv, lhs_ml, lhs_vm, lhs_vv, &
                           lhs_vl, lhs_lm, lhs_lv, lhs_ll, lhs)
!
        implicit none
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoDataMk, hhoDataGv
        integer(kind=8), intent(in) :: mapMeca(MSIZE_TDOFS_VEC)
        integer(kind=8), intent(in) :: mapVari(MSIZE_TDOFS_SCAL)
        integer(kind=8), intent(in) :: mapLagv(MSIZE_CELL_SCAL)
        type(HHO_matrix), intent(in) :: lhs_mv, lhs_ml, lhs_mm, lhs_ll, lhs_vm
        type(HHO_matrix), intent(in) :: lhs_vv, lhs_vl, lhs_lm, lhs_lv
        type(HHO_matrix), intent(out) :: lhs
!--------------------------------------------------------------------------------------------------
! Assembly LHS for GRAD_VARI
!
!--------------------------------------------------------------------------------------------------
        integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs, gv_cbs, gv_fbs, gv_total_dofs
        integer(kind=8) :: i_row, i_col, total_dofs
!
        call hhoMecaDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs)
        call hhoTherDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs)
!
        total_dofs = mk_total_dofs+gv_total_dofs+gv_cbs
        call lhs%initialize(total_dofs, total_dofs, 0.d0)
!
! --- Bloc meca
        do i_row = 1, mk_total_dofs
            do i_col = 1, mk_total_dofs
                lhs%m(mapMeca(i_row), mapMeca(i_col)) = lhs_mm%m(i_row, i_col)
            end do
            do i_col = 1, gv_total_dofs
                lhs%m(mapMeca(i_row), mapVari(i_col)) = lhs_mv%m(i_row, i_col)
            end do
            do i_col = 1, gv_cbs
                lhs%m(mapMeca(i_row), mapLagv(i_col)) = lhs_ml%m(i_row, i_col)
            end do
        end do
!
! --- Bloc vari
        do i_row = 1, gv_total_dofs
            do i_col = 1, mk_total_dofs
                lhs%m(mapVari(i_row), mapMeca(i_col)) = lhs_vm%m(i_row, i_col)
            end do
            do i_col = 1, gv_total_dofs
                lhs%m(mapVari(i_row), mapVari(i_col)) = lhs_vv%m(i_row, i_col)
            end do
            do i_col = 1, gv_cbs
                lhs%m(mapVari(i_row), mapLagv(i_col)) = lhs_vl%m(i_row, i_col)
            end do
        end do
!
! --- Bloc lagr
        do i_row = 1, gv_cbs
            do i_col = 1, mk_total_dofs
                lhs%m(mapLagv(i_row), mapMeca(i_col)) = lhs_lm%m(i_row, i_col)
            end do
            do i_col = 1, gv_total_dofs
                lhs%m(mapLagv(i_row), mapVari(i_col)) = lhs_lv%m(i_row, i_col)
            end do
            do i_col = 1, gv_cbs
                lhs%m(mapLagv(i_row), mapLagv(i_col)) = lhs_ll%m(i_row, i_col)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_gv(this, hhoCell, hhoDataMk, hhoDataGv, hhoCS)
!
        implicit none
!
        class(HHO_GV_State), intent(inout) :: this
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoDataMk, hhoDataGv
        type(HHO_Compor_State), intent(in) :: hhoCS
!
! --------------------------------------------------------------------------------------------------
!
!  initialize HHO_GV_STATE
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs
        integer(kind=8) :: gv_cbs, gv_fbs, gv_total_dofs, total_dofs
        real(kind=8) :: tmp_prev(MSIZE_TDOFS_MIX), tmp_incr(MSIZE_TDOFS_MIX)
        aster_logical :: forc_noda, pilo
!
        forc_noda = hhoCS%option == "FORC_NODA"
        pilo = hhoCS%option(1:4) == "PILO"
        if (hhoCS%option .ne. "RIGI_MECA") then
            if (hhoCS%typmod(2) .ne. "GRADVARI") then
                ASSERT(ASTER_FALSE)
            end if
!
            call hhoMecaDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs)
            call hhoTherDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs)
            total_dofs = mk_total_dofs+gv_total_dofs+gv_cbs
!
            if (forc_noda) then
                call readVector('PDEPLAR', total_dofs, tmp_prev)
            else
                call readVector('PDEPLMR', total_dofs, tmp_prev)
                if (pilo) then
                    call readVector('PDDEPLR', total_dofs, tmp_incr)
                else
                    call readVector('PDEPLPR', total_dofs, tmp_incr)
                end if
            end if
!
            call hhoExtrField(hhoCell, hhoDataMk, hhoDataGv, &
                              tmp_prev, vari=this%vari_prev, lagv=this%lagv_prev)
            if (.not. forc_noda) then
                call hhoExtrField(hhoCell, hhoDataMk, hhoDataGv, &
                                  tmp_incr, vari=this%vari_incr, lagv=this%lagv_incr)
            end if
        end if
!
! --- compute in T+
!
        call dcopy_1(gv_total_dofs, this%vari_prev, this%vari_curr)
        call dcopy_1(gv_cbs, this%lagv_prev, this%lagv_curr)
!
        if (.not. forc_noda) then
            call daxpy_1(gv_total_dofs, 1.d0, this%vari_incr, this%vari_curr)
            call daxpy_1(gv_cbs, 1.d0, this%lagv_incr, this%lagv_curr)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine free_gv(this)
!
        implicit none
!
        class(HHO_GV_State), intent(inout) :: this
!
! --------------------------------------------------------------------------------------------------
!
!  clean HHO_GV_STATE
! --------------------------------------------------------------------------------------------------
!
        call this%grad%free()
        call this%stab%free()
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalcOpGv(hhoCell, hhoDataMk, hhoDataGv, l_largestrains, hhoMecaState, hhoGvState)
!
        implicit none
!
        type(HHO_Cell), intent(inout) :: hhoCell
        type(HHO_Data), intent(in) :: hhoDataMk, hhoDataGv
        aster_logical, intent(in) :: l_largestrains
        type(HHO_Meca_State), intent(inout) :: hhoMecaState
        type(HHO_GV_State), intent(inout) :: hhoGvState
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Mechanic
!
! Compute operators for mechanic
!
! --------------------------------------------------------------------------------------------------
!
! In  hhoCell         : hho Cell
! In hhoData          : information about the HHO formulation
! In l_largestrains   : large strains ?
! Out gradfull        : full gradient for mechanics
! Out stab            : stabilization for mechanics
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: gv_cbs, gv_fbs, gv_total_dofs, gv_gbs, iret, iad
        aster_logical :: l_stab
!
        call tecach('NNN', 'PCHHOST', 'L', iret, iad=iad)
        l_stab = iret == 0
!
        if (ASTER_TRUE) then
!
            if (l_stab) then
                call hhoReloadPreCalcMeca(hhoCell, hhoDataMk, l_largestrains, hhoMecaState%grad, &
                                          hhoMecaState%stab)
            else
                call hhoReloadPreCalcMeca(hhoCell, hhoDataMk, l_largestrains, hhoMecaState%grad)
            end if
!
            if (hhoDataMk%grad_degree() == hhoDataGv%grad_degree()) then
                call hhoTherNLDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs, gv_gbs)

                call hhoGVState%grad%initialize(gv_gbs, gv_total_dofs, 0.d0)
                call hhoGVState%grad%read('PCHHOGT', ASTER_FALSE)
                if (l_stab) then
                    call hhoGVState%stab%initialize(gv_total_dofs, gv_total_dofs, 0.d0)
                    call hhoGVState%stab%read('PCHHOST', ASTER_TRUE)
                end if
            else
                call hhoCalcOpTher(hhoCell, hhoDataGv, hhoGVState%grad, hhoGVState%stab)
            end if
!
        else
            call hhoCalcOpMeca(hhoCell, hhoDataMk, l_largestrains, hhoMecaState%grad, &
                               hhoMecaState%stab)
            call hhoCalcOpTher(hhoCell, hhoDataGv, hhoGVState%grad, hhoGVState%stab)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsLargeMV(hhoCell, dPK1_dv, weight, BSCEval, gv_cbs, &
                                    mk_gbs, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: dPK1_dv(3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gv_cbs, mk_gbs
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (dPK1_dv:cphi, gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In dPK1_dv      : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_dPK1_dv_ij
        integer(kind=8) :: i, j, k, row, gbs_cmp, offset, axis_dim
! --------------------------------------------------------------------------------------------------
!
        offset = AT%ncols-gv_cbs+1
        axis_dim = 0
        if (hhoCell%l_axis) then
            axis_dim = 1
        end if
!
! -------- Compute scalar_product of (dPK1_dv:cphi, gphi)_T
        gbs_cmp = mk_gbs/(hhoCell%ndim*hhoCell%ndim+axis_dim)
!
        row = 1
        do i = 1, hhoCell%ndim
            do j = 1, hhoCell%ndim
                qp_dPK1_dv_ij = weight*dPK1_dv(i, j)
                do k = 1, gbs_cmp
                    call daxpy_1(gv_cbs, qp_dPK1_dv_ij*BSCEval(k), BSCEval, AT%m(row, offset:))
                    row = row+1
                end do
            end do
        end do
!
        if (hhoCell%l_axis) then
            qp_dPK1_dv_ij = weight*dPK1_dv(3, 3)
            do k = 1, gbs_cmp
                call daxpy_1(gv_cbs, qp_dPK1_dv_ij*BSCEval(k), BSCEval, AT%m(row, offset:))
                row = row+1
            end do
        end if
!
        ASSERT(mk_gbs == row-1)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsLargeVM(hhoCell, dsv_dF, weight, BSCEval, gv_cbs, &
                                    mk_gbs, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: dsv_dF(3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gv_cbs, mk_gbs
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (cphi, dsv_dF:gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In dsv_dF       : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_dsv_dF_ij
        integer(kind=8) :: i, j, k, col, gbs_cmp, offset, axis_dim
! --------------------------------------------------------------------------------------------------
!
! --------- Eval (dsv_dF : sgphi)_T
        axis_dim = 0
        if (hhoCell%l_axis) then
            axis_dim = 1
        end if
        gbs_cmp = mk_gbs/(hhoCell%ndim*hhoCell%ndim+axis_dim)
        offset = AT%nrows-gv_cbs+1
!
! -------- Compute scalar_product of (cphi, dsv_dF:gphi)_T
        col = 1
        do i = 1, hhoCell%ndim
            do j = 1, hhoCell%ndim
                qp_dsv_dF_ij = weight*dsv_dF(i, j)
                do k = 1, gbs_cmp
                    call daxpy_1(gv_cbs, qp_dsv_dF_ij*BSCEval(k), BSCEval, AT%m(offset:, col))
                    col = col+1
                end do
            end do
        end do
!
        if (hhoCell%l_axis) then
            qp_dsv_dF_ij = weight*dsv_dF(3, 3)
            do k = 1, gbs_cmp
                call daxpy_1(gv_cbs, qp_dsv_dF_ij*BSCEval(k), BSCEval, AT%m(offset:, col))
                col = col+1
            end do
        end if
!
        ASSERT(mk_gbs == col-1)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsLargeML(hhoCell, dPK1_dl, weight, BSCEval, gv_cbs, &
                                    mk_gbs, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: dPK1_dl(3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gv_cbs, mk_gbs
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (dPK1_dl:cphi, gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In dPK1_dl      : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_dPK1_dl_ij
        integer(kind=8) :: i, j, k, row, gbs_cmp, axis_dim
! --------------------------------------------------------------------------------------------------
!
! -------- Compute scalar_product of (dPK1_dl:cphi, gphi)_T
        axis_dim = 0
        if (hhoCell%l_axis) then
            axis_dim = 1
        end if
        gbs_cmp = mk_gbs/(hhoCell%ndim*hhoCell%ndim+axis_dim)
!
        row = 1
        do i = 1, hhoCell%ndim
            do j = 1, hhoCell%ndim
                qp_dPK1_dl_ij = weight*dPK1_dl(i, j)
                do k = 1, gbs_cmp
                    call daxpy_1(gv_cbs, qp_dPK1_dl_ij*BSCEval(k), BSCEval, AT%m(row, :))
                    row = row+1
                end do
            end do
        end do
!
        if (hhoCell%l_axis) then
            qp_dPK1_dl_ij = weight*dPK1_dl(3, 3)
            do k = 1, gbs_cmp
                call daxpy_1(gv_cbs, qp_dPK1_dl_ij*BSCEval(k), BSCEval, AT%m(row, :))
                row = row+1
            end do
        end if
!
        ASSERT(mk_gbs == row-1)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsLargeLM(hhoCell, dsl_dF, weight, BSCEval, gv_cbs, &
                                    mk_gbs, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: dsl_dF(3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gv_cbs, mk_gbs
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (cphi, dsl_dF:gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In dsl_dF       : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_dsl_dF_ij
        integer(kind=8) :: i, j, k, col, gbs_cmp, axis_dim
! --------------------------------------------------------------------------------------------------
!
! --------- Eval (dsl_dF : sgphi)_T
        axis_dim = 0
        if (hhoCell%l_axis) then
            axis_dim = 1
        end if
        gbs_cmp = mk_gbs/(hhoCell%ndim*hhoCell%ndim+axis_dim)
!
! -------- Compute scalar_product of (C_sgphi(j), gphi(j))_T
        col = 1
        do i = 1, hhoCell%ndim
            do j = 1, hhoCell%ndim
                qp_dsl_dF_ij = weight*dsl_dF(i, j)
                do k = 1, gbs_cmp
                    call daxpy_1(gv_cbs, qp_dsl_dF_ij*BSCEval(k), BSCEval, AT%m(:, col))
                    col = col+1
                end do
            end do
        end do
!
        if (hhoCell%l_axis) then
            qp_dsl_dF_ij = weight*dsl_dF(3, 3)
            do k = 1, gbs_cmp
                call daxpy_1(gv_cbs, qp_dsl_dF_ij*BSCEval(k), BSCEval, AT%m(:, col))
                col = col+1
            end do
        end if
!
        ASSERT(mk_gbs == col-1)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsSmallMV(hhoCell, dSig_dv, weight, BSCEval, gv_cbs, &
                                    mk_gbs_cmp, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: dSig_dv(6)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gv_cbs, mk_gbs_cmp
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (dSig_dv:cphi, gsphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In dSig_dv      : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Acphi(6)
        integer(kind=8) :: i, k, col, deca
! --------------------------------------------------------------------------------------------------
!
! -------- Compute scalar_product of (qp_Acphi, sgphi)_T
! -------- (RAPPEL: the composents of the gradient are saved by G11, G22, G33, G12, G13, G23)
        col = AT%ncols-gv_cbs+1
        do k = 1, gv_cbs
! --------- Eval (dSig_dv : scphi)_T
            qp_Acphi = weight*dSig_dv*BSCEval(k)
            deca = 1
            do i = 1, hhoCell%ndim
                call daxpy_1(mk_gbs_cmp, qp_Acphi(i), BSCEval, AT%m(deca:, col))
                deca = deca+mk_gbs_cmp
            end do
!
! ---- non-diagonal terms
            select case (hhoCell%ndim)
            case (3)
                do i = 1, 3
                    call daxpy_1(mk_gbs_cmp, qp_Acphi(3+i), BSCEval, AT%m(deca:, col))
                    deca = deca+mk_gbs_cmp
                end do
            case (2)
                call daxpy_1(mk_gbs_cmp, qp_Acphi(4), BSCEval, AT%m(deca:, col))
                deca = deca+mk_gbs_cmp
                if (hhoCell%l_axis) then
! -------- Compute scalar_product of (qp_Acphi, Pi^k_T)_T
                    call daxpy_1(mk_gbs_cmp, qp_Acphi(3), BSCEval, AT%m(deca:, col))
                    deca = deca+mk_gbs_cmp
                end if
            case default
                ASSERT(ASTER_FALSE)
            end select
!
            col = col+1
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsSmallVM(hhoCell, dsv_dEps, weight, BSCEval, gv_cbs, &
                                    mk_gbs_sym, mk_gbs_cmp, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: dsv_dEps(6)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gv_cbs, mk_gbs_sym, mk_gbs_cmp
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (cphi, dsv_dEps:sgphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In dsv_dEps     : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In mk_gbs_cmp   : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Agphi(MSIZE_CELL_MAT), qp_dsv_dEps(6)
        integer(kind=8) :: i, offset, deca
! --------------------------------------------------------------------------------------------------
!
! -------- Compute (cphi, dsv_dEps:sgphi)
! -------- (RAPPEL: the composents of the gradient are saved by G11, G22, G33, G12, G13, G23)
        qp_dsv_dEps = weight*dsv_dEps
        deca = 1
        do i = 1, hhoCell%ndim
            call daxpy_1(mk_gbs_cmp, qp_dsv_dEps(i), BSCEval, qp_Agphi(deca))
            deca = deca+mk_gbs_cmp
        end do
!
! ---- non-diagonal terms
        select case (hhoCell%ndim)
        case (3)
            do i = 1, 3
                call daxpy_1(mk_gbs_cmp, qp_dsv_dEps(3+i), BSCEval, qp_Agphi(deca))
                deca = deca+mk_gbs_cmp
            end do
        case (2)
            call daxpy_1(mk_gbs_cmp, qp_dsv_dEps(4), BSCEval, qp_Agphi(deca))
            deca = deca+mk_gbs_cmp
            if (hhoCell%l_axis) then
                call daxpy_1(mk_gbs_cmp, qp_dsv_dEps(3), BSCEval, qp_Agphi(deca))
                deca = deca+mk_gbs_cmp
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
        ASSERT(deca-1 == mk_gbs_sym)
!
        offset = AT%nrows-gv_cbs+1
!
! -------- Compute scalar_product of (qp_Aphi, sgphi)_T
! -------- (RAPPEL: the composents of the gradient are saved by G11, G22, G33, G12, G13, G23)
        do i = 1, mk_gbs_sym
            call daxpy_1(gv_cbs, qp_Agphi(i), BSCEval, AT%m(offset:, i))
        end do
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsSmallML(hhoCell, dSig_dl, weight, BSCEval, gv_cbs, &
                                    mk_gbs_cmp, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: dSig_dl(6)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gv_cbs, mk_gbs_cmp
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (dSig_dl:cphi, gsphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In dSig_dl      : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Acphi(6)
        integer(kind=8) :: i, k, deca
! --------------------------------------------------------------------------------------------------
!
! -------- Compute scalar_product of (dSig_dl:cphi, gsphi)_T
! -------- (RAPPEL: the composents of the gradient are saved by G11, G22, G33, G12, G13, G23)
        do k = 1, gv_cbs
! --------- Eval (dSig_dl : scphi)_T
            qp_Acphi = weight*dSig_dl*BSCEval(k)
            deca = 1
            do i = 1, hhoCell%ndim
                call daxpy_1(mk_gbs_cmp, qp_Acphi(i), BSCEval, AT%m(deca:, k))
                deca = deca+mk_gbs_cmp
            end do
!
            ! ---- non-diagonal terms
            select case (hhoCell%ndim)
            case (3)
                do i = 1, 3
                    call daxpy_1(mk_gbs_cmp, qp_Acphi(3+i), BSCEval, AT%m(deca:, k))
                    deca = deca+mk_gbs_cmp
                end do
            case (2)
                call daxpy_1(mk_gbs_cmp, qp_Acphi(4), BSCEval, AT%m(deca:, k))
                deca = deca+mk_gbs_cmp
                if (hhoCell%l_axis) then
                    ! -------- Compute scalar_product of (qp_Acphi, Pi^k_T(ur/r))_T
                    call daxpy_1(mk_gbs_cmp, qp_Acphi(3), BSCEval, AT%m(deca:, k))
                    deca = deca+mk_gbs_cmp
                end if
            case default
                ASSERT(ASTER_FALSE)
            end select
        end do
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsSmallLM(hhoCell, dsl_dEps, weight, BSCEval, gv_cbs, &
                                    mk_gbs_sym, mk_gbs_cmp, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: dsl_dEps(6)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gv_cbs, mk_gbs_sym, mk_gbs_cmp
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (cphi, dsl_dEps:gsphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In dsl_dEps     : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Agphi(MSIZE_CELL_MAT), qp_dsl_dEps(6)
        integer(kind=8) :: i, deca
! --------------------------------------------------------------------------------------------------
!
! -------- Compute (dsl_dEps : gsphi)
! -------- (RAPPEL: the composents of the gradient are saved by G11, G22, G33, G12, G13, G23)
        qp_dsl_dEps = weight*dsl_dEps
        deca = 1
        do i = 1, hhoCell%ndim
            call daxpy_1(mk_gbs_cmp, qp_dsl_dEps(i), BSCEval, qp_Agphi(deca))
            deca = deca+mk_gbs_cmp
        end do
!
! ---- non-diagonal terms
        select case (hhoCell%ndim)
        case (3)
            do i = 1, 3
                call daxpy_1(mk_gbs_cmp, qp_dsl_dEps(3+i), BSCEval, qp_Agphi(deca))
                deca = deca+mk_gbs_cmp
            end do
        case (2)
            call daxpy_1(mk_gbs_cmp, qp_dsl_dEps(4), BSCEval, qp_Agphi(deca))
            deca = deca+mk_gbs_cmp
            if (hhoCell%l_axis) then
                call daxpy_1(mk_gbs_cmp, qp_dsl_dEps(3), BSCEval, qp_Agphi(deca))
                deca = deca+mk_gbs_cmp
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
        ASSERT(deca-1 == mk_gbs_sym)
!
! -------- Compute scalar_product of (qp_Acphi, sgphi)_T
! -------- (RAPPEL: the composents of the gradient are saved by G11, G22, G33, G12, G13, G23)
        do i = 1, mk_gbs_sym
            call daxpy_1(gv_cbs, qp_Agphi(i), BSCEval, AT%m(:, i))
        end do
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine forc_noda_stress(hhoCS, ndim, kpg, F_curr, PK1_curr, &
                                Cauchy_curr, sig_vari, sig_lagv, sig_gv)
!
        implicit none
!
        type(HHO_Compor_State), intent(inout) :: hhoCS
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: F_curr(3, 3)
        real(kind=8), intent(out) :: PK1_curr(3, 3), Cauchy_curr(6)
        real(kind=8), intent(out) :: sig_vari, sig_lagv, sig_gv(3)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Get stress for FORC_NODA
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: neu, neg, ntot
        real(kind=8) :: sigPrev(11)
! --------------------------------------------------------------------------------------------------
!
        neu = 2*ndim
        neg = 2+ndim
        ntot = neu+neg
!
        sigPrev(1:ntot) = hhoCS%sig_prev((kpg-1)*ntot+1:kpg*ntot)
        sig_vari = sigPrev(neu+1)
        sig_lagv = sigPrev(neu+2)
        sig_gv = 0.d0
        sig_gv(1:ndim) = sigPrev(neu+2+1:ntot)
!
        call tranfoMatToSym(ndim, sigPrev(1:neu), Cauchy_curr)
!
        PK1_curr = 0.d0
        if (hhoCS%l_largestrain) then
            call sigtopk1(ndim, Cauchy_curr, F_curr, PK1_curr)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function hhoCalcStabCoeffGV(fami)
!
        implicit none
!
        character(len=8) :: fami
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  GRAD_VARI - Evaluate stabilzation coefficient
! --------------------------------------------------------------------------------------------------
!
! --- Local variables
!
        integer(kind=8) :: jvMaterc, jvMaterCode
        integer(kind=8) :: iok(1)
        real(kind=8) :: vale(1)
!
        call jevech('PMATERC', 'L', jvMaterc)
        jvMaterCode = zi(jvMaterc-1+1)
!
! ----- C_GRAD_VARI is constant on the cell
        call rcvalb(fami, 1, 1, '+', jvMaterCode, &
                    ' ', 'NON_LOCAL', &
                    0, ' ', [0.d0], &
                    1, ['C_GRAD_VARI'], vale, iok, 1)
        hhoCalcStabCoeffGV = 10.d0*vale(1)
!
        if (hhoCalcStabCoeffGV <= 0.0) then
            hhoCalcStabCoeffGV = 1.d0
        end if
        ASSERT(hhoCalcStabCoeffGV > 0.d0)
!
    end function
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoDataGVInit(hhoDataGv)
!
        implicit none
!
        type(HHO_Data), intent(out) :: hhoDataGv

!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Initi data for GV
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: grad_deg, face_deg, cell_deg
!
        ASSERT(lteatt('TYPMOD2', 'HHO') .or. lteatt('TYPMOD2', 'HHO_GRAD'))
!
        if (lteatt('FORMULATION', 'HHO_LINE')) then
            face_deg = 1
            cell_deg = 1
            grad_deg = 1
        elseif (lteatt('FORMULATION', 'HHO_VLINE')) then
            face_deg = 0
            cell_deg = 0
            grad_deg = 0
        elseif (lteatt('FORMULATION', 'HHO_QUAD')) then
            face_deg = 2
            cell_deg = 2
            grad_deg = 2
        elseif (lteatt('FORMULATION', 'HHO_VQUAD')) then
            face_deg = 1
            cell_deg = 1
            grad_deg = 1
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Init
!
        call hhoDataGv%initialize(face_deg, cell_deg, grad_deg, 0.0, ASTER_FALSE, ASTER_FALSE)
!
        if (hhoDataGv%debug()) then
            call hhoDataGv%print()
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoExtrField(hhoCell, hhoDataMk, hhoDataGv, &
                            field, depl, vari, lagv)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoDataMk, hhoDataGv
        real(kind=8), intent(in) :: field(MSIZE_TDOFS_MIX)
        real(kind=8), intent(out), optional :: depl(MSIZE_TDOFS_VEC)
        real(kind=8), intent(out), optional :: vari(MSIZE_TDOFS_SCAL)
        real(kind=8), intent(out), optional :: lagv(MSIZE_CELL_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Extract each component
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs
        integer(kind=8) :: gv_cbs, gv_fbs, gv_total_dofs, total_dofs
        integer(kind=8) :: num_tot, num_gv, iFace, iDof, num_mk
! --------------------------------------------------------------------------------------------------
!
        call hhoMecaDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs)
        call hhoTherDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs)
        total_dofs = mk_total_dofs+gv_total_dofs+gv_cbs
!
        num_tot = 0
        num_gv = 0
        num_mk = 0
        do iFace = 1, hhoCell%nbfaces
            if (present(depl)) then
                do iDof = 1, mk_fbs
                    num_tot = num_tot+1
                    num_mk = num_mk+1
                    depl(num_mk) = field(num_tot)
                end do
            else
                num_tot = num_tot+mk_fbs
            end if
!
            if (present(vari)) then
                do iDof = 1, gv_fbs
                    num_tot = num_tot+1
                    num_gv = num_gv+1
                    vari(num_gv) = field(num_tot)
                end do
            else
                num_tot = num_tot+gv_fbs
            end if
        end do
!
        if (present(depl)) then
            do iDof = 1, mk_cbs
                num_tot = num_tot+1
                num_mk = num_mk+1
                depl(num_mk) = field(num_tot)
            end do
        else
            num_tot = num_tot+mk_cbs
        end if
!
        if (present(vari)) then
            do iDof = 1, gv_cbs
                num_tot = num_tot+1
                num_gv = num_gv+1
                vari(num_gv) = field(num_tot)
            end do
        else
            num_tot = num_tot+gv_cbs
        end if
!
        if (present(lagv)) then
            do iDof = 1, gv_cbs
                num_tot = num_tot+1
                lagv(iDof) = field(num_tot)
            end do
        else
            num_tot = num_tot+gv_cbs
        end if
!
        ASSERT(num_tot == total_dofs)
!
    end subroutine
!
end module
