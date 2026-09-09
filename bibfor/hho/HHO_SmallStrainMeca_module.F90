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

module HHO_SmallStrainMeca_module
!
    use Behaviour_module
    use Behaviour_type
    use FE_algebra_module
    use HHO_algebra_module
    use HHO_basis_module
    use HHO_compor_module
    use HHO_eval_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use HHO_utils_module
    use HHO_L2proj_module
!
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/dmatmc.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/nbsigm.h"
#include "asterfort/nmcomp.h"
#include "blas/daxpy.h"
#include "blas/dger.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - mechanics
!
! Module for small deformations with hho
!
! --------------------------------------------------------------------------------------------------
!
    public :: hhoSmallStrainLCMeca, tranfoMatToSym, hhoMatrElasMeca
    public :: hhoComputeRhsSmall, hhoComputeLhsSmall, hhoAssembleLhsSmall
    public :: hhoComputeCgphi, tranfoSymToMat, tranfoTensToSym
    private :: hhoComputeLhsSmallAxis
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoSmallStrainLCMeca(hhoCell, hhoData, hhoQuadCellRigi, hhoCS, gradrec, &
                                    time_prev, time_curr, depl_prev, depl_incr, &
                                    lhs, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellRigi
        type(HHO_Compor_State), intent(inout) :: hhoCS
        type(HHO_matrix), intent(in) :: gradrec
        real(kind=8), intent(in) :: time_prev
        real(kind=8), intent(in) :: time_curr
        real(kind=8), intent(in) :: depl_prev(MSIZE_TDOFS_VEC)
        real(kind=8), intent(in) :: depl_incr(MSIZE_TDOFS_VEC)
        type(HHO_matrix), intent(inout) :: lhs
        real(kind=8), intent(inout) :: rhs(MSIZE_TDOFS_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the local contribution for mechanics with small deformations
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   InOut hhoCS        : hho compor state
!   In gradrec      : local gradient reconstruction
!   In time_prev    : previous time T-
!   In time_curr    : current time T+
!   In depl_prev    : displacement at T-
!   In depl_incr    : increment of displacement between T- and T+
!   Out lhs         : local contribution (lhs)
!   Out rhs         : local contribution (rhs)
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: ksp = 1
        type(HHO_basis_cell) :: hhoBasisCell
        real(kind=8) :: E_prev_coeff(MSIZE_CELL_MAT), E_incr_coeff(MSIZE_CELL_MAT)
        real(kind=8) :: dsidep(6, 6), E_prev(6), E_incr(6), Cauchy_curr(6), Cauchy_prev(6)
        real(kind=8) :: coorpg(3), weight
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL), bT(MSIZE_CELL_MAT)
        type(HHO_matrix) :: AT
        integer(kind=8) :: cbs, fbs, faces_dofs, total_dofs, gbs, kpg, gbs_cmp, gbs_sym, cbs_cmp
        integer(kind=8) :: cod(MSIZE_QP_CELL), gbs_axis
        aster_logical :: l_lhs, l_rhs
! --------------------------------------------------------------------------------------------------
!
        cod = 0
! ------ number of dofs
        call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                           gbs, gbs_sym, gbs_axis)
        faces_dofs = total_dofs-cbs
        gbs_cmp = (gbs-gbs_axis)/(hhoCell%ndim*hhoCell%ndim)
        cbs_cmp = cbs/hhoCell%ndim
!
        bT = 0.d0
        dsidep = 0.d0
        E_prev_coeff = 0.d0
        E_incr_coeff = 0.d0
        Cauchy_curr = 0.d0
        l_lhs = L_MATR(hhoCS%option)
        l_rhs = L_VECT(hhoCS%option)

        if (l_lhs) then
            call AT%initialize(gbs_sym, gbs_sym, 0.d0)
        end if

! ----- Prepare external state variables (geometry)
        call behaviourPrepESVAGeomHHO(hhoCell, hhoQuadCellRigi, hhoCS%BEHInteg)

! ----- init basis
        call hhoBasisCell%initialize(hhoCell)

! ----- compute E_prev = gradrec_sym * depl_prev
        call gradrec%dot(depl_prev, E_prev_coeff)

! ----- compute E_incr = gradrec_sym * depl_incr
        call gradrec%dot(depl_incr, E_incr_coeff)
!
! ----- Loop on quadrature point
!
        do kpg = 1, hhoQuadCellRigi%nbQuadPoints
            coorpg(1:3) = hhoQuadCellRigi%points(1:3, kpg)
            weight = hhoQuadCellRigi%weights(kpg)
! --------- Eval basis function at the quadrature point
            call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
! --------- Eval deformations
            E_prev = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_prev_coeff)
!
            E_incr = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_incr_coeff)
!
! --------- tranform sigm in symmetric form
            call tranfoMatToSym(hhoCell%ndim, &
                                hhoCS%sig_prev((kpg-1)*hhoCS%nbsigm+1:kpg*hhoCS%nbsigm), &
                                Cauchy_prev)

! --------- Set main parameters for behaviour (on point)
            call behaviourSetParaPoin(kpg, ksp, hhoCS%BEHInteg)

! --------- Integrate
            call nmcomp(hhoCS%BEHInteg, &
                        hhoCell%ndim, hhoCS%option, hhoCS%typmod, &
                        time_prev, time_curr, &
                        hhoCS%compor, hhoCS%carcri, hhoCS%multComp, &
                        6, E_prev, E_incr, &
                        6, Cauchy_prev, &
                        hhoCS%vari_prev((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                        Cauchy_curr, &
                        hhoCS%vari_curr((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                        36, dsidep, &
                        cod(kpg))
!
            if (cod(kpg) .eq. 1) then
                goto 999
            end if
!
! --------- For new prediction and nmisot.F90
            if (L_PRED(hhoCS%option)) then
                Cauchy_curr = 0.d0
            end if
!
            if (L_SIGM(hhoCS%option)) then
! -------- tranform Cauchy_curr in symmetric form
                call tranfoSymToMat(hhoCell%ndim, Cauchy_curr, &
                                    hhoCS%sig_curr((kpg-1)*hhoCS%nbsigm+1:kpg*hhoCS%nbsigm))
            end if
!
            if (l_rhs) then
                call hhoComputeRhsSmall(hhoCell, Cauchy_curr, weight, BSCEval, gbs_cmp, bT)
            end if
!
            if (l_lhs) then
                call hhoComputeLhsSmall(hhoCell, dsidep, hhoCS%matsym, weight, BSCEval, &
                                        gbs_sym, gbs_cmp, AT)
            end if
        end do
!
! ----- compute rhs += Gradrec**T * bT
        if (l_rhs) then
            call hho_dgemv_T(1.d0, gradrec, bT, 1.d0, rhs)
        end if
!
        if (l_lhs) then
            call hhoAssembleLhsSmall(hhoCS, gradrec, AT, lhs)
        end if
!
! print*, "AT", hhoNorm2Mat(AT(1:gbs_sym,1:gbs_sym))
! print*, "bT", norm2(bT)
! print*, "KT", hhoNorm2Mat(lhs(1:total_dofs,1:total_dofs))
! print*, "fT", norm2(rhs)
!
999     continue
!
! ---- Return code summary
!
        call codere(cod, hhoQuadCellRigi%nbQuadPoints, hhoCS%codret)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatrElasMeca(hhoCell, hhoData, hhoQuadCellRigi, hhoCS, gradrec, &
                               time_curr, lhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellRigi
        type(HHO_Compor_State), intent(inout) :: hhoCS
        type(HHO_matrix), intent(in) :: gradrec
        real(kind=8), intent(in) :: time_curr
        type(HHO_matrix), intent(inout) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute matrix for RIGI_MECA
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   In hhoCS        : hho compor state
!   In gradrec      : local gradient reconstruction
!   In time_curr    : current time T+
!   Out lhs         : local contribution (lhs)
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: ksp = 1
        type(HHO_basis_cell) :: hhoBasisCell
        real(kind=8) :: dsidep(6, 6), dsidep3D(6, 6)
        real(kind=8) :: coorpg(3), weight
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
        type(HHO_matrix) :: AT
        integer(kind=8) :: cbs, fbs, total_dofs, faces_dofs, gbs, kpg, gbs_cmp, gbs_sym, nb_sig
        integer(kind=8) :: cbs_cmp, gbs_axis
!
! --------------------------------------------------------------------------------------------------
!
! ----- number of dofs
        call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                           gbs, gbs_sym, gbs_axis)
        faces_dofs = total_dofs-cbs
        gbs_cmp = (gbs-gbs_axis)/(hhoCell%ndim*hhoCell%ndim)
        cbs_cmp = cbs/hhoCell%ndim
!
        dsidep = 0.d0
        nb_sig = nbsigm()
!
        if (hhoCS%option /= "RIGI_MECA") then
            ASSERT(ASTER_FALSE)
        end if

        call AT%initialize(gbs_sym, gbs_sym, 0.d0)
!
! ----- init basis
        call hhoBasisCell%initialize(hhoCell)

! ----- Loop on quadrature points
        do kpg = 1, hhoQuadCellRigi%nbQuadPoints
            coorpg(1:3) = hhoQuadCellRigi%points(1:3, kpg)
            weight = hhoQuadCellRigi%weights(kpg)

! --------- Set main parameters for behaviour (on point)
            call behaviourSetParaPoin(kpg, ksp, hhoCS%BEHInteg)

! --------- Eval basis function at the quadrature point
            call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)

! --------- Compute behaviour
            call dmatmc(hhoCS%BEHInteg%materPara, '+', time_curr, &
                        nb_sig, dsidep)
            call tranfoTensToSym(nb_sig, dsidep, dsidep3D)
!
            call hhoComputeLhsSmall(hhoCell, dsidep3D, ASTER_TRUE, weight, BSCEval, gbs_sym, &
                                    gbs_cmp, AT)
        end do
!
! ----- compute lhs += gradrec**T * AT * gradrec
!
        call hhoAssembleLhsSmall(hhoCS, gradrec, AT, lhs)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeRhsSmall(hhoCell, stress, weight, BSCEval, gbs_cmp, bT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: stress(6)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gbs_cmp
        real(kind=8), intent(inout) :: bT(MSIZE_CELL_MAT)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product bT += (stress, sgphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In stress       : stress tensor (XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ)
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   Out bT          : contribution of bt
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_stress(6)
        integer(kind=8) :: i, deca
! --------------------------------------------------------------------------------------------------
!
        qp_stress = weight*stress
! -------- Compute scalar_product of (stress, sgphi)_T
! -------- (RAPPEL: the composents of the gradient are saved by G11, G22, G33, G12, G13, G23)
        deca = 0
        do i = 1, hhoCell%ndim
            call daxpy_1(gbs_cmp, qp_stress(i), BSCEval, bT(deca+1))
            deca = deca+gbs_cmp
        end do
!
! ---- non-diagonal terms
        select case (hhoCell%ndim)
        case (3)
            do i = 1, 3
                call daxpy_1(gbs_cmp, qp_stress(3+i), BSCEval, bT(deca+1))
                deca = deca+gbs_cmp
            end do
        case (2)
            call daxpy_1(gbs_cmp, qp_stress(4), BSCEval, bT(deca+1))
            deca = deca+gbs_cmp
!
            if (hhoCell%l_axis) then
                call daxpy_1(gbs_cmp, qp_stress(3), BSCEval, bT(deca+1))
                deca = deca+gbs_cmp
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
    subroutine hhoComputeLhsSmall(hhoCell, module_tang, matsym, weight, BSCEval, gbs_sym, &
                                  gbs_cmp, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: module_tang(6, 6)
        aster_logical, intent(in) :: matsym
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gbs_sym
        integer(kind=8), intent(in) :: gbs_cmp
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (gphi, module_tang:gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Cgphi(6, MSIZE_CELL_MAT)
        integer(kind=8) :: icmp, jcol, ib, irow, gbs_axis
! --------------------------------------------------------------------------------------------------
!
! --------- Eval (C : sgphi)_T
        call hhoComputeCgphi(hhoCell, module_tang, BSCEval, gbs_cmp, weight, &
                             qp_Cgphi)
        gbs_axis = 0
        if (hhoCell%l_axis) then
            gbs_axis = gbs_cmp
        end if
!
! -------- Compute scalar_product of (sgphi(irow), C_sgphi(jcol))_T
        do jcol = 1, gbs_sym-gbs_axis
            irow = 1
! ---------- diagonal term
            do icmp = 1, hhoCell%ndim
                do ib = 1, gbs_cmp
                    AT%m(irow, jcol) = AT%m(irow, jcol)+qp_Cgphi(icmp, jcol)*BSCEval(ib)
                    irow = irow+1
                    if (matsym .and. irow > jcol) then
                        go to 100
                    end if
                end do
            end do
!
! --------- non-diagonal terms
            select case (hhoCell%ndim)
            case (3)
                do icmp = 4, 6
                    do ib = 1, gbs_cmp
                        AT%m(irow, jcol) = AT%m(irow, jcol)+qp_Cgphi(icmp, jcol)*BSCEval(ib)
                        irow = irow+1
                        if (matsym .and. irow > jcol) then
                            go to 100
                        end if
                    end do
                end do
            case (2)
                do ib = 1, gbs_cmp
                    AT%m(irow, jcol) = AT%m(irow, jcol)+qp_Cgphi(4, jcol)*BSCEval(ib)
                    irow = irow+1
                    if (matsym .and. irow > jcol) then
                        go to 100
                    end if
                end do
            case default
                ASSERT(ASTER_FALSE)
            end select
!
100         continue
        end do
!
        if (hhoCell%l_axis) then
            call hhoComputeLhsSmallAxis(hhoCell, matsym, qp_Cgphi, BSCEval, gbs_sym, &
                                        gbs_cmp, AT)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsSmallAxis(hhoCell, matsym, qp_Cgphi, BSCEval, gbs_sym, &
                                      gbs_cmp, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        aster_logical, intent(in) :: matsym
        real(kind=8), intent(in) :: qp_Cgphi(6, MSIZE_CELL_MAT)
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gbs_sym
        integer(kind=8), intent(in) :: gbs_cmp
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (cphi/r, module_tang:cphi/r)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: icmp, jcol, ib, beginAxis, endAxis, irow
        blas_int :: b_incx, b_lda, b_n
! --------------------------------------------------------------------------------------------------
!
        ASSERT(hhoCell%ndim == 2)
        ASSERT(hhoCell%l_axis)
!
        beginAxis = gbs_sym-gbs_cmp+1
        endAxis = gbs_sym
!
! -------- Compute scalar_product of (Pi^k_T, module_tang(3,3):Pi^k_T)_T
        b_n = to_blas_int(gbs_cmp)
        b_incx = to_blas_int(1)
        b_lda = to_blas_int(gbs_cmp)
        call dger(b_n, b_n, 1.0, qp_Cgphi(3, beginAxis:endAxis), b_incx, BSCEval, b_incx, &
                  AT%m(beginAxis:endAxis, beginAxis:endAxis), b_lda)
!
! ---------- extra-diagonal term (sgphi(irow), module_tang(icmp,3)*Pi^k_T(jcol))
        do jcol = beginAxis, endAxis
            irow = 1
            do icmp = 1, 2
                do ib = 1, gbs_cmp
                    AT%m(irow, jcol) = AT%m(irow, jcol)+qp_Cgphi(icmp, jcol)*BSCEval(ib)
                    irow = irow+1
                end do
            end do
!
            do ib = 1, gbs_cmp
                AT%m(irow, jcol) = AT%m(irow, jcol)+qp_Cgphi(4, jcol)*BSCEval(ib)
                irow = irow+1
            end do
        end do
!
        if (.not. matsym) then
! ---------- extra-diagonal term (Pi^k_T(irow), module_tang(3, jcol):sgphi(jcol))
            do irow = beginAxis, endAxis
                ib = irow-beginAxis+1
                do jcol = 1, beginAxis-1
                    AT%m(irow, jcol) = AT%m(irow, jcol)+qp_Cgphi(3, jcol)*BSCEval(ib)
                end do
            end do
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoAssembleLhsSmall(hhoCS, gradrec, AT, lhs)
!
        implicit none
!
        type(HHO_Compor_State), intent(in) :: hhoCS
        type(HHO_matrix), intent(in) :: gradrec
        type(HHO_matrix), intent(inout) ::  AT
        type(HHO_matrix), intent(inout) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics - assemble LHS
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_matrix) :: TMP
        integer(kind=8) :: gbs_sym, total_dofs
! --------------------------------------------------------------------------------------------------
!
!
        gbs_sym = gradrec%nrows
        total_dofs = gradrec%ncols
!
! ----- compute lhs += gradrec**T * AT * gradrec
!
! ----- Copy symetric part of AT
        if (hhoCS%matsym) then
            call AT%copySymU()
        end if
        call TMP%initialize(gbs_sym, total_dofs, 0.d0)
! ----- step1: TMP = AT * gradrec
        call hho_dgemm_NN(1.d0, AT, gradrec, 0.d0, TMP)
!
! ----- step2: lhs += gradrec**T * TMP
        call hho_dgemm_TN(1.d0, gradrec, TMP, 1.d0, lhs)
!
        call TMP%free()
        call AT%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeCgphi(hhoCell, module_tang, BSCEval, gbs_cmp, weight, &
                               Cgphi)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        integer(kind=8), intent(in) :: gbs_cmp
        real(kind=8), intent(in) :: module_tang(6, 6)
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(out) :: Cgphi(6, MSIZE_CELL_MAT)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product qp_weight * (module_tang, sgphi)_T
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto_plastic moduli
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp     : size of BSCEval
!   In weight       : quadrature weight
!   Out Agphi       : matrix of scalar product
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, col, k
        real(kind=8) :: qp_C(6, 6)
! --------------------------------------------------------------------------------------------------
!
        Cgphi = 0.d0
        qp_C = weight*module_tang
        col = 1
!
        select case (hhoCell%ndim)
        case (3)
            do i = 1, 6
                do k = 1, gbs_cmp
                    call daxpy_1(6, BSCEval(k), qp_C(1:6, i), Cgphi(1:6, col))
                    col = col+1
                end do
            end do
        case (2)
! ---------- diagonal terms
            do i = 1, 2
                do k = 1, gbs_cmp
                    Cgphi(1:4, col) = qp_C(1:4, i)*BSCEval(k)
                    col = col+1
                end do
            end do
! ---- non-diagonal terms
            do k = 1, gbs_cmp
                Cgphi(1:4, col) = qp_C(1:4, 4)*BSCEval(k)
                col = col+1
            end do
            if (hhoCell%l_axis) then
! ---------- diagonal terms
                do k = 1, gbs_cmp
                    Cgphi(1:4, col) = qp_C(1:4, 3)*BSCEval(k)
                    col = col+1
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
    subroutine tranfoMatToSym(ndim, mat, mat_sym)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: mat(*)
        real(kind=8), intent(out) :: mat_sym(6)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   tranform a matrix to matrix symmetrix form
!   In ndim         : dimension of the problem
!   In matrix       : symmetrix matrix to transform (XX, YY, ZZ, XY, XZ, YZ)
!   Out mat_sym     : matrix in form (XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ)
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
! --------------------------------------------------------------------------------------------------
!
        mat_sym = 0.d0
!
        select case (ndim)
        case (3)
            mat_sym(1:3) = mat(1:3)
            mat_sym(4:6) = mat(4:6)*rac2
        case (2)
            mat_sym(1:3) = mat(1:3)
            mat_sym(4) = mat(4)*rac2
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
    subroutine tranfoSymToMat(ndim, mat_sym, mat)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(out) :: mat(*)
        real(kind=8), intent(in) :: mat_sym(6)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   tranform a matrix to matrix symmetrix form
!   In ndim        : dimension of the problem
!   In mat_sym     : matrix in form (XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ)
!   Out matrix     : symmetrix matrix to transform (XX, YY, ZZ, XY, XZ, YZ)
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), parameter :: un_rac2 = 1.d0/sqrt(2.d0)
! --------------------------------------------------------------------------------------------------
!
        select case (ndim)
        case (3)
            mat(1:3) = mat_sym(1:3)
            mat(4:6) = mat_sym(4:6)*un_rac2
        case (2)
            mat(1:3) = mat_sym(1:3)
            mat(4) = mat_sym(4)*un_rac2
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
    subroutine tranfoTensToSym(nb_sig, dsidep, dsidep3D)
!
        implicit none
!
        integer(kind=8), intent(in) :: nb_sig
        real(kind=8), intent(in) :: dsidep(nb_sig, nb_sig)
        real(kind=8), intent(out) :: dsidep3D(6, 6)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   tranform a tensor from nb_sig to 6 and add symetric notation
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, j
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
!
        select case (nb_sig)
        case (6)
            dsidep3D(1:3, 1:3) = dsidep(1:3, 1:3)
            do i = 4, 6
                do j = 1, 3
                    dsidep3D(i, j) = dsidep(i, j)*rac2
                    dsidep3D(j, i) = dsidep(j, i)*rac2
                end do
                do j = 4, 6
                    dsidep3D(i, j) = dsidep(i, j)*2.d0
                    dsidep3D(j, i) = dsidep(j, i)*2.d0
                end do
            end do
        case (4)
            dsidep3D(1:3, 1:3) = dsidep(1:3, 1:3)
            do j = 1, 3
                dsidep3D(4, j) = dsidep(4, j)*rac2
                dsidep3D(j, 4) = dsidep(j, 4)*rac2
            end do
            dsidep3D(4, 4) = dsidep(4, 4)*2.d0
            dsidep3D(5:6, 1:6) = 0.d0
            dsidep3D(1:6, 5:6) = 0.d0
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
end module
