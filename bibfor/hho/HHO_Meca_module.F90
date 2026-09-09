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
module HHO_Meca_module
!
    use Behaviour_module
    use Behaviour_type
    use FE_algebra_module
    use HHO_algebra_module
    use HHO_basis_module
    use HHO_compor_module
    use HHO_Dirichlet_module
    use HHO_eval_module
    use HHO_gradrec_module, only: hhoGradRecSymMat, hhoGradRecFullMatFromVec, &
                                  hhoGradRecSymFullMatFromVec, hhoGradRecVec, &
                                  hhoGradRecFullMat, hhoGradRecSymFullMat
    use HHO_init_module
    use HHO_LargeStrainMeca_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_SmallStrainMeca_module
    use HHO_stabilization_module, only: hhoStabVec, hdgStabVec, hhoStabSymVec
    use HHO_type
    use HHO_utils_module
    use MaterialPara_module
    use MaterialPara_type
    use NonLin_Datastructure_type
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/calcul.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/infniv.h"
#include "asterfort/isfonc.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/readVector.h"
#include "asterfort/sigtopk1.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/dsyr.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - mechanics
!
! Specific routines for mechanics
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: HHO_Meca_State
    public :: hhoMecaInit, hhoLocalContribMeca, hhoCalcStabCoeffMeca
    public :: hhoCalcOpMeca, hhoReloadPreCalcMeca, hhoLocalForcNoda
    public :: YoungModulus, hhoLocalMassMeca
!
! --------------------------------------------------------------------------------------------------
!
    type HHO_Meca_State
!
        aster_logical :: l_debug = ASTER_FALSE
! ----- Time : prev, curr, incr
        real(kind=8) :: time_prev = 0.d0
        real(kind=8) :: time_curr = 0.d0
        real(kind=8) :: time_incr = 0.d0
! ----- Displacement
        real(kind=8), dimension(MSIZE_TDOFS_VEC) :: depl_prev = 0.d0
        real(kind=8), dimension(MSIZE_TDOFS_VEC) :: depl_curr = 0.d0
        real(kind=8), dimension(MSIZE_TDOFS_VEC) :: depl_incr = 0.d0
! ----- Gradient reconstruction and stabilisation
        type(HHO_matrix) :: grad
        type(HHO_matrix) :: stab
! ----- member function
    contains
        procedure, pass :: initialize => initialize_meca
        procedure, pass :: free => free_meca
!
    end type HHO_Meca_State
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMecaInit(model, list_load, list_func_acti, hhoField)
!
        implicit none
!
        character(len=24), intent(in) :: model
        character(len=19), intent(in) :: list_load
        integer(kind=8), intent(in) :: list_func_acti(*)
        type(HHO_Field), intent(out) :: hhoField
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Non-linear mechanics
!
! Initializations
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  list_load        : name of datastructure for list of loads
! In  list_func_acti   : list of active functionnalities
! Out hhoField         : fields for HHO
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
        call infniv(ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'HHO2_8')
        end if
!
! --- Prepare fields for Dirichlet loads
!
        hhoField%fieldCineFunc = '&&HHOMEC.CINEFUNC'
        hhoField%fieldCineVale = '&&HHOMEC.CINEVALE'
!
        if (isfonc(list_func_acti, 'DIRI_CINE')) then
            call hhoDiriFuncPrepare(model(1:8), list_load, hhoField)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoReloadPreCalcMeca(hhoCell, hhoData, l_largestrains, gradfull, stab)
!
        implicit none
!
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_Cell), intent(in) :: hhoCell
        aster_logical, intent(in) :: l_largestrains
        type(HHO_matrix), intent(out) :: gradfull
        type(HHO_matrix), intent(out), optional :: stab
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - Reload Precomputation of operators
!
! In  hhoCell         : hho Cell
! In hhoData          : information about the HHO formulation
! In l_largestrains   : large strains ?
! Out gradfull        : full gradient for mechanics
! Out stab            : stabilization for mechanics
! --------------------------------------------------------------------------------------------------
!
! --- Local variables
!
        integer(kind=8) :: cbs, fbs, total_dofs, gbs
        type(HHO_matrix) :: gradfullvec, stabvec
!
        call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)
!
! -------- Reload gradient
!
        call gradfullvec%initialize(gbs, total_dofs)
        call gradfullvec%read('PCHHOGT', ASTER_FALSE)
!
        if (l_largestrains) then
            call hhoGradRecFullMatFromVec(hhoCell, hhoData, gradfullvec, gradfull)
        else
            call hhoGradRecSymFullMatFromVec(hhoCell, hhoData, gradfullvec, gradfull)
        end if
!
        call gradfullvec%free()
!
! -------- Reload stabilization
        if (present(stab)) then
            call stabvec%initialize(total_dofs, total_dofs)
            call stabvec%read('PCHHOST', ASTER_TRUE)
            call MatScal2Vec(hhoCell, hhoData, stabvec, stab)
            call stabvec%free()
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLocalContribMeca(hhoCell, hhoData, hhoQuadCellRigi, hhoMecaState, hhoCS, &
                                   lhs, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellRigi
        type(HHO_Meca_State), intent(in) :: hhoMecaState
        type(HHO_Compor_State), intent(inout) :: hhoCS
        type(HHO_matrix), intent(out) :: lhs
        real(kind=8), intent(out) :: rhs(MSIZE_TDOFS_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the local contribution for mechanics
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   In hhoMecaState : mechanical state
!   IO hhoCS        : behaviour
!   Out lhs         : local contribution (lhs)
!   Out rhs         : local contribution (rhs)
! --------------------------------------------------------------------------------------------------
!
        aster_logical :: l_rigi_meca, l_vari
        integer(kind=8) :: cbs, fbs, total_dofs

! ----- Verif compor
        ASSERT(.not. hhoCS%c_plan)
        l_rigi_meca = (hhoCS%option == "RIGI_MECA")
        l_vari = L_VARI(hhoCS%option)

! ----- Number of dofs
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)

! ----- Initialization
        call lhs%initialize(total_dofs, total_dofs, 0.0)
        rhs = 0.d0
        if (.not. l_vari) then
!           not accessed but expected to be (lgpg, *)
            AS_ALLOCATE(vr=hhoCS%vari_curr, size=max(1, hhoCS%lgpg))
        end if

!
        if (hhoCS%l_largestrain) then
!
! --- large strains and use gradient
!
            call hhoLargeStrainLCMeca(hhoCell, hhoData, hhoQuadCellRigi, hhoCS, hhoMecaState%grad, &
                                      hhoMecaState%time_prev, hhoMecaState%time_curr, &
                                      hhoMecaState%depl_prev, hhoMecaState%depl_curr, &
                                      lhs, rhs)
        else
!
! --- small strains and use symmetric gradient
!
            if (l_rigi_meca) then
                call hhoMatrElasMeca(hhoCell, hhoData, hhoQuadCellRigi, hhoCS, hhoMecaState%grad, &
                                     hhoMecaState%time_curr, lhs)
                hhoCS%codret = 0
            else
                call hhoSmallStrainLCMeca(hhoCell, hhoData, hhoQuadCellRigi, hhoCS, &
                                          hhoMecaState%grad, &
                                          hhoMecaState%time_prev, hhoMecaState%time_curr, &
                                          hhoMecaState%depl_prev, hhoMecaState%depl_incr, &
                                          lhs, rhs)
            end if
        end if
!
        if (.not. l_vari) then
            AS_DEALLOCATE(vr=hhoCS%vari_curr)
        end if
!
! --- test integration of the behavior
!
        if (hhoCS%codret .ne. 0) goto 999
!
! --- add stabilization
!
        call hhoCalcStabCoeffMeca(hhoData, hhoCS%BEHInteg%materPara%schemePara%fami, &
                                  hhoMecaState%time_curr, hhoQuadCellRigi)
!
        if (L_VECT(hhoCS%option)) then
            call hho_dsymv_U(hhoData%coeff_stab(), hhoMecaState%stab, hhoMecaState%depl_curr, &
                             1.d0, rhs)
        end if
!
        if (L_MATR(hhoCS%option)) then
            call lhs%add(hhoMecaState%stab, hhoData%coeff_stab())
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
    subroutine hhoCalcOpMeca(hhoCell, hhoData, l_largestrains, gradfull, stab)
!
        implicit none
!
        type(HHO_Cell), intent(inout) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        aster_logical, intent(in) :: l_largestrains
        type(HHO_matrix), intent(out) :: gradfull
        type(HHO_matrix), intent(out), optional :: stab
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
        type(HHO_matrix) :: gradrec_scal
!        real(kind=8) :: gradrec_sym(MSIZE_CELL_VEC, MSIZE_TDOFS_VEC)
!
! --------------------------------------------------------------------------------------------------
!
        call hhoInitFacesOfCell(hhoCell)
!
        if (l_largestrains) then
!
! ----- Compute Gradient reconstruction
            call hhoGradRecFullMat(hhoCell, hhoData, gradfull)
        else
!
! ----- Compute Symmetric Gradient reconstruction
            call hhoGradRecSymFullMat(hhoCell, hhoData, gradfull)
        end if
!
! ----- Compute Stabilizatiion
        if (present(stab)) then
            if (hhoData%cell_degree() <= hhoData%face_degree()) then
                call hhoGradRecVec(hhoCell, hhoData, gradrec_scal)
                call hhoStabVec(hhoCell, hhoData, gradrec_scal, stab)
                call gradrec_scal%free()
!               call hhoGradRecSymMat(hhoCell, hhoData, gradrec_sym)
!               call hhoStabSymVec(hhoCell, hhoData, gradrec_sym, stab)
!               call gradrec_sym%free()
            else if (hhoData%cell_degree() == (hhoData%face_degree()+1)) then
                call hdgStabVec(hhoCell, hhoData, stab)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function YoungModulus(fami, imate, time, hhoQuad) result(coeff)
!
        implicit none
!
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: imate
        real(kind=8), intent(in) :: time
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8) :: coeff
!
! --------------------------------------------------------------------------------------------------
!
!   Compute the average Young modulus
!   In fami         : familly of quadrature points (of hhoQuadCellRigi)
!   In npg          : number of quadrature points
!   In imate        : materiau code
! --------------------------------------------------------------------------------------------------
!
        type(Behaviour_Integ) :: BEHInteg
        character(len=16) :: elas_keyword
        integer(kind=8) :: elas_id, ipg
        real(kind=8) :: e
!
        coeff = 0.d0
! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
!
        call get_elas_id(imate, elas_id, elas_keyword)
!
        do ipg = 1, hhoQuad%nbQuadPoints
            BEHInteg%behavESVA%behavESVAGeom%coorElga(ipg, 1:3) = hhoQuad%points(1:3, ipg)
            call get_elas_para(fami, imate, '+', ipg, 1, elas_id, elas_keyword, &
                               e_=e, time=time, BEHInteg=BEHInteg)
            coeff = coeff+e
        end do
!
        coeff = coeff/real(hhoQuad%nbQuadPoints, kind=8)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalcStabCoeffMeca(hhoData, fami, time, hhoQuad)
!
        implicit none
!
        type(HHO_Data), intent(inout) :: hhoData
        character(len=8) :: fami
        real(kind=8), intent(in) :: time
        type(HHO_Quadrature), intent(in) :: hhoQuad
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - Evaluate stabilzation coefficient
!
! In hhoData          : information about the HHO formulation
! --------------------------------------------------------------------------------------------------
!
! --- Local variables
!
        integer(kind=8) :: jmate, imate
!
        if (hhoData%adapt()) then
            call jevech('PMATERC', 'L', jmate)
            imate = zi(jmate-1+1)
            call hhoData%setCoeffStab(10.d0*YoungModulus(fami, imate, time, hhoQuad))
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_meca(this, hhoCell, hhoDataMk, hhoComporState, hhoDataGv)
!
        implicit none
!
        class(HHO_Meca_State), intent(inout) :: this
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoDataMk
        type(HHO_Data), intent(in), optional :: hhoDataGv
        type(HHO_Compor_State), intent(in) :: hhoComporState
!
! --------------------------------------------------------------------------------------------------
!
!  initialize HHO_MECA_STATE
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iinstm, iinstp, iret, num_tot, num_mk
        integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs, iFace, iDof
        integer(kind=8) :: gv_cbs, gv_fbs, gv_total_dofs, total_dofs
        real(kind=8) :: tmp_prev(MSIZE_TDOFS_MIX), tmp_incr(MSIZE_TDOFS_MIX)
        aster_logical :: forc_noda, pilo
!
        forc_noda = hhoComporState%option == "FORC_NODA"
        pilo = hhoComporState%option(1:4) == "PILO"
        if (hhoComporState%option .ne. "RIGI_MECA" &
            .and. hhoComporState%option .ne. "REFE_FORC_NODA") then
            if (.not. forc_noda .and. .not. pilo) then
                call jevech('PINSTMR', 'L', iinstm)
                call jevech('PINSTPR', 'L', iinstp)
                this%time_curr = zr(iinstp)
                this%time_prev = zr(iinstm)
                this%time_incr = this%time_curr-this%time_prev
            end if
!
            call hhoMecaDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs)
!
            if (hhoComporState%typmod(2) == "HHO") then
                if (forc_noda) then
                    call readVector('PDEPLAR', mk_total_dofs, this%depl_prev)
                else
! ----------------- get displacement in T-
                    call readVector('PDEPLMR', mk_total_dofs, this%depl_prev)

! ----------------- get increment displacement beetween T- and T+
                    if (pilo) then
                        call readVector('PDDEPLR', mk_total_dofs, this%depl_incr)
                    else
                        call readVector('PDEPLPR', mk_total_dofs, this%depl_incr)
                    end if
                end if
            else
                ! GRAD_VARI
                call hhoTherDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs)
                total_dofs = mk_total_dofs+gv_total_dofs+gv_cbs
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
                num_tot = 0
                num_mk = 0
                do iFace = 1, hhoCell%nbfaces
                    do iDof = 1, mk_fbs
                        num_tot = num_tot+1
                        num_mk = num_mk+1
                        this%depl_prev(num_mk) = tmp_prev(num_tot)
                        if (.not. forc_noda) then
                            this%depl_incr(num_mk) = tmp_incr(num_tot)
                        end if
                    end do
                    num_tot = num_tot+gv_fbs
                end do
                do iDof = 1, mk_cbs
                    num_tot = num_tot+1
                    num_mk = num_mk+1
                    this%depl_prev(num_mk) = tmp_prev(num_tot)
                    if (.not. forc_noda) then
                        this%depl_incr(num_mk) = tmp_incr(num_tot)
                    end if
                end do
            end if

! --------- compute displacement in T+
            call dcopy_1(mk_total_dofs, this%depl_prev, this%depl_curr)
            if (.not. forc_noda) then
                call daxpy_1(mk_total_dofs, 1.d0, this%depl_incr, this%depl_curr)
            end if

        else if (hhoComporState%option == "RIGI_MECA") then
            call tecach('ONO', 'PINSTR', 'L', iret, iad=iinstp)
            if (iinstp .ne. 0) then
                this%time_curr = zr(iinstp)
            end if

        else if (hhoComporState%option == "REFE_FORC_NODA") then
            !! Nothing to load
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
    subroutine free_meca(this)
!
        implicit none
!
        class(HHO_Meca_State), intent(inout) :: this
!
! --------------------------------------------------------------------------------------------------
!
!  free HHO_MECA_STATE
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
    subroutine hhoLocalMassMeca(hhoCell, hhoData, hhoQuadCellMass, fami, mass)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellMass
        character(len=8), intent(in) :: fami
        type(HHO_matrix), intent(out) :: mass
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the local mass contribution for mechanics
!   RHS = (rho * uT, vT)_T
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   In fami         : familly of quadrature points (of hhoQuadCellRigi)
!   Out mass         : local contribution (lhs)
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_matrix) :: mass_vec
        character(len=32) :: phenom
        integer(kind=8) :: cbs, fbs, total_dofs, faces_dofs
        integer(kind=8) :: jmate, ipg, icodre(3), dimMatScal
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8) :: coorpg(3), weight, rho, rho_(1), coeff
        real(kind=8) :: mass_scal(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
        blas_int :: b_incx, b_lda, b_n
!
! --- Get input fields
!
        call jevech('PMATERC', 'L', jmate)
        call rccoma(zi(jmate), 'ELAS', 1, phenom, icodre(1))
!
! --- number of dofs
!
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        faces_dofs = total_dofs-cbs
        dimMatScal = cbs/hhoCell%ndim
!
        b_n = to_blas_int(dimMatScal)
        b_incx = to_blas_int(1)
        b_lda = to_blas_int(MSIZE_CELL_SCAL)
!
! -- initialization
!
        call mass%initialize(total_dofs, total_dofs, 0.0)
        mass_scal = 0.d0
!
        call hhoBasisCell%initialize(hhoCell)
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCellMass%nbQuadPoints
            coorpg(1:3) = hhoQuadCellMass%points(1:3, ipg)
            weight = hhoQuadCellMass%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
            call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%cell_degree(), BSCEval)
!
! -------- Compute behavior
!
            call rcvalb(fami, ipg, 1, '+', zi(jmate), &
                        ' ', phenom, 0, ' ', [0.d0], &
                        1, 'RHO', rho_, icodre(1), 1)
            rho = rho_(1)
!
! -------- Compute mass
!
            coeff = rho*weight
            call dsyr('U', b_n, coeff, BSCEval, b_incx, mass_scal, b_lda)
!
        end do
!
! ----- Copy the lower part
!
        call hhoCopySymPartMat('U', mass_scal(1:dimMatScal, 1:dimMatScal))
        call MatCellScal2Vec(hhoCell, hhoData, mass_scal, mass_vec)
        call mass%copy(mass_vec, faces_dofs, faces_dofs)
        call mass_vec%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLocalForcNoda(hhoCell, hhoData, hhoQuadCellRigi, hhoMecaState, hhoCS, &
                                stress, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellRigi
        type(HHO_Meca_State), intent(in) :: hhoMecaState
        type(HHO_Compor_State), intent(in) :: hhoCS
        real(kind=8), intent(in) :: stress(*)
        real(kind=8), intent(out) :: rhs(MSIZE_TDOFS_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the local residual for a given state
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   In fami         : familly of quadrature points (of hhoQuadCellRigi)
!   Out rhs         : local contribution (rhs)
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        integer(kind=8) :: cbs, fbs, total_dofs, gbs, gbs_sym, gbs_axis
        integer(kind=8) :: ipg, ncomp, gbs_curr, gbs_cmp, cbs_cmp, faces_dofs
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL), rhs_axis(MSIZE_CELL_SCAL)
        real(kind=8) :: coorpg(3), weight
        real(kind=8) :: Cauchy_curr(6), PK1_curr(3, 3), G_curr(3, 3), F_curr(3, 3)
        real(kind=8), dimension(MSIZE_CELL_MAT) :: bT, G_curr_coeff
!
        rhs = 0.d0
        bT = 0.d0
        rhs_axis = 0.d0
        ncomp = hhoCS%nbsigm
!
! ----- init basis
!
        call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                           gbs, gbs_sym, gbs_axis)
        faces_dofs = total_dofs-cbs
        gbs_cmp = (gbs-gbs_axis)/(hhoCell%ndim*hhoCell%ndim)
        cbs_cmp = cbs/hhoCell%ndim
        call hhoBasisCell%initialize(hhoCell)
!
! --- Compute local contribution
!
        if (hhoCS%l_largestrain) then
            call hhoMecaState%grad%dot(hhoMecaState%depl_curr, G_curr_coeff)
            gbs_curr = gbs
        else
            gbs_curr = gbs_sym
        end if
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCellRigi%nbQuadPoints
            coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
            weight = hhoQuadCellRigi%weights(ipg)
!
! -------- tranform sigm in symmetric form
!
            call tranfoMatToSym(hhoCell%ndim, stress((ipg-1)*ncomp+1:ipg*ncomp), Cauchy_curr)
!
! --------- Eval basis function at the quadrature point
!
            call hhoBasisCell%BSEval(coorpg(1:3), 0, &
                                     max(hhoData%grad_degree(), hhoData%cell_degree()), &
                                     BSCEval)
!
! --------- Eval gradient at T- and T+
!
            if (hhoCS%l_largestrain) then
                G_curr = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, G_curr_coeff)
!
! --------- Eval gradient of the deformation at T- and T+
!
                call hhoCalculF(G_curr, F_curr)
!
                call sigtopk1(hhoCell%ndim, Cauchy_curr, F_curr, PK1_curr)
!
                call hhoComputeRhsLarge(hhoCell, PK1_curr, weight, BSCEval, gbs, bT)
            else
!
                call hhoComputeRhsSmall(hhoCell, Cauchy_curr, weight, BSCEval, gbs_cmp, bT)
            end if
        end do
!
        call hho_dgemv_T(1.d0, hhoMecaState%grad, bT, 0.d0, rhs)
!
! --- add stabilization
!
        call hhoCalcStabCoeffMeca(hhoData, &
                                  hhoCS%BEHInteg%materPara%schemePara%fami, 0.d0, hhoQuadCellRigi)
!
        call hho_dsymv_U(hhoData%coeff_stab(), hhoMecaState%stab, hhoMecaState%depl_curr, &
                         1.d0, rhs)
!
    end subroutine
!
end module
