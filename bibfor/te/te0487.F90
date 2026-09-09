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
subroutine te0487(nomopt, nomte)
!
    use HHO_type
    use HHO_basis_module
    use HHO_utils_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_eval_module
    use HHO_ther_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_matrix_module
    use HHO_algebra_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/readVector.h"
#include "asterfort/writeVector.h"
#include "asterfort/rccoma.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - FLUX_ELGA
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_basis_cell) :: hhoBasisCell
    type(HHO_Quadrature) :: hhoQuadCellRigi
    integer(kind=8) :: cbs, fbs, total_dofs, npg, deca, gbs, cell_offset
    integer(kind=8) :: ipg, icodre(3), jtemps, jmate
    character(len=8), parameter :: fami = 'RIGI'
    character(len=32) :: phenom
    type(HHO_matrix) :: gradrec
    real(kind=8), dimension(3*MSIZE_QP_CELL) :: flux
    real(kind=8) :: module_tang(3, 3), G_curr(3), sig_curr(3), dsig_curr(3)
    real(kind=8) :: coorpg(3), weight, time_curr, temp_eval_curr
    real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
    real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: temp_curr
    real(kind=8), dimension(MSIZE_CELL_VEC) :: G_curr_coeff
!
! --- Get element parameters
!
    call elrefe_info(fami=fami, npg=npg)
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCellRigi)
!
! --- Number of dofs
    call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)
    ASSERT(cbs <= MSIZE_CELL_SCAL)
    ASSERT(fbs <= MSIZE_FACE_SCAL)
    ASSERT(total_dofs <= MSIZE_TDOFS_SCAL)
    cell_offset = total_dofs-cbs+1
!
    if (nomopt /= "FLUX_ELGA") then
        ASSERT(ASTER_FALSE)
    end if
!
! --- Compute Operators
!
    if (hhoData%precompute()) then
!
        call hhoReloadPreCalcTher(hhoCell, hhoData, gradrec)
    else
        call hhoCalcOpTher(hhoCell, hhoData, gradrec)
    end if
!
! --- Get input fields
!
    call jevech('PMATERC', 'L', jmate)
!
    call rccoma(zi(jmate), 'THER', 1, phenom, icodre(1))
!
    call jevech('PINSTR', 'L', jtemps)
    time_curr = zr(jtemps)
!
! -- initialization
!
    flux = 0.d0
!
    G_curr_coeff = 0.d0
!
    call hhoBasisCell%initialize(hhoCell)
!
! --- compute temp in T+
!
    temp_curr = 0.d0
    call readVector('PTEMPER', total_dofs, temp_curr)
!
! ----- compute G_curr = gradrec * temp_curr
!
    call gradrec%dot(temp_curr, G_curr_coeff)
!
! ----- Loop on quadrature point
!
    deca = 1
    do ipg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
        weight = hhoQuadCellRigi%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
        call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
! --------- Eval gradient at T+
!
        G_curr = hhoEvalVecCell(hhoCell%ndim, gbs, BSCEval, G_curr_coeff)
!
! --------- Eval temperature at T+
!
        temp_eval_curr = hhoEvalScalCell2( &
                         hhoBasisCell, hhoData%cell_degree(), coorpg(1:3), &
                         temp_curr(cell_offset:))
!
! ------- Compute behavior
!
        call hhoComputeBehaviourTher(phenom, fami, ipg, hhoCell%ndim, time_curr, &
                                     jmate, coorpg, temp_eval_curr, G_curr, sig_curr, &
                                     dsig_curr, module_tang)
!
        flux(deca:deca+hhoCell%ndim) = -sig_curr(1:hhoCell%ndim)
        deca = deca+hhoCell%ndim
    end do
!
! --- Save fluxes
!
    call writeVector('PFLUXPG', hhoCell%ndim*npg, flux)
!
    call gradrec%free()
!
end subroutine
