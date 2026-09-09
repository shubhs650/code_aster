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

subroutine te0437(nomopt, nomte)
!
    use HHO_type
    use HHO_basis_module
    use HHO_eval_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_ther_module
    use HHO_utils_module
    use HHO_matrix_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/readVector.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - MTAN_THER_SOURNL
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Quadrature) :: hhoQuad
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_basis_cell) :: hhoBasisCell
    integer(kind=8) :: cbs, fbs, total_dofs, npg, isour, ipg, itime, faces_dofs
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8) :: VoluValuesQP(MSIZE_QP_CELL)
    real(kind=8) :: theta, sour, dsdt, temp_eval
    real(kind=8), dimension(MSIZE_CELL_SCAL) :: temp_T
    real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
    type(HHO_matrix) :: lhs, lhs_cell
!
! --- Get element parameters
!
    call elrefe_info(fami=fami, npg=npg)
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuad)
!
! --- Number of dofs
    call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    ASSERT(total_dofs <= MSIZE_TDOFS_SCAL)
    faces_dofs = total_dofs-cbs
!
    call hhoBasisCell%initialize(hhoCell)
!
    call jevech('PINSTR', 'L', itime)
    theta = zr(itime+2)
    ASSERT(theta < -0.5)
!
! ---- Which option ?
!
    if (nomopt .eq. 'MTAN_THER_SOURNL') then
        call jevech('PSOURNL', 'L', isour)
        if (zk8(isour) (1:7) .eq. '&FOZERO') goto 999
!
! ----- Get real value
!
        temp_T = 0.d0
        call readVector('PTEMPEI', cbs, temp_T, faces_dofs)
!
        do ipg = 1, hhoQuad%nbQuadPoints
            temp_eval = hhoEvalScalCell2(hhoBasisCell, hhoData%cell_degree(), &
                                         hhoQuad%points(1:3, ipg), temp_T)

            call foderi(zk8(isour), temp_eval, sour, dsdt)

            VoluValuesQP(ipg) = -dsdt
        end do
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
! ---- Compute mass matrix
!
    call lhs%initialize(total_dofs, total_dofs, 0.d0)
    call lhs_cell%initialize(cbs, cbs, 0.d0)
!
! ----- Loop on quadrature point
    do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval basis function at the quadrature point
        call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, &
                                 hhoData%cell_degree(), BSCEval)
! --------  Eval massMat
        call hhoComputeLhsMassTher(VoluValuesQP(ipg), hhoQuad%weights(ipg), &
                                   BSCEval, cbs, lhs_cell)
    end do
!
! ----- Copy the lower part
!
    call lhs_cell%copySymU()
    call lhs%copy(lhs_cell, faces_dofs, faces_dofs)
    call lhs_cell%free()
!
! ---- save result
!
    call lhs%write('PMATTTR', ASTER_TRUE)
    call lhs%free()
!
999 continue
!
end subroutine
