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
subroutine te0445(nomopt, nomte)
!
    use HHO_type
    use HHO_utils_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_ther_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_matrix_module
    use FE_algebra_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
#include "blas/daxpy.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - CHAR_THER_EVOL
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Quadrature) :: hhoQuadCellRigi, hhoQuadCellMass
    integer(kind=8) :: cbs, fbs, total_dofs, npg_rigi, npg_mass, itemps
    character(len=8), parameter :: fami_rigi = 'RIGI', fami_mass = 'MASS'
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_matrix) :: gradfull, stab
    real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: rhs_rigi, rhs_mass, rhs
    real(kind=8) :: theta, dtime
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
! --- Get element parameters
!
    call elrefe_info(fami=fami_rigi, npg=npg_rigi)
    call elrefe_info(fami=fami_mass, npg=npg_mass)
!
! --- Number of dofs
    call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    ASSERT(total_dofs <= MSIZE_TDOFS_SCAL)
!
    if (nomopt /= "CHAR_THER_EVOL") then
        ASSERT(ASTER_FALSE)
    end if
!
! --- Initialize quadrature for the rigidity
!
    call hhoQuadCellRigi%initCell(hhoCell, npg_rigi)
    call hhoQuadCellMass%initCell(hhoCell, npg_mass)
!
! --- Compute Operators
!
    if (hhoData%precompute()) then
        call hhoReloadPreCalcTher(hhoCell, hhoData, gradfull, stab)
    else
        call hhoCalcOpTher(hhoCell, hhoData, gradfull, stab)
    end if
!
! --- Compute local rigidity contribution
!
    call hhoLocalRigiTher(hhoCell, hhoData, hhoQuadCellRigi, nomopt, gradfull, &
                          stab, fami_rigi, rhs=rhs_rigi)
!
! --- Compute local mass contribution
!
    call hhoLocalMassTher(hhoCell, hhoData, hhoQuadCellMass, fami_mass, rhs=rhs_mass)
!
! --- Compute rhs = 1/dt * rhs_mass - (1-theta) * rhs_rigi
!
    call jevech('PINSTR', 'L', itemps)
    dtime = zr(itemps+1)
    theta = zr(itemps+2)
!
    rhs = 0.d0
    call daxpy_1(cbs, 1.d0/dtime, rhs_mass(total_dofs-cbs+1:total_dofs), &
                 rhs(total_dofs-cbs+1:total_dofs))
    call daxpy_1(total_dofs, -(1.d0-theta), rhs_rigi, rhs)
!
! --- Save rhs
!
    call writeVector('PVECTTR', total_dofs, rhs)
!
    call gradfull%free()
    call stab%free()
!
end subroutine
