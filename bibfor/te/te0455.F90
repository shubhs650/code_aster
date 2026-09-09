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
subroutine te0455(option, nomte)
!
    use Behaviour_module, only: behaviourOption
    use HHO_compor_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_matrix_module
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use HHO_utils_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/nmtstm.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: nomte, option
!
! --------------------------------------------------------------------------------------------------
!
! HHO
! Mechanics - STAT_NON_LINE
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI", typmod2 = "HHO"
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
    type(HHO_Quadrature) :: hhoQuadCellRigi
    type(HHO_Compor_State) :: hhoCS
    integer(kind=8) :: cbs, fbs, total_dofs
    integer(kind=8) :: jmatt, npg
    integer(kind=8) :: jcret
    aster_logical :: l_largestrains, lMatr, lVect, lSigm, lVari, matsym
    real(kind=8) :: rhs(MSIZE_TDOFS_VEC)
    type(HHO_matrix) :: lhs
!
! --------------------------------------------------------------------------------------------------
!
    if (option /= "RIGI_MECA_TANG" .and. option /= "RIGI_MECA_ELAS" .and. &
        option /= "RIGI_MECA" .and. option /= "FULL_MECA" .and. option /= "RAPH_MECA") then
        ASSERT(ASTER_FALSE)
    end if

! - Get element parameters
    call elrefe_info(fami=fami, npg=npg)

! - Get HHO data on the modelisation
    call hhoInfoInitCell(hhoCell, hhoData)

! - Number of dofs
    call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    ASSERT(total_dofs <= MSIZE_TDOFS_VEC)

! - Initialize quadrature for the rigidity
    call hhoQuadCellRigi%initCell(hhoCell, npg)

! - Type of finite element
    call hhoCS%initialize(fami, option, hhoCell%ndim, hhoCell%barycenter, typmod2)
    call hhoMecaState%initialize(hhoCell, hhoData, hhoCS)
!
    if (option .ne. "RIGI_MECA") then
        l_largestrains = hhoCS%l_largestrain
        call behaviourOption(option, hhoCS%compor, lMatr, lVect, lVari, &
                             lSigm, hhoCS%codret)

    else
        l_largestrains = ASTER_FALSE
        lMatr = ASTER_TRUE
        lSigm = ASTER_FALSE
        lVect = ASTER_FALSE
    end if

! - Compute Operators
    if (hhoData%precompute()) then
        call hhoReloadPreCalcMeca(hhoCell, hhoData, l_largestrains, hhoMecaState%grad, &
                                  hhoMecaState%stab)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, l_largestrains, hhoMecaState%grad, &
                           hhoMecaState%stab)
    end if

! - Compute local contribution
    call hhoLocalContribMeca(hhoCell, hhoData, hhoQuadCellRigi, hhoMecaState, hhoCS, &
                             lhs, rhs)

! - Save return code
    if (lSigm) then
        call jevech('PCODRET', 'E', jcret)
        zi(jcret) = hhoCS%codret
    end if

! - Save rhs
    if (lVect) then
        call writeVector('PVECTUR', total_dofs, rhs)
    end if

! - Save of lhs
    if (lMatr) then
        if (option .ne. "RIGI_MECA") then
            call nmtstm(hhoCS%carcri, jmatt, matsym)
        else
            matsym = ASTER_TRUE
        end if
!
        if (matsym) then
            call lhs%write('PMATUUR', ASTER_TRUE)
        else
            call lhs%write('PMATUNS', ASTER_FALSE)
        end if
    end if
!
    call lhs%free()
    call hhoMecaState%free()
!
end subroutine
