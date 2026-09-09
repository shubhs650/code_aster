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
subroutine te0485(option, nomte)
!
    use Behaviour_module, only: behaviourOption
    use HHO_compor_module
    use HHO_GV_module
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
! Mechanics - STAT_NON_LINE - GRAD_VARI
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI", typmod2 = "GRADVARI"
    type(HHO_Data) :: hhoDataMk, hhoDataGv
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
    type(HHO_GV_State) :: hhoGVState
    type(HHO_Compor_State) :: hhoCS
    type(HHO_Quadrature) :: hhoQuadCellRigi
    integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs
    integer(kind=8) :: gv_cbs, gv_fbs, gv_total_dofs, total_dofs
    integer(kind=8) :: jmatt, npg, jcret
    aster_logical :: lMatr, lVect, lSigm, lVari, matsym
    real(kind=8) :: rhs(MSIZE_TDOFS_MIX)
    type(HHO_matrix) :: lhs
!
! --------------------------------------------------------------------------------------------------
!
    if (option /= "RIGI_MECA_TANG" .and. &
        option /= "FULL_MECA" .and. &
        option /= "FORC_NODA" .and. &
        option /= "RAPH_MECA") then
        ASSERT(ASTER_FALSE)
    end if

! - Get element parameters
    call elrefe_info(fami=fami, npg=npg)

! - Get HHO data on the modelisation
    call hhoInfoInitCell(hhoCell, hhoDataMk)
    call hhoDataGVInit(hhoDataGv)

! - Number of dofs
    call hhoMecaDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs)
    call hhoTherDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs)
    total_dofs = mk_total_dofs+gv_total_dofs+gv_cbs

! - Initialize quadrature for the rigidity
    call hhoQuadCellRigi%initCell(hhoCell, npg)

! - Properties of behaviour
    call hhoCS%initialize(fami, option, hhoCell%ndim, hhoCell%barycenter, typmod2)

! - Initialize displacement, vari, ...
    call hhoMecaState%initialize(hhoCell, hhoDataMk, hhoCS, hhoDataGv)
    call hhoGVState%initialize(hhoCell, hhoDataMk, hhoDataGv, hhoCS)

! - Compute Operators
    call hhoCalcOpGv(hhoCell, hhoDataMk, hhoDataGv, hhoCS%l_largestrain, &
                     hhoMecaState, hhoGvState)

! - Compute local contribution
    call hhoGradVariLC(hhoCell, hhoDataMk, hhoDataGv, hhoQuadCellRigi, &
                       hhoMecaState, hhoCS, hhoGVState, lhs, rhs)

    call behaviourOption(option, hhoCS%compor, &
                         lMatr, lVect, &
                         lVari, lSigm)

! - Save return code
    if (lSigm) then
        call jevech('PCODRET', 'E', jcret)
        zi(jcret) = hhoCS%codret
    end if

! - Save rhs
    if (lVect .or. option == "FORC_NODA") then
        call writeVector('PVECTUR', total_dofs, rhs)
    end if

! - Save of lhs
    if (lMatr) then
        call nmtstm(hhoCS%carcri, jmatt, matsym)
        if (matsym) then
            call lhs%write('PMATUUR', ASTER_TRUE)
        else
            call lhs%write('PMATUNS', ASTER_FALSE)
        end if
    end if
!
    call lhs%free()
    call hhoMecaState%free()
    call hhoGVState%free()
!
end subroutine
