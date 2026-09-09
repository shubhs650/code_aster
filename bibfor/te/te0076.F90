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
subroutine te0076(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_eval_module
    use coorSyst_module, only: hasOrieField
    implicit none
!
#include "asterfort/nlcomp.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/utmess.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: THER_*
!
! Options: RIGI_THER
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    character(len=8), parameter :: famiR = "RIGI"
    integer(kind=8) :: kpg, jvMaterc, jvInstr
    integer(kind=8) :: propCode(1)
    character(len=16) :: therKeyword
    real(kind=8) :: time
    real(kind=8) :: rigi(MAX_BS_CG, MAX_BS_CG)
    real(kind=8) :: valQPK(3, 3, MAX_QP)
!
! --------------------------------------------------------------------------------------------------
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, famiR)
    call FEBasis%initCell(FECell)
!
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PINSTR', 'L', jvInstr)
    time = zr(jvInstr)
!
    call rccoma(zi(jvMaterc), 'THER', 1, therKeyword, propCode(1))
    if (therKeyword == "THER_ORTH") then
        if (.not. hasOrieField()) then
            call utmess('F', 'THERMIQUE1_3')
        end if
    end if
!
    valQPK = 0.d0
    do kpg = 1, FEQuadCell%nbQuadPoints
        call nlcomp(therKeyword, famiR, kpg, jvMaterc, FECell%ndim, FEQuadCell%points(1:3, kpg), &
                    time, 0.d0, valQPK(1:3, 1:3, kpg))
    end do
!
    call FEStiffJacoScal(FEQuadCell, FEBasis, valQPK, rigi)
    call writeMatrix("PMATTTR", FEBasis%size, FEBasis%size, ASTER_TRUE, rigi)
!
end subroutine
