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
subroutine te0078(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_rhs_module
    use FE_eval_module
    use coorSyst_module, only: hasOrieField
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/nlcomp.h"
#include "asterfort/utmess.h"
#include "asterfort/writeVector.h"
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
! Options: CHAR_THER_EVOL
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadRigi, FEQuadMass
    type(FE_basis) :: FEBasis
!
    integer(kind=8), parameter :: nbProp = 1
    integer(kind=8) :: propCode(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'RHO_CP'/)
    character(len=16) :: therKeyword
    real(kind=8) :: valQPM(MAX_QP), tpg, dtpg(3), flux(3), BGSEval(3, MAX_BS_CG)
    real(kind=8) :: resi_f(MAX_BS_CG), resi_m(MAX_BS_CG), resi(MAX_BS_CG)
    real(kind=8) :: cp, propVale(1), Kglo(3, 3), time, deltat, theta
    integer(kind=8) :: kpg, jvMaterc, jvInstr
    real(kind=8), pointer :: temp(:) => null()
    character(len=8), parameter :: famiR = "RIGI"
    character(len=8), parameter :: famiM = "MASS"
!
! --------------------------------------------------------------------------------------------------
!
    call FECell%init()
    call FEBasis%initCell(FECell)
    call FEQuadMass%initCell(FECell, famiM)
    call FEQuadRigi%initCell(FECell, famiR)
!
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PINSTR', 'L', jvInstr)
    call jevech('PTEMPER', 'L', vr=temp)
!
    time = zr(jvInstr)
    deltat = zr(jvInstr+1)
    theta = zr(jvInstr+2)
!
    call rccoma(zi(jvMaterc), 'THER', 1, therKeyword, propCode(1))
    if (therKeyword == "THER_ORTH") then
        if (.not. hasOrieField()) then
            call utmess('F', 'THERMIQUE1_3')
        end if
    end if
!
    resi_f = 0.d0
    do kpg = 1, FEQuadRigi%nbQuadPoints
        BGSEval = FEBasis%grad(FEQuadRigi%points_param(1:3, kpg), FEQuadRigi%jacob(1:3, 1:3, kpg))
!
        dtpg = FEEvalGradVec(FEBasis, temp, FEQuadRigi%points_param(1:3, kpg), BGSEval)
!
        call nlcomp(therKeyword, famiR, kpg, jvMaterc, FECell%ndim, FEQuadRigi%points(1:3, kpg), &
                    time, 0.d0, Kglo, dtp_=dtpg, fluglo_=flux)
!
        call FEStiffResiScalAdd(FEBasis, BGSEval, FEQuadRigi%weights(kpg), flux, resi_f)
    end do
!
    do kpg = 1, FEQuadMass%nbQuadPoints
        call rcvalb(famiM, kpg, 1, '+', &
                    zi(jvMaterc), ' ', therKeyword, &
                    1, 'INST', [time], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        cp = propVale(1)
!
        tpg = FEEvalFuncRScal(FEBasis, temp, FEQuadMass%points_param(1:3, kpg))
        ValQPM(kpg) = cp*tpg
    end do
!
    call FeMakeRhsScal(FEQuadMass, FEBasis, ValQPM, resi_m)
!
    resi = (theta-1.0d0)*resi_f+resi_m/deltat
!
    call writeVector('PVECTTR', FEBasis%size, resi)
!
end subroutine
