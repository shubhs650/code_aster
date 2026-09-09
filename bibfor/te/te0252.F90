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
subroutine te0252(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_rhs_module
    use FE_eval_module
    use coorSyst_module, only: hasOrieField
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcfode.h"
#include "asterfort/rcvalb.h"
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
! Options: MASS_THER_RESI
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
    integer(kind=8), parameter :: nbProp = 1
    integer(kind=8) :: propCode(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'CHALHYDR'/)
    character(len=16) :: therKeyword, relaName
    real(kind=8) :: valQP(MAX_QP), tpgi, r8bid
    real(kind=8) :: resi(MAX_BS_CG)
    real(kind=8) :: propVale(1)
    integer(kind=8) :: kpg, jvMaterc
    integer(kind=8) :: ifon(6), nbDof
    aster_logical :: aniso
    character(len=16), pointer :: compor(:) => null()
    real(kind=8), pointer :: hydrgp(:) => null()
    real(kind=8), pointer :: temper(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, "MASS")
    call FEBasis%initCell(FECell)
    nbDof = FEBasis%size
!
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PTEMPEI', 'L', vr=temper)
    call jevech('PMATERC', 'L', jvMaterc)
!
    relaName = compor(RELA_NAME)
    if (relaName(1:5) .eq. 'THER_') then
        call rccoma(zi(jvMaterc), 'THER', 1, therKeyword, propCode(1))
        aniso = ASTER_FALSE
        if (therKeyword(1:12) .eq. 'THER_NL_ORTH') then
            aniso = ASTER_TRUE
        end if
        call ntfcma(relaName, zi(jvMaterc), aniso, ifon)
        if (relaName(1:9) .eq. 'THER_HYDR') then
            call jevech('PHYDRPR', 'L', vr=hydrgp)
            call rcvalb('FPG1', 1, 1, '+', &
                        zi(jvMaterc), ' ', 'THER_HYDR', &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, propCode, 1)
        end if
    end if
!
    valQP = 0.0
    do kpg = 1, FEQuadCell%nbQuadPoints
        tpgi = FEEvalFuncRScal(FEBasis, temper, FEQuadCell%points_param(1:3, kpg))
!
        if (relaName(1:5) .eq. 'THER_') then
            call rcfode(ifon(1), tpgi, valQP(kpg), r8bid)
            if (relaName(1:9) .eq. 'THER_HYDR') then
                valQP(kpg) = valQP(kpg)-propVale(1)*hydrgp(kpg)
            end if
        else if (relaName(1:5) .eq. 'SECH_') then
            valQP(kpg) = tpgi
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
    call FeMakeRhsScal(FEQuadCell, FEBasis, valQP, resi)
!
    call writeVector("PRESIDU", nbDof, resi)
!
end subroutine
