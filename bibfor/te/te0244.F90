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
subroutine te0244(option, nomte)
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
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/ntcomp.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcdiff.h"
#include "asterfort/rcfode.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
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
! Options: CHAR_THER_EVOLNI
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
    character(len=16), parameter :: propName(nbProp) = (/'CHALHYDR'/)
    real(kind=8) :: propVale(nbProp)
    character(len=16) :: therKeyword, relaName
    real(kind=8) :: valQPM(MAX_QP), BGSEval(3, MAX_BS_CG)
    real(kind=8) :: valQPMP(MAX_QP)
    real(kind=8) :: resi_f(MAX_BS_CG), resi_m(MAX_BS_CG), resi(MAX_BS_CG)
    real(kind=8) :: resi_mp(MAX_BS_CG), resi_p(MAX_BS_CG), dfluxglo(3)
    real(kind=8) ::  deltat, theta, chal(1), diff, Kglo(3, 3)
    real(kind=8) :: beta, dbeta, tpg, dtpg(3), flux(3), sechpg, dsechpg(3)
    integer(kind=8) :: kpg, jvMaterc, jvCamass, ifon(6), jvInstr, iret
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: lhyd, aniso
    real(kind=8), pointer :: temper(:) => null()
    real(kind=8), pointer :: hydrpg(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call FECell%init()
    call FEBasis%initCell(FECell)
    call FEQuadMass%initCell(FECell, "MASS")
    call FEQuadRigi%initCell(FECell, "RIGI")
!
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PINSTR', 'L', jvInstr)
    call jevech('PTEMPER', 'L', vr=temper)
    call jevech('PCOMPOR', 'L', vk16=compor)
    relaName = compor(RELA_NAME)
!
    deltat = zr(jvInstr+1)
    theta = zr(jvInstr+2)
    jvCamass = 0

    if (relaName(1:5) .eq. 'THER_') then
        call rccoma(zi(jvMaterc), 'THER', 1, therKeyword, propCode(1))
        aniso = ASTER_FALSE
        if (therKeyword(1:12) .eq. 'THER_NL_ORTH') then
            aniso = ASTER_TRUE
        end if
        call ntfcma(relaName, zi(jvMaterc), aniso, ifon)
        if (aniso) then
            if (.not. hasOrieField(jvCamass)) then
                call utmess('F', 'THERMIQUE1_3')
            end if
        end if
    end if
!
    resi_f = 0.d0
    do kpg = 1, FEQuadRigi%nbQuadPoints
        BGSEval = FEBasis%grad(FEQuadRigi%points_param(1:3, kpg), FEQuadRigi%jacob(1:3, 1:3, kpg))
!
        if (relaName(1:5) .eq. 'THER_') then
            tpg = FEEvalFuncRScal(FEBasis, temper, FEQuadRigi%points_param(1:3, kpg))
            dtpg = FEEvalGradVec(FEBasis, temper, FEQuadRigi%points_param(1:3, kpg), BGSEval)
            call ntcomp(relaName, FECell%ndim, tpg, dtpg, &
                        FEQuadRigi%points(1:3, kpg), aniso, ifon, flux, Kglo, dfluxglo)
        else if (relaName(1:5) .eq. 'SECH_') then
            sechpg = FEEvalFuncRScal(FEBasis, temper, FEQuadRigi%points_param(1:3, kpg))
            dsechpg = FEEvalGradVec(FEBasis, temper, FEQuadRigi%points_param(1:3, kpg), BGSEval)
            call rcvarc(' ', 'TEMP', '-', 'RIGI', kpg, 1, tpg, iret)
            if (iret .ne. 0) call utmess('F', 'THERMIQUE1_2')
            call rcdiff(zi(jvMaterc), relaName, tpg, sechpg, diff)
            flux = diff*dsechpg
        else
            ASSERT(ASTER_FALSE)
        end if
        call FEStiffResiScalAdd(FEBasis, BGSEval, FEQuadRigi%weights(kpg), flux, resi_f)
    end do
!
    if (relaName(1:9) .eq. 'THER_HYDR') then
        lhyd = ASTER_TRUE
        call jevech('PHYDRPM', 'L', vr=hydrpg)
        call rcvalb('FPG1', 1, 1, '+', &
                    zi(jvMaterc), ' ', 'THER_HYDR', &
                    0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 1)
    else
        lhyd = ASTER_FALSE
    end if
!
    do kpg = 1, FEQuadMass%nbQuadPoints
        tpg = FEEvalFuncRScal(FEBasis, temper, FEQuadMass%points_param(1:3, kpg))
        if (relaName(1:5) .eq. 'THER_') then
            call rcfode(ifon(1), tpg, beta, dbeta)
            if (lhyd) then
                valQPMP(kpg) = (dbeta*tpg-propVale(1)*hydrpg(kpg))
                valQPM(kpg) = (beta-propVale(1)*hydrpg(kpg))
            else
                valQPMP(kpg) = dbeta*tpg
                valQPM(kpg) = beta
            end if
        else if (relaName(1:5) .eq. 'SECH_') then
            valQPM(kpg) = tpg
            valQPMP(kpg) = tpg
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
    call FeMakeRhsScal(FEQuadMass, FEBasis, ValQPM, resi_m)
    call FeMakeRhsScal(FEQuadMass, FEBasis, ValQPMP, resi_mp)
!
    resi = (theta-1.0d0)*resi_f+resi_m/deltat
    resi_p = (theta-1.0d0)*resi_f+resi_mp/deltat
!
    call writeVector('PVECTTR', FEBasis%size, resi)
    call writeVector('PVECTTI', FEBasis%size, resi_p)
!
end subroutine
