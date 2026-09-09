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
subroutine te0243(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_eval_module
    use coorSyst_module, only: hasOrieField
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/ntcomp.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcdiff.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/writeVector.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: THER_*
!
! Options: RIGI_THER_TANG
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
    integer(kind=8) :: jvCamass, nbProp
    parameter(nbProp=3)
    integer(kind=8) :: propCode(nbProp)
    character(len=32) :: therKeyword
    real(kind=8) :: tpg, dtpg(3), diff, fluglo(3), Kglo(3, 3)
    real(kind=8) :: sechpg, dsechpg(3)
    real(kind=8) :: resi(MAX_BS_CG), rigi(MAX_BS_CG, MAX_BS_CG), dfluxglo(3), resi_p(MAX_BS_CG)
    real(kind=8) :: BGSEval(3, MAX_BS_CG), BSEval(MAX_BS_CG), eps, tempSave, delta
    real(kind=8), pointer :: flux(:) => null()
    real(kind=8), pointer :: temper(:) => null()

    integer(kind=8) ::  kpg, ifon(6)
    integer(kind=8) ::  jvMaterc, j, i_dof, iret
    character(len=16) :: relaName
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: aniso, l_rhs, l_diff
!
! --------------------------------------------------------------------------------------------------
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
!
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PTEMPEI', 'L', vr=temper)
    call jevech('PCOMPOR', 'L', vk16=compor)

    relaName = compor(RELA_NAME)
    l_rhs = option == "RAPH_THER"
    l_diff = .false.
    eps = 1.e-8

    if (l_rhs) then
        call jevech('PFLUXPR', 'E', vr=flux)
    end if
!
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
    resi = 0.d0
    rigi = 0.d0
    do kpg = 1, FEQuadCell%nbQuadPoints
        BSEval = FEBasis%func(FEQuadCell%points_param(1:3, kpg))
        BGSEval = FEBasis%grad(FEQuadCell%points_param(1:3, kpg), FEQuadCell%jacob(1:3, 1:3, kpg))
!
        if (relaName(1:5) .eq. 'THER_') then
            tpg = FEEvalFuncRScal(FEBasis, temper, FEQuadCell%points_param(1:3, kpg))
            dtpg = FEEvalGradVec(FEBasis, temper, FEQuadCell%points_param(1:3, kpg), BGSEval)
            call ntcomp(relaName, FECell%ndim, tpg, dtpg, &
                        FEQuadCell%points(1:3, kpg), aniso, ifon, fluglo, Kglo, dfluxglo)
            if (l_rhs) then
                flux(FECell%ndim*(kpg-1)+1:FECell%ndim*(kpg-1)+FECell%ndim) = -fluglo(1:FECell%ndim)
            end if
        else if (relaName(1:5) .eq. 'SECH_') then
            sechpg = FEEvalFuncRScal(FEBasis, temper, FEQuadCell%points_param(1:3, kpg))
            dsechpg = FEEvalGradVec(FEBasis, temper, FEQuadCell%points_param(1:3, kpg), BGSEval)
            call rcvarc(' ', 'TEMP', '+', 'RIGI', kpg, 1, tpg, iret)
            if (iret .ne. 0) call utmess('F', 'THERMIQUE1_2')
            call rcdiff(zi(jvMaterc), relaName, tpg, sechpg, diff)
            fluglo = diff*dsechpg
            Kglo = 0.d0
            do j = 1, FECell%ndim
                Kglo(j, j) = diff
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
        if (l_rhs .or. l_diff) then
            call FEStiffResiScalAdd(FEBasis, BGSEval, FEQuadCell%weights(kpg), fluglo, resi)
        end if
        ! if (.not. l_diff) then
        if (.not. l_rhs) then
            call FEStiffJacoScalAdd(FEBasis, BGSEval, FEQuadCell%weights(kpg), Kglo, rigi)
            if (relaName(1:5) .eq. 'THER_') then
                call FEMassStiffJacoScalAdd(BSEval, BGSEval, FEQuadCell%weights(kpg), &
                                            dfluxglo, rigi)
            end if
        end if
    end do

    if (l_diff .and. .not. l_rhs) then
        rigi = 0.d0
        do i_dof = 1, FEBasis%size
            ! save value
            tempSave = temper(i_dof)
            if (abs(temper(i_dof)) .lt. eps) then
                temper(i_dof) = tempSave+eps
            else
                temper(i_dof) = (1.d0+eps)*tempSave
            end if
            delta = temper(i_dof)-tempSave
            resi_p = 0.d0
            do kpg = 1, FEQuadCell%nbQuadPoints
                BGSEval = FEBasis%grad(FEQuadCell%points_param(1:3, kpg), &
                                       FEQuadCell%jacob(1:3, 1:3, kpg))
                !
                if (relaName(1:5) .eq. 'THER_') then
                    tpg = FEEvalFuncRScal(FEBasis, temper, FEQuadCell%points_param(1:3, kpg))
                    dtpg = &
                        FEEvalGradVec(FEBasis, temper, FEQuadCell%points_param(1:3, kpg), BGSEval)
                    call ntcomp(relaName, FECell%ndim, tpg, dtpg, &
                                FEQuadCell%points(1:3, kpg), aniso, ifon, fluglo, Kglo, dfluxglo)
                else if (relaName(1:5) .eq. 'SECH_') then
                    sechpg = FEEvalFuncRScal(FEBasis, temper, FEQuadCell%points_param(1:3, kpg))
                    dsechpg = FEEvalGradVec(FEBasis, temper, &
                                            FEQuadCell%points_param(1:3, kpg), BGSEval)
                    call rcvarc(' ', 'TEMP', '+', 'RIGI', kpg, 1, tpg, iret)
                    if (iret .ne. 0) call utmess('F', 'THERMIQUE1_2')
                    call rcdiff(zi(jvMaterc), relaName, tpg, sechpg, diff)
                    fluglo = diff*dsechpg
                else
                    ASSERT(ASTER_FALSE)
                end if
                call FEStiffResiScalAdd(FEBasis, BGSEval, FEQuadCell%weights(kpg), fluglo, resi_p)
            end do
            rigi(:, i_dof) = (resi_p-resi)/delta
            ! restore value
            temper(i_dof) = tempSave
        end do
        write (6, *) '--------------'
        do i_dof = 1, FEBasis%size
            write (6, *) 'rigi(', i_dof, ',:)=', rigi(i_dof, 1:FEBasis%size)
        end do
    end if
    !
    if (l_rhs) then
        call writeVector("PRESIDU", FEBasis%size, resi)
    else
        call writeMatrix("PMATTSR", FEBasis%size, FEBasis%size, ASTER_FALSE, rigi)
    end if
end subroutine
