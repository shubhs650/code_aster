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
! ==================================================================================================
!
! Module for management of material of plates
!
! ==================================================================================================
!
module plateMaterial_module
! ==================================================================================================
    use MaterialPara_type
    use plate_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: chckMultiLayer
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! chckMultiLayer
!
! Check properties of multi-layered plates
!
! --------------------------------------------------------------------------------------------------
    subroutine chckMultiLayer(materPara, plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(in) :: materPara
        type(plateCara_Para), intent(in) :: plateCara
! ----- Local
        character(len=16), parameter :: elasKeyword = 'ELAS_COQMU'
        character(len=8), parameter :: fami = 'FPG1'
        character(len=1), parameter :: poum = '+'
        integer(kind=8), parameter :: kpg = 1, ksp = 1
        real(kind=8), parameter :: r8bid = 0.d0
        integer(kind=8) :: nbLayer, iLayer, iret, jvNbspIn
        real(kind=8) :: epTotaSum, epTota, epLayer
        character(len=3) :: iLayerStr
        integer(kind=8), parameter :: nbProp = 1
        character(len=16) :: propName(1)
        integer(kind=8) :: propCode(1)
        real(kind=8) :: propVale(1)
!   ------------------------------------------------------------------------------------------------
!
        call tecach('NNO', 'PNBSP_I', 'L', iret, iad=jvNbspIn)
        ASSERT(iret .eq. 0)
        iLayer = 0
        epTotaSum = 0.d0
        epTota = plateCara%thick
        nbLayer = plateCara%nbLayer
        ASSERT(materPara%elasID .eq. ELAS_COMPOSITE)

5       continue
        iLayer = iLayer+1

! ----- Get thickness of current layer
        call codent(iLayer, 'G', iLayerStr)
        propName(1) = 'C'//iLayerStr//'_V1'
        call rcvalb(fami, kpg, ksp, poum, &
                    materPara%jvMaterCode, ' ', elasKeyword, &
                    0, ' ', [r8bid], &
                    nbProp, propName, propVale, &
                    propCode, 0)
        if (propCode(1) .eq. 0) then
            epLayer = propVale(1)
            epTotaSum = epTotaSum+epLayer
            goto 5
        end if

! ----- Check total thickness
        if ((iLayer-1) .ne. nbLayer) then
            call utmess('F', 'PLATE1_51', ni=2, vali=[iLayer-1, nbLayer])
        end if
        if (abs(epTota-epTotaSum)/epTota .gt. 1.d-2) then
            call utmess('F', 'PLATE1_52', nr=2, valr=[epTotaSum, epTota])
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module plateMaterial_module
