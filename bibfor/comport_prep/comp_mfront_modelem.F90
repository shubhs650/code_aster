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
subroutine comp_mfront_modelem(elemTypeName, l_mfront_cp, &
                               modelMGIS, cplaMGIS, codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/teattr.h"
!
    character(len=16), intent(in) :: elemTypeName
    aster_logical, intent(in) :: l_mfront_cp
    integer(kind=8), intent(out) :: modelMGIS
    character(len=16), intent(out) :: cplaMGIS
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Select type of modelisation for MFront - On selected element
!
! --------------------------------------------------------------------------------------------------
!
! In  elemTypeName     : type of finite element
! In  l_mfront_cp      : .true. if plane stress is possible for this MFront behaviour
! Out modelMGIS        : finite element support for MFront
! Out cplaMGIS         : stress plane hypothesis (for Deborst)
! Out codret           : code for error
!                        0 - OK
!                        1 - Error - Not same finite element
!                        2 - Error - No MFront modelisation allowed on this element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=16) :: principal, model_type
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    modelMGIS = MGIS_MODEL_UNSET
    cplaMGIS = 'VIDE'

! - Get attributes on finite element
    call teattr('C', 'TYPMOD', model_type, iret, typel=elemTypeName)
    call teattr('C', 'PRINCIPAL', principal, iret, typel=elemTypeName)

! - Select modelisation for MFront
    if (principal .eq. 'OUI') then
        if (model_type .eq. '3D') then
            modelMGIS = MGIS_MODEL_TRIDIMENSIONAL
        elseif (model_type .eq. 'C_PLAN') then
            if (l_mfront_cp) then
                modelMGIS = MGIS_MODEL_PLANESTRESS
                cplaMGIS = 'ANALYTIQUE'
            else
                modelMGIS = MGIS_MODEL_AXISYMMETRICAL
                cplaMGIS = 'DEBORST'
            end if
        elseif (model_type .eq. 'D_PLAN') then
            modelMGIS = MGIS_MODEL_PLANESTRAIN
        elseif (model_type .eq. 'PLAN') then
            modelMGIS = MGIS_MODEL_PLANESTRAIN
        elseif (model_type .eq. 'AXIS') then
            modelMGIS = MGIS_MODEL_AXISYMMETRICAL
        elseif (model_type .eq. '1D') then
            modelMGIS = MGIS_MODEL_AXISYMMETRICAL
            cplaMGIS = 'DEBORST'
        else
            codret = 2
        end if
    end if
!
end subroutine
