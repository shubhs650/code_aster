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
subroutine getExternalStrainModel(defoComp, strainMGIS)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
!
    character(len=16), intent(in) :: defoComp
    integer(kind=8), intent(out) :: strainMGIS
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get model of strains for external programs (MFRONT)
!
! --------------------------------------------------------------------------------------------------
!
! In  defoComp         : model of strain (DEFORMATION keyword)
! Out strainMGIS       : model of (large) strains
!
! --------------------------------------------------------------------------------------------------
!
    strainMGIS = MGIS_STRAIN_UNSET

!   Obsolete - for trace
    ASSERT(defoComp .ne. 'GROT_GDEP')
    ASSERT(defoComp .ne. 'SIMO_MIEHE')
    if (defoComp .eq. 'PETIT' .or. &
        defoComp .eq. 'PETIT_REAC' .or. &
        defoComp .eq. 'GDEF_LOG') then
        strainMGIS = MGIS_STRAIN_SMALL
    else if (defoComp .eq. 'GREEN_LAGRANGE') then
        strainMGIS = MGIS_STRAIN_F
    else
        ASSERT(ASTER_FALSE)
    end if

end subroutine
