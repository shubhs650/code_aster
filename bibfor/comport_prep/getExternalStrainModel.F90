! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine getExternalStrainModel(defo_comp, strain_model)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: defo_comp
    integer, intent(out) :: strain_model
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get model of strains for external programs (MFRONT)
!
! --------------------------------------------------------------------------------------------------
!
! In  defo_comp        : value of DEFORMATION keyword
! Out strain_model     : model of (large) strains
!                        1 - small strains
!                        2 - Simo-Miehe
!                        3 - GreenLagrange
! --------------------------------------------------------------------------------------------------
    strain_model = MFRONT_STRAIN_UNSET

! ----- Indicator for large strains

!   Obsolete - for trace
    ASSERT(defo_comp .ne. 'GROT_GDEP')
    ASSERT(defo_comp .ne. 'SIMO_MIEHE')

!   for GDEF_LOG, prelog/poslog are called in te*
    if (defo_comp .eq. 'PETIT' .or. &
        defo_comp .eq. 'PETIT_REAC' .or. &
        defo_comp .eq. 'GDEF_LOG') then
        strain_model = MFRONT_STRAIN_SMALL
    else if ( defo_comp .eq. 'GREEN_LAGRANGE' ) then
        strain_model = MFRONT_STRAIN_F            
    else
        ASSERT(ASTER_FALSE)
    end if

end subroutine
