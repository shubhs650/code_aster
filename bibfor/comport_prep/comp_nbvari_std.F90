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
subroutine comp_nbvari_std(relaComp, defoComp, typeCpla, &
                           kitComp, postIter, reguVisc, postIncr, &
                           nbVari, numeLaw)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lcinfo.h"
#include "asterc/lcdiscard.h"
#include "asterfort/comp_meca_code.h"
!
    character(len=16), intent(in) :: relaComp, defoComp, typeCpla
    character(len=16), intent(in) :: kitComp(4), postIter, reguVisc, postIncr
    integer(kind=8), intent(out) :: nbVari, numeLaw
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Get number of internal state variables for standard constitutive laws
!
! --------------------------------------------------------------------------------------------------
!
! In  relaComp         : behaviour (RELATION keyword)
! In  defoComp         : model of strain (DEFORMATION keyword)
! In  typeCpla         : plane stress method (analytical or De Borst algorithm)
! In  kitComp          : KIT behaviour
! In  postIter         : type of post_treatment at each Newton iteration (POST_ITER keyword)
! In  reguVisc         : keyword for viscuous regularization (REGU_VISC keyword)
! In  postIncr         : type of post-treatment at end of time step (POST_INCR keyword)
! Out nbVari           : number of internal state variables
! Out numeLaw          : index of subroutine for behaviour
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: compCodePY
    integer(kind=8) :: idummy
!
! --------------------------------------------------------------------------------------------------
!
    nbVari = 0
    numeLaw = 0

! - Coding composite comportment (Python)
    call comp_meca_code(relaComp, defoComp, typeCpla, kitComp, &
                        postIter, reguVisc, postIncr, &
                        compCodePY)

! - Get number of total internal state variables and index of law
    call lcinfo(compCodePY, numeLaw, nbVari, idummy)

! - End of encoding
    call lcdiscard(compCodePY)
!
end subroutine
