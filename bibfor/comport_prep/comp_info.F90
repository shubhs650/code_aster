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
subroutine comp_info(modelZ, compor)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/comp_meca_pvar.h"
#include "asterfort/dismoi.h"
#include "asterfort/imvari.h"
#include "asterfort/jedetc.h"
!
    character(len=*), intent(in) :: modelZ
    character(len=19), intent(in) :: compor
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Print informations about COMPORTEMENT keyword
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  compor           : name of <CARTE> COMPOR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: comporInfo = '&&NMDOCC.INFO'
    character(len=19) :: modelFED
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelFED)

! - Prepare informations about internal variables
    call comp_meca_pvar(modelFED, compor, comporInfo)

! - Print informations about internal variables
    call imvari(comporInfo)

! - Clean
    call jedetc('V', comporInfo, 1)
!
end subroutine
