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
subroutine comp_read_typmodElem(elemTypeNume, l_mfront_cp, modelMGIS)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/comp_mfront_modelem.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in):: elemTypeNume
    aster_logical, intent(in) :: l_mfront_cp
    integer(kind=8), intent(out) :: modelMGIS
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Find dimension and type of modelisation for MFront
!
! --------------------------------------------------------------------------------------------------
!
! In  elemTypeNume     : type of finite element
! In  l_mfront_cp      : .true. if analytical plane stress
! Out modelMGIS        : type of modelisation MFront
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: codret
    character(len=16) :: elemTypeName, cplaMGIS
!
! --------------------------------------------------------------------------------------------------
!
    modelMGIS = MGIS_MODEL_UNSET
    cplaMGIS = 'VIDE'
    call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
    call comp_mfront_modelem(elemTypeName, l_mfront_cp, &
                             modelMGIS, cplaMGIS, codret)
    if (codret .eq. 2) then
        call utmess('F', 'COMPOR4_14', si=modelMGIS, &
                    sk="MGISBehaviourFort.h")
    end if
    ASSERT(codret .eq. 0)
!
end subroutine
