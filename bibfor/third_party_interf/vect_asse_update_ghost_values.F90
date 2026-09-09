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
subroutine vect_asse_update_ghost_values(vasse, nume_equa)
#include "asterf_types.h"
    implicit none
#include "asterf_config.h"
#include "asterf_debug.h"
#include "asterf.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vector_update_ghost_values.h"
#include "jeveux.h"
!
    character(len=19), intent(inout) :: vasse
    character(len=19), intent(in) :: nume_equa
#if defined(ASTER_HAVE_MPI)
!
    real(kind=8), pointer :: vale(:) => null()
!
    character(len=8) :: noma
!
!----------------------------------------------------------------------
!
!   Communicates the values of the ghosts DOFs on a FieldOnNodes
!
!----------------------------------------------------------------------
!
    call jemarq()

    call dismoi('NOM_MAILLA', vasse, 'CHAMP', repk=noma)
    call jeveuo(vasse//'.VALE', 'E', vr=vale)
!
    call vector_update_ghost_values(vale, nume_equa, "SEND")
!
    call jedema()
#endif
!
end subroutine
