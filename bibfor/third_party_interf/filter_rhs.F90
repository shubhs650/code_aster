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

subroutine filter_rhs(rsolu, nume_equa)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vector_update_ghost_values.h"
    real(kind=8), intent(inout) :: rsolu(*)
    character(len=19), intent(in) :: nume_equa
!-----------------------------------------------------------------------
    integer(kind=8) :: rang, nloc, ieq, iexi
    mpi_int :: mrank
    integer(kind=8), dimension(:), pointer :: delg => null()
    integer(kind=8), dimension(:), pointer :: nequ => null()
    integer(kind=8), dimension(:), pointer :: pddl => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    call jeexin(nume_equa//'.PDDL', iexi)
    if (iexi .ne. 0) then
        call jeveuo(nume_equa//'.NEQU', 'L', vi=nequ)
        call jeveuo(nume_equa//'.DELG', 'L', vi=delg)
        call jeveuo(nume_equa//'.PDDL', 'L', vi=pddl)
        call asmpi_info(rank=mrank)
        rang = to_aster_int(mrank)
        nloc = nequ(1)
        do ieq = 1, nloc
            if (delg(ieq) .ge. 0d0 .and. pddl(ieq) .ne. rang) then
                rsolu(ieq) = 0.d0
            end if
        end do
        call vector_update_ghost_values(rsolu, nume_equa, "BIDIR")
    end if
!
    call jedema()
end subroutine
