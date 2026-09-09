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
subroutine vector_zero_late_values(vector, nume_equa)
#include "asterf_types.h"
    implicit none
#include "asterf_config.h"
#include "asterf.h"
#include "jeveux.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
    real(kind=8), intent(inout) :: vector(*)
    character(len=19), intent(in) :: nume_equa
#if defined(ASTER_HAVE_MPI)
    integer(kind=8) :: ieq, nddl, rang, nbproc
    mpi_int :: mrank, msize
    integer(kind=8), dimension(:), pointer :: delg => null()
    integer(kind=8), dimension(:), pointer :: pddl => null()
    integer(kind=8), dimension(:), pointer :: nequ => null()
!
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
    call jemarq()
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    call jeveuo(nume_equa//'.PDDL', 'L', vi=pddl)
    call jeveuo(nume_equa//'.NEQU', 'L', vi=nequ)
    nddl = nequ(1)
!
    call jeveuo(nume_equa//'.DELG', 'L', vi=delg)
    do ieq = 1, nddl
        if ((delg(ieq) .lt. 0) .and. (pddl(ieq) .ne. rang)) then
            vector(ieq) = 0.d0
        end if
    end do
!
    call jedema()
#endif
!
end subroutine
