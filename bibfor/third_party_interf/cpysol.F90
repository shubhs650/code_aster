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
subroutine cpysol(nomat, numddl, rsolu, debglo, vecpet)
    use, intrinsic :: iso_c_binding
    implicit none
#include "asterf_config.h"
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterf_petsc.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/crnustd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrconl.h"
#include "asterfort/vector_update_ghost_values.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
#ifdef ASTER_HAVE_PETSC
    PetscInt :: debglo
#else
    integer(kind=4) :: debglo
#endif
    real(kind=8) :: rsolu(*), vecpet(*)
    character(len=14) :: numddl
    character(len=19) :: nomat
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
!
    integer(kind=8) :: rang, nbproc, lmat, iaux, nloc, iret, numglo
    integer(kind=8) :: nuno1, nucmp1, step
    aster_logical :: ldebug
    integer(kind=8), pointer :: v_nuls(:) => null()
    integer(kind=8), pointer :: v_deeg(:) => null()
    integer(kind=8), pointer :: nequ(:) => null()
    integer(kind=8), pointer :: pddl(:) => null()
    integer(kind=8), pointer :: nulg(:) => null()
    integer(kind=8), save :: nstep = 0
!
    mpi_int :: mrank, msize
!
    character(len=8) :: k8bid
    character(len=19) :: nomlig, comm_name, tag_name, joints, nume_equa
    character(len=24) :: domj, recv, send, gcom, pgid
    character(len=32) :: nojoine, nojoinr
!----------------------------------------------------------------------
!
    call jemarq()
!
    step = -1
    ldebug = ASTER_FALSE .and. step == nstep
    nstep = nstep+1
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    DEBUG_MPI('cpysol', rang, nbproc)
!
    nume_equa = numddl//".NUME"
    call jeveuo(nume_equa//'.NULG', 'L', vi=nulg)
    call jeveuo(nume_equa//'.PDDL', 'L', vi=pddl)
    call jeveuo(nume_equa//'.NEQU', 'L', vi=nequ)
    nloc = nequ(1)
!
    do iaux = 1, nloc
        if (pddl(iaux) .eq. rang) then
            numglo = nulg(iaux)
            rsolu(iaux) = vecpet(numglo-debglo+1)
        end if
    end do
!
    call vector_update_ghost_values(rsolu, nume_equa, "SEND")
!
! -- REMISE A L'ECHELLE DES LAGRANGES DANS LA SOLUTION
    call jeveuo(nomat//'.&INT', 'L', lmat)
    call mrconl('MULT', lmat, 0_8, 'R', rsolu, 1_8)
!
! -- debug
    if (ldebug) then
        print *, "DEBUG IN CPYSOL"
        call jeexin(nume_equa//'.NULS', iret)
        if (iret == 0) then
            call crnustd(numddl)
        end if
        call jeveuo(nume_equa//'.NULS', 'L', vi=v_nuls)
        call jeveuo(nume_equa//'.DEEG', 'L', vi=v_deeg)
        do iaux = 1, nloc
            if (pddl(iaux) .eq. rang) then
                nuno1 = v_deeg(2*(iaux-1)+1)
                nucmp1 = v_deeg(2*(iaux-1)+2)
                write (30+rang, *) nuno1, nucmp1, v_nuls(iaux), rsolu(iaux)
                !,zi(jnulg - 1 + iaux)
            end if
        end do
        flush (30+rang)
    end if

    call jedema()
#endif
!
end subroutine
