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
subroutine vector_comm_ghost_values_i(vector, mesh, mode)
#include "asterf_types.h"
    implicit none
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_sendrecv_i.h"
#include "asterf_config.h"
#include "asterf_debug.h"
#include "asterf.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/create_graph_comm.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    integer(kind=8), intent(inout) :: vector(*)
    character(len=8), intent(in) :: mesh
    character(len=*), intent(in) :: mode
#if defined(ASTER_HAVE_MPI)
!
    integer(kind=8) :: rang, nbproc, numpro, jjointr, jjointe
    integer(kind=8) :: lgenvo, lgrecep, jvaleue, jvaleur, iaux, jaux
    integer(kind=8) :: iret, numloc, numpr2, iret1, iret2
    integer(kind=8) :: nb_comm, domj_i
    integer(kind=8) :: jvaleue2, jvaleur2
    aster_logical :: l_parallel_mesh, lbidir
    integer(kind=8), pointer :: v_comm(:) => null()
    integer(kind=8), pointer :: v_tag(:) => null()
    integer(kind=8), pointer :: v_dom(:) => null()
    integer(kind=8), pointer :: v_gco(:) => null()
    integer(kind=4), pointer :: v_pgid(:) => null()
!
    mpi_int :: n4r, n4e, tag4, numpr4
    mpi_int :: mrank, msize, mpicou
!
    character(len=8) :: k8bid
    character(len=19) :: comm_name, tag_name, joints
    character(len=24) :: domj, recv, send, gcom, pgid
    character(len=32) :: nojoine, nojoinr
!
!----------------------------------------------------------------------
!
!   Communicates the values of ghost nodes (of a mesh)
!    2 modes:
!       - "SEND": one directionnal send from owner procs to others
!       - "BIDIR": bidirectionnal send (values ​​are summed)
!
!----------------------------------------------------------------------
!
    call jemarq()

    l_parallel_mesh = isParallelMesh(mesh)
    ASSERT(l_parallel_mesh)
!
    call asmpi_comm('GET', mpicou)
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    DEBUG_MPI('vector_comm_ghost_values_i', rang, nbproc)
!
    if (mode .eq. "SEND") then
        lbidir = .false.
    else if (mode .eq. "BIDIR") then
        lbidir = .true.
    else
        ASSERT(.false.)
    end if
!
!   -- Build the comm grpah
    comm_name = '&&CPYSOL.COMM'
    tag_name = '&&CPYSOL.TAG'
    joints = mesh//'.JOIN'
    domj = joints//".DOMJ"
    send = joints//".SEND"
    recv = joints//".RECV"
    gcom = joints//".GCOM"
    pgid = joints//".PGID"
    call create_graph_comm(mesh, "MAILLAGE_P", nb_comm, comm_name, tag_name)
    call jeveuo(comm_name, 'L', vi=v_comm)
    call jeveuo(tag_name, 'L', vi=v_tag)
    if (nb_comm > 0) then
        call jeveuo(domj, 'L', vi=v_dom)
        call jeveuo(gcom, 'L', vi=v_gco)
        call jeveuo(pgid, 'L', vi4=v_pgid)
        mpicou = int(v_gco(1), 4)
    end if
!
    do iaux = 1, nb_comm
        domj_i = v_comm(iaux)
        numpro = v_dom(domj_i)
        numpr2 = v_pgid(numpro+1)
        nojoine = jexnum(send, domj_i)
        call jeexin(nojoine, iret1)
        nojoinr = jexnum(recv, domj_i)
        call jeexin(nojoinr, iret2)
        lgrecep = 0
        lgenvo = 0
        if ((iret1+iret2) .ne. 0) then
            if (iret1 .ne. 0) then
                nojoine = jexnum(send, domj_i)
                call jelira(nojoine, 'LONMAX', lgenvo, k8bid)
                lgenvo = lgenvo/2
            end if
            if (iret2 .ne. 0) then
                nojoinr = jexnum(recv, domj_i)
                call jelira(nojoinr, 'LONMAX', lgrecep, k8bid)
                lgrecep = lgrecep/2
            end if
            ASSERT((lgenvo+lgrecep) .gt. 0)
!
            call wkvect('&&CPYSOL.TMP1E', 'V V I', max(1_8, lgenvo), jvaleue)
            call wkvect('&&CPYSOL.TMP1R', 'V V I', max(1_8, lgrecep), jvaleur)
            if (lbidir) then
                call wkvect('&&CPYSOL.TMP2E', 'V V I', max(1_8, lgrecep), jvaleue2)
                call wkvect('&&CPYSOL.TMP2R', 'V V I', max(1_8, lgenvo), jvaleur2)
            end if

            if (lgenvo > 0) then
                call jeveuo(nojoine, 'L', jjointe)
                do jaux = 0, lgenvo-1
                    numloc = zi(jjointe+2*jaux)
                    zi(jvaleue+jaux) = vector(numloc)
                end do
            end if

            if (lgrecep > 0) then
                call jeveuo(nojoinr, 'L', jjointr)
            end if
            if (lbidir) then
                if (lgrecep > 0) then
                    do jaux = 0, lgrecep-1
                        numloc = zi(jjointr+2*jaux)
                        zi(jvaleue2+jaux) = vector(numloc)
                    end do
                end if
            end if
!
            n4e = to_mpi_int(lgenvo)
            n4r = to_mpi_int(lgrecep)
            tag4 = to_mpi_int(v_tag(iaux))
            numpr4 = to_mpi_int(numpr2)
            call asmpi_sendrecv_i(zi(jvaleue), n4e, numpr4, tag4, &
                                  zi(jvaleur), n4r, numpr4, tag4, mpicou)
            if (lbidir) then
                call asmpi_sendrecv_i(zi(jvaleue2), n4r, numpr4, tag4, &
                                      zi(jvaleur2), n4e, numpr4, tag4, mpicou)
            end if

            if (.not. lbidir) then
                if (lgrecep > 0) then
                    do jaux = 0, lgrecep-1
                        numloc = zi(jjointr+2*jaux)
                        vector(numloc) = zi(jvaleur+2*jaux)
                    end do
                end if
            else
                if (lgrecep > 0) then
                    do jaux = 0, lgrecep-1
                        numloc = zi(jjointr+2*jaux)
                        vector(numloc) = vector(numloc)+zi(jvaleur+jaux)
                    end do
                end if
                if (lgenvo > 0) then
                    do jaux = 0, lgenvo-1
                        numloc = zi(jjointe+2*jaux)
                        vector(numloc) = vector(numloc)+zi(jvaleur2+jaux)
                    end do
                end if
            end if
            call jedetr('&&CPYSOL.TMP1E')
            call jedetr('&&CPYSOL.TMP1R')
            if (lbidir) then
                call jedetr('&&CPYSOL.TMP2E')
                call jedetr('&&CPYSOL.TMP2R')
            end if
        end if
    end do
!
    call jedetr(comm_name)
    call jedetr(tag_name)

    call jedema()
#endif
!
end subroutine
