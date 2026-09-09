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
subroutine vect_asse_from_petsc(vasse, nume_equa, vecpet, scaling, ilocal)
#include "asterf_types.h"
#include "asterf_petsc.h"
    use aster_petsc_module
    implicit none
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_recv_r.h"
#include "asterc/asmpi_send_r.h"
#include "asterc/asmpi_sendrecv_r.h"
#include "asterf_config.h"
#include "asterf_debug.h"
#include "asterf.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/create_graph_comm.h"
#include "asterfort/crnustd.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
#ifdef ASTER_HAVE_PETSC
    Vec, intent(in) :: vecpet
#else
    integer(kind=4), intent(in) :: vecpet
#endif
    character(len=19), intent(inout) :: vasse
    character(len=19), intent(in) :: nume_equa
    real(kind=8), intent(in) :: scaling
    integer(kind=8) :: ilocal
#if defined(ASTER_HAVE_MPI) && defined(ASTER_HAVE_PETSC)
!
    integer(kind=8) :: rang, nbproc, numpro, jjointr, jjointe
    integer(kind=8) :: lgenvo, lgrecep, jvaleue, jvaleur, iaux, jaux, jnulg
    integer(kind=8) :: jprddl, jnequ, nloc, nlili, ili, iret, ijoin
    integer(kind=8) :: numglo, nuno1, nucmp1, numloc, numpr2
    integer(kind=8) :: iret1, iret2, jjoine, nbnoee, idprn1, idprn2, nec
    integer(kind=8) :: jjoinr, jnujoi1, jnujoi2, nbnoer, nddll, neq
    integer(kind=8) :: numnoe, nb_comm, gd, ieq, domj_i, nbLigrTot, iLigrT
    aster_logical :: ldebug, l_parallel_mesh, l_local
    integer(kind=8), pointer :: v_nuls(:) => null()
    integer(kind=8), pointer :: v_deeg(:) => null()
    integer(kind=8), pointer :: v_comm(:) => null()
    integer(kind=8), pointer :: v_tag(:) => null()
    integer(kind=8), pointer :: v_dom(:) => null()
    integer(kind=8), pointer :: v_gco(:) => null()
    integer(kind=4), pointer :: v_pgid(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    integer(kind=8), pointer :: v_lilt(:) => null()
!
    mpi_int :: n4r, n4e, tag4, numpr4
    mpi_int :: mrank, msize, mpicou
!
    character(len=8) :: k8bid, noma
    character(len=19) :: nomlig, comm_name, tag_name, joints
    character(len=24) :: domj, recv, send, gcom, pgid
    character(len=32) :: nojoine, nojoinr
!
    PetscScalar, pointer :: xx(:) => null()
    PetscErrorCode ::  ierr
    VecScatter :: ctx
    PetscInt :: low, high
    Vec :: vecgth
!----------------------------------------------------------------------
!
!   Import a PETSc vector into a FieldOnNodes
!
!----------------------------------------------------------------------
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS PRNO DES S.D. LIGREL
!     REPERTORIEES DANS LE CHAMP LILI DE NUME_DDL ET A LEURS ADRESSES
!     ZZPRNO(ILI,NUNOEL,1) = NUMERO DE L'EQUATION ASSOCIEES AU 1ER DDL
!                            DU NOEUD NUNOEL DANS LA NUMEROTATION LOCALE
!                            AU LIGREL ILI DE .LILI
!     ZZPRNO(ILI,NUNOEL,2) = NOMBRE DE DDL PORTES PAR LE NOEUD NUNOEL
!     ZZPRNO(ILI,NUNOEL,2+1) = 1ER CODE
!     ZZPRNO(ILI,NUNOEL,2+NEC) = NEC IEME CODE
!
#define zzprno(ili, nunoel, l) zi(idprn1 - 1 + zi(idprn2 + ili - 1) + (nunoel - 1)*(nec + 2) + l-1)
!
    call jemarq()
!
    ldebug = ASTER_FALSE

    call dismoi('NOM_MAILLA', vasse, 'CHAMP', repk=noma)
    call jeveuo(vasse//'.VALE', 'E', vr=vale)

    l_parallel_mesh = isParallelMesh(noma)

!  is the mesh distributed across processors ?
    if (.not. l_parallel_mesh) then

!       -- Reconstruction of the solution on each process
        call VecScatterCreateToAll(vecpet, ctx, vecgth, ierr)
        ASSERT(ierr .eq. 0)
        call VecScatterBegin(ctx, vecpet, vecgth, INSERT_VALUES, SCATTER_FORWARD, &
                             ierr)
        ASSERT(ierr .eq. 0)
        call VecScatterEnd(ctx, vecpet, vecgth, INSERT_VALUES, SCATTER_FORWARD, &
                           ierr)
        ASSERT(ierr .eq. 0)
        call VecScatterDestroy(ctx, ierr)
        ASSERT(ierr .eq. 0)
!
!       -- Copy to the field
        call VecGetArray(vecgth, xx, ierr)
        ASSERT(ierr .eq. 0)
!
        call jelira(vasse//'.VALE', 'LONMAX', neq)
        do ieq = 1, neq
            vale(ieq) = xx(ieq)
        end do
!
        call VecRestoreArray(vecgth, xx, ierr)
        ASSERT(ierr .eq. 0)
!
!       -- Cleanup
        call VecDestroy(vecgth, ierr)
        ASSERT(ierr .eq. 0)

    else

!
        call asmpi_comm('GET', mpicou)
!
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
        DEBUG_MPI('vect_asse_from_petsc', rang, nbproc)
!
!
        call jeveuo(nume_equa//'.NULG', 'L', jnulg)
        call jeveuo(nume_equa//'.PDDL', 'L', jprddl)
        call jeveuo(nume_equa//'.NEQU', 'L', jnequ)
        call jeveuo(nume_equa//'.LILT', 'L', vi=v_lilt)
        call jelira(nume_equa//'.LILT', 'LONMAX', ival=nbLigrTot)
        nloc = zi(jnequ)

        call VecGetOwnershipRange(vecpet, low, high, ierr)
        ASSERT(ierr .eq. 0)
        ASSERT((ilocal == 1) .or. (ilocal == 0))
        l_local = (ilocal .eq. 1)
!
!       -- Retrieve the fortran array
        call VecGetArray(vecpet, xx, ierr)
        ASSERT(ierr .eq. 0)
!
        if (l_local) then
            do iaux = 0, nloc-1
                vale(iaux+1) = xx(iaux+1)
            end do
        else
            do iaux = 0, nloc-1
                if (zi(jprddl+iaux) .eq. rang) then
                    numglo = zi(jnulg+iaux)
                    vale(iaux+1) = xx(numglo-low+1)
                end if
            end do
!
!           -- Build the comm grpah
            comm_name = '&&CPYSOL.COMM'
            tag_name = '&&CPYSOL.TAG'
            call dismoi("JOINTS", nume_equa, "NUME_EQUA", repk=joints, arret="F")
            domj = joints//".DOMJ"
            send = joints//".SEND"
            recv = joints//".RECV"
            gcom = joints//".GCOM"
            pgid = joints//".PGID"
            call create_graph_comm(nume_equa, "NUME_EQUA", nb_comm, comm_name, tag_name)
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
                nojoinr = jexnum(recv, domj_i)
                nojoine = jexnum(send, domj_i)
                call jeexin(nojoine, iret1)
                call jeexin(nojoinr, iret2)
                lgrecep = 0
                lgenvo = 0
                if ((iret1+iret2) .ne. 0) then
                    if (iret1 .ne. 0) then
                        call jelira(nojoine, 'LONMAX', lgenvo, k8bid)
                    end if
                    if (iret2 .ne. 0) then
                        call jelira(nojoinr, 'LONMAX', lgrecep, k8bid)
                    end if
                    ASSERT((lgenvo+lgrecep) .gt. 0)
!
                    call wkvect('&&CPYSOL.TMP1E', 'V V R', max(1_8, lgenvo), jvaleue)
                    call wkvect('&&CPYSOL.TMP1R', 'V V R', max(1_8, lgrecep), jvaleur)

                    if (lgenvo > 0) then
                        call jeveuo(nojoine, 'L', jjointe)
                        do jaux = 0, lgenvo-1
                            numloc = zi(jjointe+jaux)
                            ASSERT(zi(jprddl+numloc-1) .eq. rang)
                            zr(jvaleue+jaux) = vale(numloc)
                        end do
                    end if
!
                    n4e = to_mpi_int(lgenvo)
                    n4r = to_mpi_int(lgrecep)
                    tag4 = to_mpi_int(v_tag(iaux))
                    numpr4 = to_mpi_int(numpr2)
                    call asmpi_sendrecv_r(zr(jvaleue), n4e, numpr4, tag4, &
                                          zr(jvaleur), n4r, numpr4, tag4, mpicou)

                    if (lgrecep > 0) then
                        call jeveuo(nojoinr, 'L', jjointr)
                        do jaux = 0, lgrecep-1
                            numloc = zi(jjointr+jaux)
                            ASSERT(zi(jprddl+numloc-1) .eq. numpro)
                            vale(numloc) = zr(jvaleur+jaux)
                        end do
                    end if
                    call jedetr('&&CPYSOL.TMP1E')
                    call jedetr('&&CPYSOL.TMP1R')
                end if
            end do
!
            call jedetr(comm_name)
            call jedetr(tag_name)

!   Retrieve .PRNO and .NUME adresses
            call jeveuo(nume_equa//'.PRNO', 'E', idprn1)
            call jeveuo(jexatr(nume_equa//'.PRNO', 'LONCUM'), 'L', idprn2)

!   !!! Check that no Super Elements exist in the model
            call dismoi('NUM_GD_SI', nume_equa, 'NUME_EQUA', repi=gd)
            nec = nbec(gd)
            call jelira(nume_equa//'.PRNO', 'NMAXOC', nlili, k8bid)
            do iLigrT = 2, nbLigrTot
                if (v_lilt(iLigrT) .eq. -1) then
                    cycle
                else
                    ili = v_lilt(iLigrT)
                end if
                call jenuno(jexnum(nume_equa//'.LILI', ili), nomlig)
                call create_graph_comm(nomlig, "LIGREL", nb_comm, comm_name, tag_name)
                if (nb_comm > 0) then
                    call jeveuo(comm_name, 'L', vi=v_comm)
                    call jeveuo(tag_name, 'L', vi=v_tag)
                    call dismoi("JOINTS", nomlig, "LIGREL", repk=joints, arret="F")
                    domj = joints//".DOMJ"
                    send = joints//".SEND"
                    recv = joints//".RECV"
                    gcom = joints//".GCOM"
                    pgid = joints//".PGID"
                    call jeveuo(gcom, 'L', vi=v_gco)
                    call jeveuo(pgid, 'L', vi4=v_pgid)
                    mpicou = int(v_gco(1), 4)
                    call jeveuo(domj, 'L', vi=v_dom)
                    do ijoin = 1, nb_comm
                        domj_i = v_comm(ijoin)
                        numpro = v_dom(domj_i)
                        numpr2 = v_pgid(numpro+1)
                        numpr4 = to_mpi_int(numpr2)
                        tag4 = to_mpi_int(v_tag(ijoin))
                        nojoine = jexnum(send, domj_i)
                        nojoinr = jexnum(recv, domj_i)

                        call jeexin(nojoine, iret1)
                        if (iret1 .ne. 0) then
                            call jeveuo(nojoine, 'L', jjoine)
                            call jelira(nojoine, 'LONMAX', nbnoee, k8bid)
                            call wkvect('&&CRNUGL.NUM_DDL_GLOB_E', 'V V R', nbnoee, jnujoi1)
                            do jaux = 1, nbnoee
                                numnoe = -zi(jjoine+jaux-1)
                                nddll = zzprno(ili, numnoe, 1)
                                zr(jnujoi1+jaux-1) = vale(nddll)
                            end do
                            n4e = to_mpi_int(nbnoee)
                        end if

                        call jeexin(nojoinr, iret2)
                        if (iret2 .ne. 0) then
                            call jeveuo(nojoinr, 'L', jjoinr)
                            call jelira(nojoinr, 'LONMAX', nbnoer, k8bid)
                            call wkvect('&&CRNUGL.NUM_DDL_GLOB_R', 'V V R', nbnoer, jnujoi2)
                            n4r = to_mpi_int(nbnoer)
                        end if

                        if (rang .lt. numpro) then
                            if (iret1 .ne. 0) then
                                call asmpi_send_r(zr(jnujoi1), n4e, numpr4, tag4, mpicou)
                            end if
                            if (iret2 .ne. 0) then
                                call asmpi_recv_r(zr(jnujoi2), n4r, numpr4, tag4, mpicou)
                            end if
                        else if (rang .gt. numpro) then
                            if (iret2 .ne. 0) then
                                call asmpi_recv_r(zr(jnujoi2), n4r, numpr4, tag4, mpicou)
                            end if
                            if (iret1 .ne. 0) then
                                call asmpi_send_r(zr(jnujoi1), n4e, numpr4, tag4, mpicou)
                            end if
                        end if

                        if (iret2 .ne. 0) then
                            do jaux = 1, nbnoer
                                numnoe = -zi(jjoinr+jaux-1)
                                nddll = zzprno(ili, numnoe, 1)
                                vale(nddll) = zr(jnujoi2+jaux-1)
                            end do
                        end if
                        call jedetr('&&CRNUGL.NUM_DDL_GLOB_E')
                        call jedetr('&&CRNUGL.NUM_DDL_GLOB_R')
                    end do
                end if
                call jedetr(comm_name)
                call jedetr(tag_name)
            end do
        end if
!
!
        call VecRestoreArray(vecpet, xx, ierr)
        ASSERT(ierr .eq. 0)

!
! -- debug
        if (ldebug) then
            print *, "DEBUG IN vect_asse_from_petsc"
            call jeexin(nume_equa//'.NULS', iret)
            if (iret == 0) then
                call crnustd(nume_equa)
            end if
            call jeveuo(nume_equa//'.NULS', 'L', vi=v_nuls)
            call jeveuo(nume_equa//'.DEEG', 'L', vi=v_deeg)
            do iaux = 1, nloc
                if (zi(jprddl+iaux-1) .eq. rang) then
                    nuno1 = v_deeg(2*(iaux-1)+1)
                    nucmp1 = v_deeg(2*(iaux-1)+2)
                    write (30+rang, *) nuno1, nucmp1, v_nuls(iaux), vale(iaux)
                    !,zi(jnulg - 1 + iaux)
                end if
            end do
            flush (30+rang)
        end if

    end if

!   Scale the Lagrange multipliers
    call jelira(vasse//'.VALE', 'LONMAX', neq)
    call jeveuo(nume_equa//'.DEEQ', 'L', vi=deeq)
    do ieq = 1, neq
        if (deeq(2*(ieq-1)+2) .le. 0) then
            vale(ieq) = scaling*vale(ieq)
        end if
    end do

    call jedema()
#endif
!
end subroutine
