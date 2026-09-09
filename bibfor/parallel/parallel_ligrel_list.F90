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
subroutine parallel_ligrel_list(numeEquZ, base)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "asterc/asmpi_allgather_i.h"
#include "asterc/asmpi_allgatherv_char24.h"
#include "asterc/asmpi_comm.h"
!
    character(len=*), intent(in) :: numeEquZ
    character(len=1), intent(in) :: base
!
! --------------------------------------------------------------------------------------------------
!
! Le but de cette routine est de classer dans le meme ordre les ligrels de charge
! en utilisant leur identifiant unique.
! De ce fait, en communicant la numerotation, on aura les bons raccords en vis-a-vis
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iLigr, iret, nbLigrTot, iProc, nbLigr, shift
    integer(kind=8) :: nbProc, rank, count, iNume, hashNb, iLigrT, ier
    integer(kind=8), allocatable :: v_nbLigr(:)
    integer(kind=8), pointer :: v_lilt(:) => null()
    character(len=8) :: typeLagr, typeLagrC, mesh, model
    character(len=19) :: numeEqua
    character(len=24) :: modeLoc, idenRela, hash
    character(len=19) :: ligrelName, joints
    character(len=24), pointer :: v_hash(:) => null()
    character(len=24), pointer :: v_refn(:) => null()
    character(len=24), pointer :: v_hash_list(:) => null()
    character(len=24), allocatable :: v_recv(:)
    aster_logical ::  lParallelMesh
    mpi_int :: mpicou, countSend, mrank, msize
    mpi_int, parameter :: mpi_one = to_mpi_int(1)
    mpi_int, allocatable :: v_count(:)
    mpi_int, allocatable :: v_displ(:)
!
! --------------------------------------------------------------------------------------------------
!
    numeEqua = numeEquZ

    call jeveuo(numeEqua//'.REFN', 'L', vk24=v_refn)
    mesh = v_refn(1)
    model = v_refn(3)
    call jelira(numeEqua//'.PRNO', 'NMAXOC', ival=nbLigr)
    lParallelMesh = isParallelMesh(mesh)
    if (.not. lParallelMesh) then
        call wkvect(numeEqua//'.LILT', base//' V I', nbLigr, vi=v_lilt)
        do iLigr = 1, nbLigr
            v_lilt(iLigr) = iLigr
        end do
    else
        call asmpi_comm('GET', mpicou)
        call asmpi_info(rank=mrank, size=msize)
        rank = to_aster_int(mrank)
        nbProc = to_aster_int(msize)

        if (nbLigr > 1) then
            call jecreo('&&TMP.HASHTABLE', 'V N K24')
            call jeecra('&&TMP.HASHTABLE', 'NOMMAX', nbLigr-1)
            call wkvect('&&TMP.HASHLIST', 'V V K24', nbLigr-1, vk24=v_hash_list)
            shift = 1
            do iLigr = 2, nbLigr
                call jenuno(jexnum(numeEqua//'.LILI', iLigr), ligrelName)
                call dismoi('JOINTS', ligrelName, 'LIGREL', repk=joints, arret='F')
!               Dans le cas ou le premier ligrel est &MAILLA et le deuxieme est le ligrel de modele
!               on doit oublier 2 ligrels pour obtenir des ligrels de charges
                if (ligrelName(1:8) .eq. model) then
                    if (iLigr .eq. 2) then
                        shift = 2
                        cycle
                    else
                        ASSERT(.false.)
                    end if
                end if
                hash = joints//'.HASH'
                call jeexin(hash, ier)
                if (ier .ne. 0) then
                    call jeveuo(hash, 'L', vk24=v_hash)
                    call jecroc(jexnom('&&TMP.HASHTABLE', v_hash(1)))
                    v_hash_list(iLigr-shift) = v_hash(1)
                else
!                   Si joints//'.HASH' est absent, on est dans le cas d'une charge séquentielle
!                   Son nom suffira a l'identifier car elle ne communiquera pas
                    v_hash_list(iLigr-shift) = ligrelName
                end if
            end do

            allocate (v_nbLigr(nbProc))
            allocate (v_displ(nbProc+1))
            allocate (v_count(nbProc))
            call asmpi_allgather_i([nbLigr-shift], mpi_one, v_nbLigr, mpi_one, mpicou)
            v_displ = 0
            count = 0
            do iProc = 1, nbProc
                count = count+v_nbLigr(iProc)
                v_count(iProc) = to_mpi_int(v_nbLigr(iProc))
                v_displ(iProc+1) = to_mpi_int(count)
            end do
            allocate (v_recv(count))
            countSend = to_mpi_int(nbLigr-shift)
            call asmpi_allgatherv_char24(v_hash_list, countSend, v_recv, v_count, v_displ, &
                                         mpicou)
            deallocate (v_nbLigr)
            deallocate (v_count)
            deallocate (v_displ)

            if (count .ne. 0) then
                call jecreo('&&TMP.HASHTABLETOT', 'V N K24')
                call jeecra('&&TMP.HASHTABLETOT', 'NOMMAX', count)
                do iLigr = 1, count
                    call jenonu(jexnom('&&TMP.HASHTABLETOT', v_recv(iLigr)), iNume)
                    if (iNume .eq. 0) then
                        call jecroc(jexnom('&&TMP.HASHTABLETOT', v_recv(iLigr)))
                    end if
                end do
                call jelira('&&TMP.HASHTABLETOT', 'NOMUTI', ival=hashNb)
            else
                hashNb = 0
            end if
            deallocate (v_recv)

            if (shift .eq. 2 .and. (hashNb+shift) .eq. 2) then
                call wkvect(numeEqua//'.LILT', base//' V I', 2, vi=v_lilt)
                v_lilt(1) = 1
                v_lilt(2) = 2
            else
                call wkvect(numeEqua//'.LILT', base//' V I', hashNb+shift, vi=v_lilt)
                v_lilt(1) = 1
                call jenuno(jexnum(numeEqua//'.LILI', 1), ligrelName)
                if (shift .eq. 2) then
                    v_lilt(2) = 2
                end if
                do iLigr = 1, hashNb
                    call jenuno(jexnum('&&TMP.HASHTABLETOT', iLigr), hash)
                    call jenonu(jexnom('&&TMP.HASHTABLE', hash), iNume)
                    if (iNume .ne. 0) then
                        call jenuno(jexnum(numeEqua//'.LILI', iNume+shift), ligrelName)
                        v_lilt(iLigr+shift) = iNume+shift
                    else
                        v_lilt(iLigr+shift) = -1
                    end if
                end do
                call jedetr('&&TMP.HASHTABLETOT')
            end if
            call jedetr('&&TMP.HASHLIST')
            call jedetr('&&TMP.HASHTABLE')
        else
            call wkvect(numeEqua//'.LILT', base//' V I', 1, vi=v_lilt)
            v_lilt(1) = 1
        end if
    end if
!
end subroutine
