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

subroutine vector_ghosts_comm(mesh, vect_size, v_vect)
    implicit none
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_sendrecv_i.h"
#include "asterf_config.h"
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/create_graph_comm.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(inout) :: vect_size
    integer(kind=8), intent(inout) :: v_vect(*)

#ifdef ASTER_HAVE_MPI
!
    real(kind=8) :: r_cmp_nb
    integer(kind=8) :: rang, nbproc, iaux, nddll, node_nb, value_nb, cmp_nb
    integer(kind=8) :: nb_comm, icmp
    integer(kind=8) :: numpro, jjoinr, jjoine, nbnoee, jaux, numno1, numno2
    integer(kind=8) :: jenvoi1, lgenve1, lgenvr1, poscom
    integer(kind=8) :: nddl, numpr2
    integer(kind=8) :: nbnoer, jrecep1, curpos
    integer(kind=8) :: nb_ddl_envoi, domj_i
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: jnujoi2, iret
    mpi_int :: mrank, mnbproc, mpicou, tag4, numpr4, n4e, n4r
    integer(kind=8), pointer :: v_comm(:) => null()
    integer(kind=8), pointer :: v_tag(:) => null()
    integer(kind=8), pointer :: v_dom(:) => null()
    integer(kind=8), pointer :: v_gco(:) => null()
    integer(kind=4), pointer :: v_pgid(:) => null()
!
    character(len=8) :: k8bid
    character(len=19) :: comm_name, tag_name, meshj
!
!----------------------------------------------------------------------
!
    call jemarq()
!
    call infniv(ifm, niv)
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(rank=mrank, size=mnbproc)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(mnbproc)
    DEBUG_MPI('vector_ghosts_comm', rang, nbproc)
!
! -- Création du graphe de comm
    meshj = mesh//".JOIN"
    comm_name = '&CRNULG.COMM'
    tag_name = '&CRNULG.TAG'
    call create_graph_comm(mesh, "MAILLAGE_P", nb_comm, comm_name, tag_name)
!
    if (nb_comm > 0) then
        call jeveuo(meshj//'.DOMJ', 'L', vi=v_dom)
        call jeveuo(meshj//'.GCOM', 'L', vi=v_gco)
        call jeveuo(meshj//'.PGID', 'L', vi4=v_pgid)
        call jeveuo(comm_name, 'L', vi=v_comm)
        call jeveuo(tag_name, 'L', vi=v_tag)
        mpicou = to_mpi_int(v_gco(1))
    end if
!
    ! call jeveuo(vector, 'L', vi=v_vect)
    ! call jelira(vector, 'LONMAX', value_nb, k8bid)
    value_nb = vect_size
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=node_nb)
    cmp_nb = value_nb/node_nb
    r_cmp_nb = value_nb/node_nb
    if (r_cmp_nb .ne. cmp_nb) then
        ASSERT(.false.)
    end if
!
!   Il faut maintenant communiquer les valeurs sur les noeuds partages
    do iaux = 1, nb_comm
        domj_i = v_comm(iaux)
        numpro = v_dom(domj_i)
        numpr2 = v_pgid(numpro+1)
        call jeveuo(jexnum(meshj//".SEND", domj_i), 'L', jjoine)
        call jelira(jexnum(meshj//".SEND", domj_i), 'LONMAX', nbnoee, k8bid)
        call jeveuo(jexnum(meshj//".RECV", domj_i), 'L', jjoinr)
        call jelira(jexnum(meshj//".RECV", domj_i), 'LONMAX', nbnoer, k8bid)
        nbnoee = nbnoee/2
        nbnoer = nbnoer/2
!
!       DES DEUX COTES LES NOEUDS NE SONT PAS DANS LE MEME ORDRE ?
        tag4 = to_mpi_int(v_tag(iaux))
        numpr4 = to_mpi_int(numpr2)
        lgenve1 = nbnoee*(1+cmp_nb)+1
        lgenvr1 = nbnoer*(1+cmp_nb)+1
        call wkvect('&&CRNULG.NOEUD_NEC_E1', 'V V I', lgenve1, jenvoi1)
        call wkvect('&&CRNULG.NOEUD_NEC_R1', 'V V I', lgenvr1, jrecep1)
!
!       On commence par envoyer, le but final est de recevoir les numeros de ddl
!       On boucle donc sur les noeuds a recevoir
        nb_ddl_envoi = 0
        do jaux = 1, nbnoee
            poscom = (jaux-1)*(1+cmp_nb)+1
            numno1 = zi(jjoine+2*(jaux-1))
            numno2 = zi(jjoine+2*jaux-1)
            zi(jenvoi1+poscom) = numno2
            do icmp = 1, cmp_nb
                zi(jenvoi1+poscom+icmp) = v_vect(cmp_nb*(numno1-1)+icmp)
            end do
            nb_ddl_envoi = nb_ddl_envoi+cmp_nb
        end do
        zi(jenvoi1) = nb_ddl_envoi
        n4e = to_mpi_int(lgenve1)
        n4r = to_mpi_int(lgenvr1)
        call asmpi_sendrecv_i(zi(jenvoi1), n4e, numpr4, tag4, &
                              zi(jrecep1), n4r, numpr4, tag4, mpicou)

        if (zi(jrecep1) > 0) then
!
            do jaux = 1, nbnoer
                poscom = (jaux-1)*(1+cmp_nb)+1
                numno1 = zi(jrecep1+poscom)
!
                do icmp = 1, cmp_nb
                    v_vect(cmp_nb*(numno1-1)+icmp) = zi(jrecep1+poscom+icmp)
                end do
            end do
        end if
!
        call jedetr('&&CRNULG.NOEUD_NEC_E1')
        call jedetr('&&CRNULG.NOEUD_NEC_R1')
    end do
    call jedetr(comm_name)
    call jedetr(tag_name)
!
    call jedema()
#else
    character(len=8) :: kbid
    kbid = mesh
#endif
!
end subroutine
