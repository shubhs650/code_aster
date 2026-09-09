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
subroutine iremed_filtre(nomast, nomsd, base, par_seqfile)
!
    implicit none
!
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_sendrecv_i.h"
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/create_graph_comm.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/lrmtyp.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
#endif
!
    character(len=1) :: base
    character(len=8) :: nomast, nomsd
    aster_logical :: par_seqfile
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE DU MAILLAGE - FORMAT MED - FILTRES
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREE :
!       NOMAST : NOM DU MAILLAGE
!       NOMSD : NOM DE LA SD DE SORTIE
!       BASE : BASE DE CREATION DE LA SD
!       PAR_SEQFILE : BOOLEEN POUR PRECISER SI ON EST EN PARALLELE AVEC FICHIER MED UNIQUE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbtyp, nbnoeu, nbmail, iret
    integer(kind=8) :: nnotyp(MT_NTYMAX), typgeo(MT_NTYMAX)
    integer(kind=8) :: renumd(MT_NTYMAX), modnum(MT_NTYMAX), numnoa(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: nuanom(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: ino, jnoex, jnbno, jnbno1, nbnot, iproc, jno, jma, jnbma, nbmat, jmaex
    integer(kind=8) :: rang, nbproc, ityp, jtyp, jtyp2, iaux, ifm, niv, jtypg
    integer(kind=8) :: ima, nb_comm, i_comm, domj_i
    integer(kind=8) :: jjoine, jjoinr, nbnoee, nbnoer, numpro, jenvoi1, jrecep1
    integer(kind=4) :: tag4, numpr4, n4e, n4r
    character(len=8) :: nomtyp(MT_NTYMAX), k8bid
    character(len=19) :: tag_name, comm_name, joints
    character(len=24) :: domj, recv, send
    character(len=32) :: nojoie, nojoir
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: point(:) => null()
    integer(kind=8), pointer :: typma(:) => null()
    integer(kind=8), pointer :: v_comm(:) => null()
    integer(kind=8), pointer :: v_tag(:) => null()
    integer(kind=8), pointer :: v_dom(:) => null()
    real(kind=8), pointer :: coordo(:) => null()
    aster_logical, pointer :: par_seq(:) => null()
    mpi_int :: mrank, msize, world
    real(kind=8) :: start, end
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
    if (niv > 1) then
        call cpu_time(start)
        write (ifm, *) "<IREMED_FILTRE> DEBUT CREATION FILTRE"
    end if
!
    call jeexin(nomsd//'.NBNO', iret)
    if (iret .eq. 0) then
!
! 1. ==> Initialisation
!
        call asmpi_comm('GET', world)
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
!
        call wkvect(nomsd//'.NOMA', base//' V K8', 1, jma)
        zk8(jma) = nomast
        call dismoi('NB_NO_MAILLA', nomast, 'MAILLAGE', repi=nbnoeu)
        call jelira(nomast//'.TYPMAIL', 'LONMAX', nbmail)
        call jeveuo(nomast//'.COORDO    .VALE', 'L', vr=coordo)
        call jeveuo(nomast//'.CONNEX', 'L', vi=connex)
        call jeveuo(jexatr(nomast//'.CONNEX', 'LONCUM'), 'L', vi=point)
        call jeveuo(nomast//'.TYPMAIL        ', 'L', vi=typma)
!
! 2. ==> . RECUPERATION DES NB/NOMS/NBNO/NBITEM DES TYPES DE MAILLES
!          DANS CATALOGUE
!        . RECUPERATION DES TYPES GEOMETRIE CORRESPONDANT POUR MED
!        . VERIF COHERENCE AVEC LE CATALOGUE
!
        call lrmtyp(nbtyp, nomtyp, nnotyp, typgeo, renumd, &
                    modnum, nuanom, numnoa)
!
! 3. ==> Construction de la liste des noeuds interieurs
!        et de la liste contenant par proc le nombre de noeuds
!        int et le numero de depart de la numerotation
!        des noeuds dans la num globale des noeuds
!
        call jeveuo(nomast//'.NOEX', 'L', jnoex)
        call wkvect(nomsd//'.PAR_SEQ', base//' V L', 1, vl=par_seq)
        par_seq(1) = par_seqfile
        call wkvect(nomsd//'.NBNO1', base//' V I', 2*nbproc, jnbno1)
        call wkvect(nomsd//'.NBNO', base//' V I', 3, jnbno)
        call wkvect(nomsd//'.NBNOP', base//' V I', 3, iaux)
        call wkvect(nomsd//'.NOEU', base//' V I', nbnoeu, jno)
!
        nbnot = 0
        do ino = 1, nbnoeu
            if (zi(jnoex+ino-1) .eq. rang) then
                nbnot = nbnot+1
                zi(jno+ino-1) = nbnot
            else
                zi(jno+ino-1) = -1
            end if
        end do
        zi(jnbno1+2*rang) = nbnot
        call asmpi_comm_vect('MPI_SUM', 'I', nbval=2*nbproc, vi=zi(jnbno1))
        do iproc = 1, nbproc-1
            zi(jnbno1+2*iproc+1) = zi(jnbno1+2*(iproc-1))+zi(jnbno1+2*(iproc-1)+1)
        end do
        do ino = 1, nbnoeu
            if (zi(jnoex+ino-1) .eq. rang) then
                zi(jno+ino-1) = zi(jno+ino-1)+zi(jnbno1+2*rang+1)
            end if
        end do
        zi(jnbno) = zi(jnbno1+rang*2+1)+1
        zi(jnbno+1) = zi(jnbno1+rang*2)
        zi(jnbno+2) = zi(jnbno1+nbproc*2-1)+zi(jnbno1+nbproc*2-2)
!
! 4. ==> Construction du graph de communication
!
        comm_name = '&&IREMED.COMM'
        tag_name = '&&IREMED.TAG'
        call create_graph_comm(nomast, "MAILLAGE_P", nb_comm, comm_name, tag_name)
        call jeveuo(comm_name, 'L', vi=v_comm)
        call jeveuo(tag_name, 'L', vi=v_tag)
!
! 5. ==> Construction d'une numerotation globale des noeuds :
!          proc 0 en premier, proc 1 ensuite, etc.
!
        joints = nomast//".JOIN"
        domj = joints//".DOMJ"
        send = joints//".SEND"
        recv = joints//".RECV"
        if (nb_comm > 0) then
            call jeveuo(domj, 'L', vi=v_dom)
        end if
        do i_comm = 1, nb_comm
            domj_i = v_comm(i_comm)
            numpro = v_dom(domj_i)
            nojoie = jexnum(send, domj_i)
            nojoir = jexnum(recv, domj_i)
            call jeveuo(nojoie, 'L', jjoine)
            call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
            call jeveuo(nojoir, 'L', jjoinr)
            call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
            nbnoee = nbnoee/2
            nbnoer = nbnoer/2
            call wkvect('&&IRMHDF.NOEUD_NEC_E1', 'V V I', nbnoee, jenvoi1)
            call wkvect('&&IRMHDF.NOEUD_NEC_R1', 'V V I', nbnoer, jrecep1)
!
            do ino = 0, nbnoee-1
                zi(jenvoi1+ino) = zi(jno+zi(jjoine+2*ino)-1)
            end do
            n4r = to_mpi_int(nbnoer)
            n4e = to_mpi_int(nbnoee)
            tag4 = to_mpi_int(v_tag(i_comm))
            numpr4 = to_mpi_int(numpro)
            call asmpi_sendrecv_i(zi(jenvoi1), n4e, numpr4, tag4, &
                                  zi(jrecep1), n4r, numpr4, tag4, world)

            do ino = 0, nbnoer-1
                if (zi(jrecep1+ino) .ne. -1) zi(jno+zi(jjoinr+2*ino)-1) = -zi(jrecep1+ino)
            end do
            call jedetr('&&IRMHDF.NOEUD_NEC_E1')
            call jedetr('&&IRMHDF.NOEUD_NEC_R1')
        end do
!
        call jeveuo(nomast//'.MAEX', 'L', jmaex)
        call wkvect(nomsd//'.NBMA', 'V V I', 2*nbproc, jnbma)
        call wkvect(nomsd//'.MAIL', base//' V I', nbmail, jma)
        call wkvect(nomsd//'.MATY', base//' V I', MT_NTYMAX*3, jtyp)
        call wkvect(nomsd//'.MATYP', base//' V I', MT_NTYMAX*3, iaux)
        call wkvect('&&FILTRE.TYPMAILL', 'V V I', MT_NTYMAX, jtyp2)
        nbmat = 0
        do ima = 1, nbmail
            if (zi(jmaex+ima-1) == rang) then
                ityp = typma(ima)
                if (ityp .eq. MT_HEXA9) ityp = MT_HEXA8
                if (ityp .eq. MT_TETRA15) ityp = MT_TETRA10
                if (ityp .eq. MT_TETRA20) ityp = MT_TETRA4
                if (ityp .eq. MT_PYRAM19) ityp = MT_PYRAM13
                if (ityp .eq. MT_PENTA21) ityp = MT_PENTA18
                if (ityp .eq. MT_PENTA7) ityp = MT_PENTA6
                if (ityp .eq. MT_TRIA10) ityp = MT_TRIA3
                if (ityp .eq. MT_QUAD12) ityp = MT_QUAD4
                zi(jtyp2+ityp-1) = zi(jtyp2+ityp-1)+1
                zi(jma+ima-1) = zi(jtyp2+ityp-1)
                nbmat = nbmat+1
            end if
        end do

        call wkvect('&&FILTRE.TYPMAILG', 'V V I', MT_NTYMAX, jtypg)
        zi(jtypg-1+1:jtypg-1+MT_NTYMAX) = zi(jtyp2-1+1:jtyp2-1+MT_NTYMAX)
        call asmpi_comm_vect('MPI_SUM', 'I', nbval=MT_NTYMAX, vi=zi(jtypg))
        do ityp = 1, MT_NTYMAX
            if (zi(jtypg+ityp-1) .ne. 0) then
                zi(jnbma:jnbma-1+2*nbproc) = 0
                zi(jnbma+rang) = zi(jtyp2+ityp-1)
                zi(jnbma+nbproc+rang) = zi(jtyp2+ityp-1)
                call asmpi_comm_vect('MPI_SUM', 'I', nbval=2*nbproc, vi=zi(jnbma))
                do iproc = 1, nbproc-1
                    zi(jnbma+nbproc+iproc) = zi(jnbma+nbproc+iproc-1)+zi(jnbma+iproc)
                end do
                !! nombre de mailles proprio du type donne
                zi(jtyp+3*(ityp-1)+1) = zi(jnbma+rang)
                !! nombre de mailles totales du type donne
                ASSERT(zi(jnbma+2*nbproc-1) == zi(jtypg+ityp-1))
                zi(jtyp+3*(ityp-1)+2) = zi(jnbma+2*nbproc-1)
                zi(jnbma+nbproc-1) = 0
                !! numero 1ere maille du type donne
                zi(jtyp+3*(ityp-1)) = zi(jnbma+nbproc+rang-1)+1
            end if
        end do
        do ima = 1, nbmail
            if (zi(jma+ima-1) .ne. 0) then
                ityp = typma(ima)
                if (ityp .eq. MT_HEXA9) ityp = MT_HEXA8
                if (ityp .eq. MT_TETRA15) ityp = MT_TETRA10
                if (ityp .eq. MT_TETRA20) ityp = MT_TETRA4
                if (ityp .eq. MT_PYRAM19) ityp = MT_PYRAM13
                if (ityp .eq. MT_PENTA21) ityp = MT_PENTA18
                if (ityp .eq. MT_PENTA7) ityp = MT_PENTA6
                if (ityp .eq. MT_TRIA10) ityp = MT_TRIA3
                if (ityp .eq. MT_QUAD12) ityp = MT_QUAD4
                zi(jma+ima-1) = zi(jma+ima-1)+zi(jtyp+3*(ityp-1))-1
            end if
        end do
        call jedetr('&&FILTRE.TYPMAILL')
        call jedetr('&&FILTRE.TYPMAILG')
        call jedetr(nomsd//'.NBMA')
        call jedetr(nomsd//'.NBNO1')
        call jedetr(comm_name)
        call jedetr(tag_name)
    end if
!
    if (niv > 1) then
        call cpu_time(end)
        write (ifm, *) "<IREMED_FILTRE> FIN CREATION FILTRE EN ", end-start, "sec"
    end if
!
    call jedema()
!
end subroutine
