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

subroutine crnustd(numddl)
    implicit none
#include "asterc/asmpi_allgather_i.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_recv_i.h"
#include "asterc/asmpi_send_i.h"
#include "asterc/asmpi_sendrecv_i.h"
#include "asterf_config.h"
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/create_graph_comm.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    character(len=14) :: numddl
!
! --------------------------------------------------------------------------------------------------
!
!  Le but est de créer une numérotation parallèle qui soit toujours la même quelque soit le
!  nombre de sous-domaines
!
!  .NUME.DEEG(2*nbbdl) : 2*(i_ddl-1) + 1 -> numéro global du noeud physique ou tardif
!                      : 2*(i_ddl-1) + 2 -> numéro de la composante (comme DEEQ)
!
!  .NUME.NULS(nbddl) : iddl -> numéro d'équation globale (différent de .NUME.NULG)
!    qui ne dépend pas de la partition du maillage
!
! --------------------------------------------------------------------------------------------------
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
!
    integer(kind=8) :: ili, idprn1, idprn2, ntot, lonmax, nbno_prno, iLigrT
    integer(kind=8) :: nbddll, ino, iret, nbcmp, iec, iret1, iret2, jjoine
    integer(kind=8) :: numero_noeud, numero_cmp, rang, nbproc, jrefn, nbLigrTot
    integer(kind=8) :: nec, numloc, dime, nbddl_lag, i_ddl, nddl, nddlg, nddll
    integer(kind=8) :: nbno, nbno_lc, nbno_gl, nbno_max, nbddll_gl, numnoe
    integer(kind=8) :: nbddl_phys_gl, nbddl_lag_gl, i_join, jnujoi1, jnujoi2
    integer(kind=8) :: numpro, nbnoee, nbnoer, jjoinr, poscom, numno1, numno2
    integer(kind=8) :: lgenve1, lgenvr1, jencod, jenco2, lgenve2, lgenvr2
    integer(kind=8) :: jaux, nb_ddl_envoi, jrecep1, jenvoi1, jenvoi2, jrecep2
    integer(kind=8) :: nbddl, ncmpmx, iad, jcpnec, jcpne2, ico2, icmp, curpos
    integer(kind=8) :: nbno_lili_lc, nbno_lili_gl, nb_comm, domj_i, numpr2
    mpi_int :: mrank, msize, mpicou, nbno4
    mpi_int :: tag4, numpr4, n4e, n4r
    integer(kind=8), pointer :: v_noex(:) => null()
    integer(kind=8), pointer :: v_deeq(:) => null()
    integer(kind=8), pointer :: v_nequ(:) => null()
    integer(kind=8), pointer :: v_delg(:) => null()
    integer(kind=8), pointer :: v_owner(:) => null()
    integer(kind=8), pointer :: v_nulg(:) => null()
    integer(kind=8), pointer :: v_nuls(:) => null()
    integer(kind=8), pointer :: v_ddlc(:) => null()
    integer(kind=8), pointer :: v_nddl(:) => null()
    integer(kind=8), pointer :: v_gddl(:) => null()
    integer(kind=8), pointer :: v_tddl(:) => null()
    integer(kind=8), pointer :: v_deeg(:) => null()
    integer(kind=8), pointer :: v_linulg(:) => null()
    integer(kind=8), pointer :: v_comm(:) => null()
    integer(kind=8), pointer :: v_tag(:) => null()
    integer(kind=8), pointer :: v_dom(:) => null()
    integer(kind=8), pointer :: v_gco(:) => null()
    integer(kind=4), pointer :: v_pgid(:) => null()
    integer(kind=8), pointer :: v_lilt(:) => null()
    aster_logical :: lligrel_cp
!
    character(len=8) :: k8bid, mesh, nomgdr
    character(len=19) :: nomlig, tag_name, comm_name, nume_equa, meshj, joints
    character(len=24) :: owner, linulg, gcom, pgid, lilt
    character(len=32) :: nojoie, nojoir
    character(len=24), pointer :: tco(:) => null()
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
#define zzprno(ili,nunoel,l)  zi(idprn1-1+zi(idprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
!
    call jemarq()
    call asmpi_comm('GET', mpicou)
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    DEBUG_MPI('crnustd', rang, nbproc)
!
    nume_equa = numddl//".NUME"
    call jeveuo(nume_equa//'.REFN', 'L', jrefn)
    mesh = zk24(jrefn) (1:8)
    nomgdr = zk24(jrefn+1) (1:8)
    meshj = mesh//".JOIN"
    call jeveuo(meshj//'.GCOM', 'L', vi=v_gco)
    call jeveuo(meshj//'.PGID', 'L', vi4=v_pgid)
    lilt = nume_equa//'.LILT'
    call jeveuo(lilt, 'L', vi=v_lilt)
    call jelira(lilt, 'LONMAX', ival=nbLigrTot)
!
    call jeveuo(mesh//'.DIME', 'L', dime)
    call jeveuo(mesh//'.NUNOLG', 'L', vi=v_nulg)
    call jeveuo(mesh//'.NOEX', 'L', vi=v_noex)
!
!   !!! VERIFIER QU'IL N'Y A PAS DE MACRO-ELTS
!   CALCUL DU NOMBRE D'ENTIERS CODES A PARTIR DE LONMAX
    call jelira(jexnum(nume_equa//'.PRNO', 1), 'LONMAX', ntot, k8bid)
    nec = ntot/zi(dime)-2
    call wkvect('&&CRNSTD.NEC', 'V V I', nec, jencod)
    call wkvect('&&CRNSTD.NEC2', 'V V I', nec, jenco2)
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgdr), 'L', iad)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgdr), 'LONMAX', ncmpmx, k8bid)
    call wkvect('&&CRNSTD.CMP', 'V V I', ncmpmx, jcpnec)
    call wkvect('&&CRNSTD.CMP2', 'V V I', ncmpmx, jcpne2)
!
    call jeveuo(nume_equa//'.DEEQ', 'L', vi=v_deeq)
    call jeveuo(nume_equa//'.NEQU', 'L', vi=v_nequ)
    call jeveuo(nume_equa//'.DELG', 'L', vi=v_delg)
!
    nbddll = v_nequ(1)
!
    call wkvect(nume_equa//'.DEEG', 'G V I', 2*nbddll, vi=v_deeg)
    call wkvect(nume_equa//".NULS", 'G V I', nbddll, vi=v_nuls)
    v_nuls(:) = -1
!
!   RECHERCHE DES ADRESSES DU .PRNO DE .NUME
    call jeveuo(nume_equa//'.PRNO', 'E', idprn1)
    call jeveuo(jexatr(nume_equa//'.PRNO', 'LONCUM'), 'L', idprn2)
!
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nbno)
! -- On compte les noeuds locaux proprio
    nbno_lc = 0
    do ino = 1, nbno
        if (v_noex(ino) == rang) then
            nbno_lc = nbno_lc+1
        end if
    end do
!
! -- Nbr de noeud total
    nbno_gl = nbno_lc
    call asmpi_comm_vect('MPI_SUM', 'I', sci=nbno_gl)
!
! -- Nbr de noeud max par proc
    nbno_max = nbno_lc
    call asmpi_comm_vect('MPI_MAX', 'I', sci=nbno_max)
!
! -- On cherche le nombre de composante par noeud physique local proprio
    call wkvect('&&CRNSTD.NODDLL', 'V V I', 2*nbno_max, vi=v_nddl)
    v_nddl(1:2*nbno_max) = -1
    nbno_lc = 0
    numloc = 0
    do ino = 1, nbno
        if (v_noex(ino) == rang) then
            nbno_lc = nbno_lc+1
! -- Pour les noeuds proprio, on garde le num global et le nombre de ddl du noeud
            v_nddl(2*(nbno_lc-1)+1) = v_nulg(ino)
            v_nddl(2*(nbno_lc-1)+2) = zzprno(1, ino, 2)
            numloc = numloc+zzprno(1, ino, 2)
        end if
    end do
!
    call wkvect('&&CRNSTD.NLGDDL', 'V V I', 2*nbproc*nbno_max, vi=v_gddl)
! -- On envoie tout les ddls - on pourrait optimiser la mémoire
    nbno4 = to_mpi_int(2*nbno_max)
    if (nbproc .gt. 1) then
        call asmpi_allgather_i(v_nddl, nbno4, v_gddl, nbno4, mpicou)
    else
        v_gddl(:) = v_nddl(:)
    end if
!
! -- On trie les ddls par noeuds
    call wkvect('&&CRNSTD.NUMDDL', 'V V I', nbno_gl, vi=v_tddl)
    v_tddl(1:nbno_gl) = -1
    do ino = 1, nbproc*nbno_max
        numero_noeud = v_gddl(2*(ino-1)+1)
        if (numero_noeud .ne. -1) then
! numero_noeud commence à 0
            v_tddl(numero_noeud+1) = v_gddl(2*(ino-1)+2)
        end if
    end do
! -- Verif
    do ino = 1, nbno_gl
        ASSERT(v_tddl(ino) >= 0)
    end do
!
! -- On calcule le décalage
    do ino = 2, nbno_gl
        v_tddl(ino) = v_tddl(ino)+v_tddl(ino-1)
    end do
    nbddl_phys_gl = v_tddl(nbno_gl)
    do ino = nbno_gl, 2, -1
        v_tddl(ino) = v_tddl(ino-1)
    end do
    v_tddl(1) = 0
!
! -- On renumerote les noeuds physiques proprio
    call wkvect('&&CRNSTD.DDLLOC', 'V V I', nbno, vi=v_ddlc)
    do i_ddl = 1, nbddll
        numero_noeud = v_deeq((i_ddl-1)*2+1)
        numero_cmp = v_deeq((i_ddl-1)*2+2)
        if (numero_noeud .gt. 0 .and. numero_cmp .gt. 0) then
            v_deeg((i_ddl-1)*2+1) = v_nulg(numero_noeud)+1
            v_deeg((i_ddl-1)*2+2) = numero_cmp
            if (v_noex(numero_noeud) == rang) then
                v_nuls(i_ddl) = v_tddl(v_nulg(numero_noeud)+1)+v_ddlc(numero_noeud)
                v_ddlc(numero_noeud) = v_ddlc(numero_noeud)+1
            end if
        end if
    end do
!
! -- On crée le graphe de comm
    comm_name = '&&CRNUSTD.COMM'
    tag_name = '&&CRNUSTD.TAG'
    call create_graph_comm(mesh, 'MAILLAGE_P', nb_comm, comm_name, tag_name)
!
    if (nb_comm > 0) then
        call jeveuo(meshj//'.DOMJ', 'L', vi=v_dom)
        call jeveuo(comm_name, 'L', vi=v_comm)
        call jeveuo(tag_name, 'L', vi=v_tag)
        joints = zk24(jrefn-1+5)
        gcom = joints//'.GCOM'
        pgid = joints//'.PGID'
        call jeveuo(gcom, 'L', vi=v_gco)
        call jeveuo(pgid, 'L', vi4=v_pgid)
        mpicou = to_mpi_int(v_gco(1))
    end if
!
! -- On renumerote les noeuds physiques non-proprio
!    Il faut faire comme dans crnlgc.F90
    do i_join = 1, nb_comm
!
!       RECHERCHE DU JOINT
        domj_i = v_comm(i_join)
        numpro = v_dom(domj_i)
        numpr2 = v_pgid(numpro+1)
        if (numpro .ne. -1) then
            nojoie = jexnum(meshj//'.SEND', domj_i)
            nojoir = jexnum(meshj//'.RECV', domj_i)
            call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
            call jeveuo(nojoir, 'L', jjoinr)
            call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
            nbnoee = nbnoee/2
            nbnoer = nbnoer/2
!
!       DES DEUX COTES LES NOEUDS NE SONT PAS DANS LE MEME ORDRE ?
            tag4 = to_mpi_int(v_tag(i_join))
            numpr4 = to_mpi_int(numpr2)
            lgenve1 = nbnoee*(1+nec)+1
            lgenvr1 = nbnoer*(1+nec)+1
            call wkvect('&&CRNSTD.NOEUD_NEC_E1', 'V V I', lgenvr1, jenvoi1)
            call wkvect('&&CRNSTD.NOEUD_NEC_R1', 'V V I', lgenve1, jrecep1)
!
            lgenve2 = nbnoee*(1+nec)+1
            lgenvr2 = nbnoer*(1+nec)+1
!
!       On commence par envoyer, le but final est de recevoir les numeros de ddl
!       On boucle donc sur les noeuds a recevoir
            nb_ddl_envoi = 0
            do jaux = 1, nbnoer
                poscom = (jaux-1)*(1+nec)+1
!           Numero local (pour le proc courant)
                numno1 = zi(jjoinr+2*(jaux-1))
!           Numero local (pour le proc d'en face)
                numno2 = zi(jjoinr+2*jaux-1)
                zi(jenvoi1+poscom) = numno2
                do iec = 1, nec
                    zi(jenvoi1+poscom+iec) = zzprno(1, numno1, 2+iec)
                end do
!           On ne souhaite recevoir les num de ddl globaux que du processeur
!           qui les possedent (cela explique le if( v_noex(numno1).eq.numpr2 ) )
!           Malgré tout, on envoie tout au proc d'en face (numeros de noeuds et
!           entiers codes) car il s'attend a avoir une liste dimensionnee a la
!           taille du raccord (cela explique le remplissage de zi(jenvoi1+...)
                if (v_noex(numno1) .eq. numpro) then
                    nb_ddl_envoi = nb_ddl_envoi+zzprno(1, numno1, 2)
                end if
            end do
            zi(jenvoi1) = nb_ddl_envoi
            n4e = to_mpi_int(lgenvr1)
            n4r = to_mpi_int(lgenve1)
            call asmpi_sendrecv_i(zi(jenvoi1), n4e, numpr4, tag4, &
                                  zi(jrecep1), n4r, numpr4, tag4, mpicou)

            call wkvect('&&CRNSTD.NUM_DDL_GLOB_E', 'V V I', zi(jrecep1)+1, jenvoi2)
            call wkvect('&&CRNSTD.NUM_DDL_GLOB_R', 'V V I', zi(jenvoi1)+1, jrecep2)

            nbddl = 0
            if (zi(jrecep1) > 0) then
!           On dimensionne le raccord a la taille du nombre de ddl attendus
!           par le proc d'en face
                do jaux = 1, nbnoee
                    poscom = (jaux-1)*(1+nec)+1
                    numno1 = zi(jrecep1+poscom)
!               Si on ne possede pas le noeud : on envoie rien (cf.
!               if( v_noex(numno1).ne.numpr2 ) then plus bas
                    if (v_noex(numno1) .ne. rang) then
                        cycle
                    end if
!
                    nddl = zzprno(1, numno1, 1)
                    nbcmp = zzprno(1, numno1, 2)
                    if (nbcmp .eq. 0) then
                        nddlg = -1
                    else
                        nddlg = v_nuls(nddl)
                    end if
!
!               Recherche des composantes demandees
                    do iec = 1, nec
                        zi(jencod+iec-1) = zzprno(1, numno1, 2+iec)
                        zi(jenco2+iec-1) = zi(jrecep1+poscom+iec)
                    end do
                    call isdeco(zi(jencod), zi(jcpnec), ncmpmx)
                    call isdeco(zi(jenco2), zi(jcpne2), ncmpmx)
                    ico2 = 0
                    do icmp = 1, ncmpmx
!                   Comme on est sur un noeud possede par le proc courant,
!                   si le proc d'en face attend un ddl (ie composante presente),
!                   il faut necessairement que la composante soit aussi presente
!                   sur ce proc (d'ou ASSERT(zi(jcpnec+icmp-1) .ne. 0))
                        if (zi(jcpne2+icmp-1) .eq. 1) then
                            ASSERT(zi(jcpnec+icmp-1) .ne. 0)
                            ASSERT(nddlg .ne. -1)
                            zi(jenvoi2+nbddl+1) = nddlg+ico2
                            nbddl = nbddl+1
                            ico2 = ico2+1
                        end if
                    end do
                end do
                ASSERT(zi(jrecep1) .eq. nbddl)
                zi(jenvoi2) = nbddl
            end if
!
            n4e = to_mpi_int(nbddl+1)
            n4r = to_mpi_int(nb_ddl_envoi+1)
            call asmpi_sendrecv_i(zi(jenvoi2), n4e, numpr4, tag4, &
                                  zi(jrecep2), n4r, numpr4, tag4, mpicou)

            if (zi(jrecep2) > 0) then
!
                curpos = 1
                do jaux = 1, nbnoer
                    numno1 = zi(jjoinr+2*(jaux-1))
!               On ne regarde que les noeuds possedes par le proc d'en face
!               puisque c'est ce qu'il a envoye
                    if (v_noex(numno1) .ne. numpro) then
                        cycle
                    end if
                    nddll = zzprno(1, numno1, 1)
                    nbcmp = zzprno(1, numno1, 2)
                    do icmp = 0, nbcmp-1
                        ASSERT(zi(jrecep2+curpos) .ne. -1)
                        v_nuls(nddll+icmp) = zi(jrecep2+curpos)
                        curpos = curpos+1
                    end do
                end do
                ASSERT(curpos-1 .eq. nb_ddl_envoi)
            end if
!
            call jedetr('&&CRNSTD.NUM_DDL_GLOB_E')
            call jedetr('&&CRNSTD.NUM_DDL_GLOB_R')
!
            call jedetr('&&CRNSTD.NOEUD_NEC_E1')
            call jedetr('&&CRNSTD.NOEUD_NEC_R1')
        end if
    end do
    call jedetr(comm_name)
    call jedetr(tag_name)
!
! -- On compte les lagranges
    nbddl_lag = 0
    nbddl_lag_gl = 0
    do iLigrT = 2, nbLigrTot
        nbno_lili_lc = 0
        if (v_lilt(iLigrT) .ne. -1) then
            ili = v_lilt(iLigrT)
            call jeexin(jexnum(nume_equa//'.PRNO', ili), iret)
            if (iret .ne. 0) then
                call jelira(jexnum(nume_equa//'.PRNO', ili), 'LONMAX', lonmax)
                nbno_prno = lonmax/(nec+2)
                call jenuno(jexnum(nume_equa//'.LILI', ili), nomlig)

                lligrel_cp = .false.
                call jeexin(nomlig//'._TCO', iret)
                if (iret .ne. 0) then
                    call jeveuo(nomlig//'._TCO', "L", vk24=tco)
                    lligrel_cp = (tco(1) .eq. 'LIGREL_CP')
                    ASSERT(.not. lligrel_cp)
                end if

                owner = nomlig//'.PNOE'
                linulg = nomlig//'.NULG'
                call jeveuo(owner, 'L', vi=v_owner)
                call jeveuo(linulg, 'L', vi=v_linulg)
                do ino = 1, nbno_prno
                    ! Le proc est proprio du noeud
                    i_ddl = zzprno(ili, ino, 1)
                    nbcmp = zzprno(ili, ino, 2)
                    ASSERT(nbcmp .eq. 1)
                    numero_noeud = -nbddl_lag_gl+v_linulg(ino)
                    numero_cmp = v_deeq((i_ddl-1)*2+2)
                    v_deeg((i_ddl-1)*2+1) = numero_noeud
                    v_deeg((i_ddl-1)*2+2) = numero_cmp
                    if (v_owner(ino) == rang) then
                        nbno_lili_lc = nbno_lili_lc+1
                        nbddl_lag = nbddl_lag+1
                        numloc = numloc+1
                        ASSERT(nbcmp .eq. 1)
                        v_nuls(i_ddl) = nbddl_phys_gl-1-v_deeg(2*(i_ddl-1)+1)
                    end if
                end do
            end if
        end if
!
! -- Nbr de noeud de Lagrange total au ligrel
        nbno_lili_gl = nbno_lili_lc
        call asmpi_comm_vect('MPI_SUM', 'I', sci=nbno_lili_gl)
        nbddl_lag_gl = nbddl_lag_gl+nbno_lili_gl
    end do
!
!   Nombre total de Lagrange
    call asmpi_comm_vect('MPI_SUM', 'I', sci=nbddl_lag)
    ASSERT(nbddl_lag == nbddl_lag_gl)
!
!   Nombre total de degré de liberté
    nbddll_gl = numloc
    call asmpi_comm_vect('MPI_SUM', 'I', sci=nbddll_gl)
!
! -- Verif nb ddl total
    ASSERT(nbddll_gl == nbddl_phys_gl+nbddl_lag_gl)
!
! -- On complete avec les joints
    do iLigrT = 2, nbLigrTot
        nbno_lili_lc = 0
        if (v_lilt(iLigrT) .eq. -1) then
            cycle
        else
            ili = v_lilt(iLigrT)
        end if
        call jenuno(jexnum(nume_equa//'.LILI', ili), nomlig)
        call create_graph_comm(nomlig, 'LIGREL', nb_comm, comm_name, tag_name)
        if (nb_comm > 0) then
            call jeveuo(comm_name, 'L', vi=v_comm)
            call jeveuo(tag_name, 'L', vi=v_tag)
            call dismoi('JOINTS', nomlig, 'LIGREL', repk=joints, arret='F')
            call jeveuo(joints//'.DOMJ', 'L', vi=v_dom)
            call jeveuo(joints//".PGID", 'L', vi4=v_pgid)
            call jeveuo(joints//".GCOM", 'L', vi=v_gco)
            mpicou = to_mpi_int(v_gco(1))
            do i_join = 1, nb_comm
                domj_i = v_comm(i_join)
                numpro = v_dom(domj_i)
                numpr2 = v_pgid(numpro+1)
                if (numpro .ne. -1) then
                    numpr4 = to_mpi_int(numpr2)
                    tag4 = to_mpi_int(v_tag(i_join))
                    nojoie = jexnum(joints//".SEND", domj_i)
                    nojoir = jexnum(joints//".RECV", domj_i)

                    call jeexin(nojoie, iret1)
                    if (iret1 .ne. 0) then
                        call jeveuo(nojoie, 'L', jjoine)
                        call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
                        call wkvect('&&CRNSTD.NUM_DDL_GLOB_E', 'V V I', nbnoee, jnujoi1)
                        do jaux = 1, nbnoee
                            numnoe = -zi(jjoine+jaux-1)
                            ASSERT(zzprno(ili, numnoe, 2) .eq. 1)
                            nddll = zzprno(ili, numnoe, 1)
                            zi(jnujoi1+jaux-1) = v_nuls(nddll)
                        end do
                        n4e = to_mpi_int(nbnoee)
                    end if

                    call jeexin(nojoir, iret2)
                    if (iret2 .ne. 0) then
                        call jeveuo(nojoir, 'L', jjoinr)
                        call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
                        call wkvect('&&CRNSTD.NUM_DDL_GLOB_R', 'V V I', nbnoer, jnujoi2)
                        n4r = to_mpi_int(nbnoer)
                    end if

                    if (rang .lt. numpro) then
                        if (iret1 .ne. 0) then
                            call asmpi_send_i(zi(jnujoi1), n4e, numpr4, tag4, mpicou)
                        end if
                        if (iret2 .ne. 0) then
                            call asmpi_recv_i(zi(jnujoi2), n4r, numpr4, tag4, mpicou)
                        end if
                    else if (rang .gt. numpro) then
                        if (iret2 .ne. 0) then
                            call asmpi_recv_i(zi(jnujoi2), n4r, numpr4, tag4, mpicou)
                        end if
                        if (iret1 .ne. 0) then
                            call asmpi_send_i(zi(jnujoi1), n4e, numpr4, tag4, mpicou)
                        end if
                    end if

                    if (iret2 .ne. 0) then
                        do jaux = 1, nbnoer
                            numnoe = -zi(jjoinr+jaux-1)
                            ASSERT(zzprno(ili, numnoe, 2) .eq. 1)
                            nddll = zzprno(ili, numnoe, 1)
                            v_nuls(nddll) = zi(jnujoi2+jaux-1)
                        end do
                    end if
                    call jedetr('&&CRNSTD.NUM_DDL_GLOB_E')
                    call jedetr('&&CRNSTD.NUM_DDL_GLOB_R')
                end if
            end do
        end if
        call jedetr(comm_name)
        call jedetr(tag_name)
        ! endif
    end do
!
! -- Verif finale
    do i_ddl = 1, nbddll
        ASSERT(v_nuls(i_ddl) .ne. -1)
    end do
!
! -- On détruit
    call jedetr('&&CRNSTD.NLGDDL')
    call jedetr('&&CRNSTD.NODDLL')
    call jedetr('&&CRNSTD.NUMDLL')
    call jedetr('&&CRNSTD.DDLLOC')
    call jedetr('&&CRNSTD.NEC')
    call jedetr('&&CRNSTD.NEC2')
    call jedetr('&&CRNSTD.CMP')
    call jedetr('&&CRNSTD.CMP2')
    call jedetc('V', '&&CRNSTD', 1)
!
    call jedema()
#else
    character(len=14) :: k14
    k14 = numddl
#endif
!
end subroutine
