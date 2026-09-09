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

subroutine ssgngm(noma, iocc, nbgnaj)
    implicit none
#include "jeveux.h"
#include "asterfort/addGroupNode.h"
#include "asterfort/utlisi.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/gmgnre.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/existGrpMa.h"
#include "asterfort/addGrpNo.h"
#include "asterfort/cleanListOfGrpMa.h"
#include "asterfort/checkListOfGrpMa.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/vector_comm_ghost_values_i.h"
!
    character(len=8) :: noma
!     BUT: AJOUTER S'IL LE FAUT A UN MAILLAGE DES GROUP_NO DE MEME NOM
!          QUE LES GROUP_MA ET CONTENANT LES NOEUDS DE CES MAILLES.
!
!     IN: NOMA  : NOM DU MAILLAGE.
!         IOCC  : NUMERO D'OCCURENCE DU MOT CLEF CREA_GROUP_NO
!-----------------------------------------------------------------------
!
    character(len=8) :: k8b, koui
    character(len=16) :: selec
    character(len=24) :: grpma, grpno, nomgno, nomgma, grpmap
    integer(kind=8) :: nbgnaj, nbis, ii, ino, n1_save
    integer(kind=8) :: ialii1, ialii2, n, ntrou
    character(len=24), pointer :: v_gno(:) => null(), v_gma(:) => null()
    integer(kind=8), pointer :: v_trav(:) => null(), v_trav_v(:) => null()
!
! DEB-------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ialgma, ialima, ialino, ianbno, iangno
    integer(kind=8) :: ibid, ier, iocc, iret
    integer(kind=8) :: n1, nb, nbgma, nbgno, nbma, nbnoto, no
    aster_logical :: l_exi_in_grp, l_exi_in_grp_p, l_parallel_mesh, l_added_grpno
!-----------------------------------------------------------------------
    call jemarq()
    grpma = noma//'.GROUPEMA       '
    grpmap = noma//'.PAR_GRPMAI'
    grpno = noma//'.GROUPENO       '
    nbgnaj = 0
    nbgno = 0
    selec = 'TOUS'
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoto)
    l_parallel_mesh = isParallelMesh(noma)
!
    if (nbnoto .eq. 0) goto 60
!
!     ---  CAS : "TOUT_GROUP_MA"
!     --------------------------
    call getvtx('CREA_GROUP_NO', 'TOUT_GROUP_MA', iocc=iocc, scal=koui, nbret=n1)
    if (n1 .eq. 1) then
        if (l_parallel_mesh) then
            call jelira(noma//'.PAR_GRPMAI', 'NOMMAX', nbgma)
        else
            call jelira(noma//'.GROUPEMA', 'NMAXOC', nbgma)
        end if
        call wkvect('&&SSGNGM.LISTE_GMA', 'V V K24', nbgma, ialgma)
        if (l_parallel_mesh) then
            do i = 1, nbgma
                call jenuno(jexnum(grpmap, i), zk24(ialgma-1+i))
            end do
        else
            do i = 1, nbgma
                call jenuno(jexnum(grpma, i), zk24(ialgma-1+i))
            end do
            call cleanListOfGrpMa(noma, zk24(ialgma), nbgma, ASTER_TRUE, iret)
        end if
        iangno = ialgma
!
!     ---  CAS : "GROUP_MA"
!     ---------------------
    else
!
!      CRITERE DE SELECTION
        call getvtx('CREA_GROUP_NO', 'CRIT_NOEUD', iocc=iocc, scal=selec, nbret=ibid)
        call getvtx('CREA_GROUP_NO', 'GROUP_MA', iocc=iocc, nbval=0, &
                    vect=k8b, nbret=nbgma)
        call wkvect('&&SSGNGM.LISTE_GMA', 'V V K24', -nbgma, ialgma)
        call getvtx('CREA_GROUP_NO', 'GROUP_MA', iocc=iocc, nbval=-nbgma, &
                    vect=zk24(ialgma), nbret=nb)
        nbgma = nb
        call getvtx('CREA_GROUP_NO', 'NOM', iocc=iocc, nbval=0, nbret=no)
        nbgno = -no
        if (no .ne. 0) then
            if ((nbgno .ne. nbgma) .and. (nbgno .ne. 1)) then
                call utmess('F', 'MODELISA7_8')
            end if
            call wkvect('&&SSGNGM.NOM_GNO', 'V V K24', nbgno, iangno)
            call getvtx('CREA_GROUP_NO', 'NOM', iocc=iocc, nbval=nbgno, vect=zk24(iangno), &
                        nbret=no)
        else
            iangno = ialgma
        end if
        call checkListOfGrpMa(noma, zk24(ialgma), nbgma, ASTER_FALSE)
        ier = 0
        do i = 1, nbgma
            nomgma = zk24(ialgma-1+i)
            call existGrpMa(noma(1:8), nomgma, l_exi_in_grp, l_exi_in_grp_p)
            if (.not. l_exi_in_grp_p) then
                ier = ier+1
                call utmess('E', 'ELEMENTS_62', sk=nomgma)
            end if
        end do
        ASSERT(ier .eq. 0)
    end if
    if (nbgma .eq. 0) goto 60
!
    call wkvect('&&SSGNGM.LISTE_NO ', 'V V I', nbnoto, ialino)
    call wkvect('&&SSGNGM.TRAV ', 'V V I', nbnoto, vi=v_trav)
    call wkvect('&&SSGNGM.NB_NO    ', 'V V I', nbgma, ianbno)
    if (l_parallel_mesh) then
        call wkvect('&&SSGNGM.TRAV_VERIF', 'V V I', nbnoto, vi=v_trav_v)
    end if
! ---------------------------------------------------------------------
    call jelira(grpno, 'NMAXOC', nbis)
! ---------------------------------------------------------------------
!     -- ON AJOUTE LES NOUVEAUX GROUPES:
!
    do i = 1, nbgma
        nomgma = zk24(ialgma-1+i)
        call existGrpMa(noma(1:8), nomgma, l_exi_in_grp, l_exi_in_grp_p)
        if (l_parallel_mesh) then
            v_trav(:) = 0
            v_trav_v(:) = 0
            n1 = 0
        end if
        if (l_exi_in_grp) then
            call jelira(jexnom(grpma, nomgma), 'LONUTI', nbma)
            call jeveuo(jexnom(grpma, nomgma), 'L', ialima)
            call gmgnre(noma, nbnoto, v_trav, zi(ialima), nbma, &
                        zi(ialino), zi(ianbno-1+i), selec)
            n1 = zi(ianbno-1+i)
            ! FUSION DES NOEUDS ISSUES DES GROUP_MA
            if ((nbgno .eq. 1) .and. (nbgno .ne. nbgma)) then
                if (i .eq. 1) then
                    if (n1 .gt. nbis) then
                        nbis = 2*n1
                        !call jedetr('&&SSCGNO.LII1')
                        call wkvect('&&SSCGNO.LITMP1', 'V V I', nbis, ialii1)
                    end if
                    n = n1
                    do ii = 1, n
                        zi(ialii1-1+ii) = zi(ialino-1+ii)
                    end do
                end if
                if (i .ne. 1) then
                    call wkvect('&&SSCGNO.LITMP2', 'V V I', nbis, ialii2)
                    call utlisi('UNION', zi(ialii1), n, zi(ialino), n1, &
                                zi(ialii2), nbis, ntrou)
                    if (ntrou .lt. 0) then
                        nbis = -2*ntrou
                        call jedetr('&&SSCGNO.LITMP2')
                        call wkvect('&&SSCGNO.LITMP2', 'V V I', nbis, ialii2)
                        call utlisi('UNION', zi(ialii1), n, zi(ialino), n1, &
                                    zi(ialii2), nbis, ntrou)
                        call jedetr('&&SSCGNO.LITMP1')
                        call wkvect('&&SSCGNO.LITMP1', 'V V I', nbis, ialii1)
                    end if
                    n = ntrou
                    do ii = 1, n
                        zi(ialii1-1+ii) = zi(ialii2-1+ii)
                    end do
                end if
            end if
        end if

        if (l_exi_in_grp .or. (.not. l_exi_in_grp .and. l_exi_in_grp_p)) then
            if ((nbgno .eq. nbgma) .or. (nbgno .eq. 0)) then
                nomgno = zk24(iangno-1+i)
                if (.not. l_parallel_mesh) then
                    call addGrpNo(noma(1:8), nomgno, zi(ialino), n1, l_added_grpno)
                    if (l_added_grpno) then
                        nbgnaj = nbgnaj+1
                    end if
                else
                    do ino = 1, nbnoto
                        v_trav_v(ino) = v_trav(ino)
                    end do
                    call vector_comm_ghost_values_i(v_trav, noma, 'BIDIR')
                    n1_save = 0
                    do ino = 1, nbnoto
                        if (v_trav_v(ino) .eq. 0 .and. v_trav(ino) .ne. 0) then
                            zi(ialino+n1+n1_save) = ino
                            n1_save = n1_save+1
                        end if
                    end do
                    n1 = n1_save+n1
                    if (n1 .ne. 0) then
                        call addGroupNode(noma, 1)
                        call addGrpNo(noma(1:8), nomgno, zi(ialino), n1, l_added_grpno)
                        if (l_added_grpno) then
                            nbgnaj = nbgnaj+1
                        end if
                    end if
                end if
            end if
        end if
    end do

    if ((nbgno .eq. 1) .and. (nbgno .ne. nbgma)) then
        nomgno = zk24(iangno)
        ! Je ne comprends pas la raison du if precedent alors dans le cas l_parallel_mesh
        ! on plante
        ASSERT(.not. l_parallel_mesh)
        call addGrpNo(noma(1:8), nomgno, zi(ialii1), n, l_added_grpno)
        if (l_added_grpno) then
            nbgnaj = nbgnaj+1
        end if
    end if

!
60  continue
    call jedetr('&&SSGNGM.LISTE_GMA')
    call jedetr('&&SSGNGM.NOM_GNO')
    call jedetr('&&SSGNGM.LISTE_NO')
    call jedetr('&&SSGNGM.TRAV')
    call jedetr('&&SSGNGM.TRAV_VERIF')
    call jedetr('&&SSGNGM.NB_NO')
    call jedetr('&&SSCGNO.LITMP1')
    call jedetr('&&SSCGNO.LITMP2')
    call jedema()

end subroutine
