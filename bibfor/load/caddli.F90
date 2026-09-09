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
subroutine caddli(keywordfact, load, mesh, model, valeType)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/afddli.h"
#include "asterfort/aflrch.h"
#include "asterfort/assert.h"
#include "asterfort/char_excl_keyw.h"
#include "asterfort/char_impo_bloc.h"
#include "asterfort/char_read_keyw.h"
#include "asterfort/char_xfem.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/getnode.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=16), intent(in) :: keywordfact
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'DDL_IMPO/TEMP_IMPO/PRES_IMPO/SECH_IMPO'
!
! --------------------------------------------------------------------------------------------------
!
! In  keywordfact : factor keyword DDL_IMPO/TEMP_IMPO/PRES_IMPO
! In  mesh        : mesh
! In  load        : load
! In  model       : model
! In  valeType    : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: n_max_cmp = 320
    integer(kind=8) :: cmp_nb
    integer(kind=8) :: cmp_acti(n_max_cmp)
    real(kind=8) :: vale_real(n_max_cmp)
    complex(kind=8) :: vale_cplx(n_max_cmp)
    character(len=8) :: vale_func(n_max_cmp)
    character(len=16) :: cmp_name(n_max_cmp)
!
    integer(kind=8) :: iocc, ino, icmp, nume_node
    integer(kind=8) :: jdirec, jprnm, jnom, jcompt
    integer(kind=8) :: nbcmp, nbec, nbnoeu, nddli, ier
    character(len=8) :: name_node, nomg
    character(len=19) :: list_rela
    character(len=4) :: coef_type
    character(len=19) :: connex_inv, modelLigrel
    character(len=19) :: ch_xfem_stat, ch_xfem_node, ch_xfem_lnno, ch_xfem_ltno, ch_xfem_heav
    integer(kind=8) :: jnoxfl, jnoxfv
    aster_logical :: lxfem
    character(len=24) :: list_node
    integer(kind=8) :: jlino
    integer(kind=8) :: nb_node, geomDime
    aster_logical :: l_bloc, l_ocmp
    integer(kind=8) :: nb_typ_bloc
    character(len=16) :: val_t_bloc(3)
    character(len=8) :: bloc_cmp_name(39)
    integer(kind=8) :: bloc_cmp_index(39)
    integer(kind=8) :: bloc_cmp_nb
    real(kind=8) :: bloc_vale_real
    character(len=8) :: bloc_vale_fonc
    complex(kind=8) :: bloc_vale_cplx
    character(len=24) :: keywordexcl
    integer(kind=8) :: n_keyexcl, itypblc
    character(len=16) :: typblc(3), lec_typ_blc
    aster_logical :: istypblc(3)
    integer(kind=8) :: cmp_nb_depl, cmp_nb_rota, cmp_nb_fourier
    integer(kind=8) :: pointer
    aster_logical :: lcolle
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call getfac(keywordfact, nddli)
    if (nddli .eq. 0) goto 999
!
! - Initializations
!
    typblc(1) = 'DEPLACEMENT'
    typblc(2) = 'ROTATION'
    typblc(3) = 'TUYAU_FOURIER'
    list_rela = '&&CADDLI.RLLISTE'

! - Model informations
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call jeveuo(modelLigrel//'.PRNM', 'L', jprnm)
!
! - Create list of excluded keywords for using in char_read_keyw
!
    keywordexcl = '&&CADDLI.KEYWORDEXCL'
    call char_excl_keyw(keywordfact, keywordexcl, n_keyexcl)
!
! - Type of linear coefficient
!
    if (keywordfact .eq. 'DDL_IMPO') then
        coef_type = 'REEL'
    else if (keywordfact .eq. 'TEMP_IMPO') then
        coef_type = 'REEL'
    else if (keywordfact .eq. 'SECH_IMPO') then
        coef_type = 'REEL'
    else if (keywordfact .eq. 'PRES_IMPO') then
        coef_type = 'COMP'
    else
        ASSERT(.false.)
    end if
    lcolle = .false.
    call jeexin(mesh//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
!
! - Information about <GRANDEUR>
!
    if (keywordfact .eq. 'DDL_IMPO') then
        nomg = 'DEPL_R'
    else if (keywordfact .eq. 'TEMP_IMPO') then
        nomg = 'TEMP_R'
    else if (keywordfact .eq. 'SECH_IMPO') then
        nomg = 'TEMP_R'
    else if (keywordfact .eq. 'PRES_IMPO') then
        nomg = 'PRES_C'
    else
        ASSERT(.false.)
    end if
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', jnom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomg), 'LONMAX', nbcmp)
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 11)
!
! - Local coordinate system (dummy)
!
    call jelira(mesh//'.COORDO    .VALE', 'LONMAX', nbnoeu)
    nbnoeu = nbnoeu/3
    call wkvect('&&CADDLI.DIRECT', 'V V R', 3*nbnoeu, jdirec)
!
! - Xfem fields
!
    call char_xfem(mesh, model, lxfem, connex_inv, ch_xfem_stat, &
                   ch_xfem_node, ch_xfem_lnno, ch_xfem_ltno, ch_xfem_heav)
    if (lxfem) then
        call jeveuo(ch_xfem_node//'.CNSL', 'L', jnoxfl)
        call jeveuo(ch_xfem_node//'.CNSV', 'L', jnoxfv)
    end if
!
! - Loop on factor keyword
!
    do iocc = 1, nddli
!
! ----- Read mesh affectation
!
        list_node = '&&CADDLI.LIST_NODE'
        call getnode(mesh, keywordfact, iocc, ' ', list_node, &
                     nb_node)
!
! ----- No nodes (empty groups)
!
        if (nb_node .eq. 0) goto 60
        call jeveuo(list_node, 'L', jlino)
!
! ----- Detection of BLOCAGE
!
        call getvtx(keywordfact, 'BLOCAGE', iocc=iocc, nbval=3, vect=val_t_bloc, nbret=nb_typ_bloc)
        l_bloc = nb_typ_bloc .gt. 0
! ----- BLOCAGE case
!
        cmp_nb = 0
        istypblc(:) = .false.
        if (l_bloc) then
! --------- Counting components
            do itypblc = 1, nb_typ_bloc
                lec_typ_blc = val_t_bloc(itypblc) (1:lxlgut(val_t_bloc(itypblc)))
                if (typblc(1) .eq. lec_typ_blc) then
                    istypblc(1) = .true.
                else if (lec_typ_blc .eq. typblc(2)) then
                    istypblc(2) = .true.
                else if (lec_typ_blc .eq. typblc(3)) then
                    istypblc(3) = .true.
                end if
            end do
!
! --------- Data preparation for BLOCAGE
!
!            ASSERT(nb_typ_bloc.eq.1)
            call char_impo_bloc(nomg, istypblc, bloc_cmp_nb, bloc_cmp_name, bloc_cmp_index, &
                                bloc_vale_real, bloc_vale_cplx, bloc_vale_fonc)
            call wkvect('&&CADDLI.ICOMPT', 'V V I', bloc_cmp_nb, jcompt)

!
! --------- Loop on nodes
!
            do ino = 1, nb_node
                nume_node = zi(jlino-1+ino)
                name_node = int_to_char8(nume_node, lcolle, mesh, "NOEUD")
                cmp_nb = 0
                cmp_nb_depl = 0
                cmp_nb_rota = 0
                cmp_nb_fourier = 0
                pointer = 0
! ---------------- Components of DEPLACEMENT of the node to be blocked
                if (istypblc(1)) then
                    do icmp = 1, 3
                        if (exisdg(zi(jprnm-1+(nume_node-1)*nbec+1), bloc_cmp_index(icmp))) then
                            cmp_nb = cmp_nb+1
                            cmp_nb_depl = cmp_nb_depl+1
                            cmp_acti(cmp_nb) = 1
                            cmp_name(cmp_nb) = bloc_cmp_name(icmp)
                            vale_real(cmp_nb) = bloc_vale_real
                            vale_cplx(cmp_nb) = bloc_vale_cplx
                            vale_func(cmp_nb) = bloc_vale_fonc
                        end if
                    end do
                    pointer = pointer+3
                end if
! ---------------- Components of ROTATION of the node to be blocked
                if (istypblc(2)) then
                    do icmp = pointer+1, pointer+3
                        if (exisdg(zi(jprnm-1+(nume_node-1)*nbec+1), bloc_cmp_index(icmp))) then
                            cmp_nb = cmp_nb+1
                            cmp_nb_rota = cmp_nb_rota+1
                            cmp_acti(cmp_nb) = 1
                            cmp_name(cmp_nb) = bloc_cmp_name(icmp)
                            vale_real(cmp_nb) = bloc_vale_real
                            vale_cplx(cmp_nb) = bloc_vale_cplx
                            vale_func(cmp_nb) = bloc_vale_fonc
                        end if
                    end do
                    pointer = pointer+3
                end if

! ---------------- Components of TUYAU_FOURIER of the node to be blocked

                if (istypblc(3)) then
                    do icmp = pointer+1, bloc_cmp_nb
                        if (exisdg(zi(jprnm-1+(nume_node-1)*nbec+1), bloc_cmp_index(icmp))) then
                            cmp_nb = cmp_nb+1
                            cmp_nb_fourier = cmp_nb_fourier+1
                            cmp_acti(cmp_nb) = 1
                            cmp_name(cmp_nb) = bloc_cmp_name(icmp)
                            vale_real(cmp_nb) = bloc_vale_real
                            vale_cplx(cmp_nb) = bloc_vale_cplx
                            vale_func(cmp_nb) = bloc_vale_fonc
                        end if
                    end do
!!!---------------- Check if the node has ddl of TUYAU_FOURIER
                    if (cmp_nb_fourier .eq. 0) then
                        call utmess('F', 'CHARGES2_92', sk=typblc(3))
                    end if
                end if

                call afddli(model, geomDime, nbcmp, zk8(jnom), nume_node, name_node, &
                            zi(jprnm-1+(nume_node-1)*nbec+1), 0, zr(jdirec+3*(nume_node-1)), &
                            coef_type, cmp_nb, cmp_name, cmp_acti, valeType, &
                            vale_real, vale_func, vale_cplx, zi(jcompt), list_rela, &
                            lxfem, jnoxfl, jnoxfv, ch_xfem_stat, ch_xfem_lnno, &
                            ch_xfem_ltno, connex_inv, mesh, ch_xfem_heav)
            end do
!
            call jedetr('&&CADDLI.ICOMPT')
        end if
!
! ----- Read affected components and their values
!
        call char_read_keyw(keywordfact, iocc, valeType, n_keyexcl, keywordexcl, &
                            n_max_cmp, cmp_nb, cmp_name, cmp_acti, vale_real, &
                            vale_func, vale_cplx)
        l_ocmp = cmp_nb .gt. 0

!
! ----- Other cases
!
        if (l_ocmp) then
!
! --------- Counting components
!
            call wkvect('&&CADDLI.ICOMPT', 'V V I', cmp_nb, jcompt)
!
! --------- Loop on nodes
!
            do ino = 1, nb_node
                nume_node = zi(jlino-1+ino)
                name_node = int_to_char8(nume_node, lcolle, mesh, "NOEUD")
                call afddli(model, geomDime, nbcmp, zk8(jnom), nume_node, name_node, &
                            zi(jprnm-1+(nume_node-1)*nbec+1), 0, zr(jdirec+3*(nume_node-1)), &
                            coef_type, cmp_nb, cmp_name, cmp_acti, valeType, &
                            vale_real, vale_func, vale_cplx, zi(jcompt), list_rela, &
                            lxfem, jnoxfl, jnoxfv, ch_xfem_stat, ch_xfem_lnno, &
                            ch_xfem_ltno, connex_inv, mesh, ch_xfem_heav)
            end do
            do icmp = 1, cmp_nb
                if (zi(jcompt-1+icmp) .eq. 0) then
                    call utmess('F', 'CHARGES2_45', sk=cmp_name(icmp))
                end if
            end do
        end if
!
60      continue
!
        call jedetr('&&CADDLI.ICOMPT')
        call jedetr(list_node)
!
    end do

! - Final linear relation affectation
    call aflrch(list_rela, load, 'LIN')
!
    call jedetr('&&CADDLI.DIRECT')
    call jedetr(keywordexcl)
    if (lxfem) then
        call jedetr(connex_inv)
        call detrsd('CHAM_NO_S', ch_xfem_node)
        call detrsd('CHAM_ELEM_S', ch_xfem_stat)
        call detrsd('CHAM_ELEM_S', ch_xfem_lnno)
        call detrsd('CHAM_ELEM_S', ch_xfem_ltno)
    end if
!
999 continue
    call jedema()
end subroutine
