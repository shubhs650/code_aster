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
subroutine irvari(ifi, field_med, vari_elga, field_loca, ligrel, &
                  nb_cmp_sele, cmp_name_sele, partie, numpt, instan, &
                  nume_store, nbmaec, limaec, result, cara_elem, &
                  carael, lfichUniq, codret)
!
    use MGIS_module
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cescrm.h"
#include "asterfort/cesexi.h"
#include "asterfort/codent.h"
#include "asterfort/comp_meca_pvar.h"
#include "asterfort/comp_meca_uvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/impres_component_hpc.h"
#include "asterfort/irceme.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8), intent(in) :: ifi
    character(len=64), intent(in) :: field_med
    character(len=19), intent(in) :: vari_elga
    character(len=8), intent(in) :: field_loca
    character(len=19), intent(in) :: ligrel
    integer(kind=8), intent(in) :: nb_cmp_sele
    character(len=*), intent(in) :: cmp_name_sele(*)
    character(len=*), intent(in) :: partie
    integer(kind=8), intent(in) :: numpt
    real(kind=8), intent(in) :: instan
    integer(kind=8), intent(in) :: nume_store
    integer(kind=8), intent(in) :: nbmaec
    integer(kind=8), intent(in) :: limaec(*)
    character(len=8), intent(in) :: result
    character(len=8), intent(in) :: cara_elem, carael
    aster_logical, intent(in) :: lfichUniq
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Post-treatment (IMPR_RESU)
!
! Create VARI_ELGA_NOMME for name of internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_store       : index to store in results
! In  field_med        : name of MED field
! In  field_loca       : localization of field
!                        /'ELNO'/'ELGA'/'ELEM'
! In  result           : name of results datastructure
! In  ligrel           : name of ligrel of the field
! Out codret           : error code
!                        0   - Everything is OK
!                        200 - External behaviour (MFRONT/UMAT) or multifibers
!                        300 - No external state variables (everything is elastic)
!       IFI    : UNITE LOGIQUE D'IMPRESSION DU CHAMP
!       PARTIE : IMPRESSION DE LA PARTIE IMAGINAIRE OU REELLE POUR
!                UN CHAMP COMPLEXE
!       nb_cmp_sele  : NOMBRE DE COMPOSANTES A ECRIRE
!       cmp_name_sele : NOMS DES COMPOSANTES A ECRIRE
!       NUMPT  : NUMERO DE PAS DE TEMPS
!       INSTAN : VALEUR DE L'INSTANT A ARCHIVER
!       NUMORD : NUMERO D'ORDRE DU CHAMP
!       NBMAEC : NOMBRE DE MAILLES A ECRIRE (0, SI TOUTES LES MAILLES)
!       LIMAEC : LISTE DES MAILLES A ECRIRE SI EXTRAIT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iMapZone, iCell, i_pt, i_vari, i_vari_redu, i_spt
    integer(kind=8) :: nb_vari, nb_pt, nb_spt, nb_vari_zone
    integer(kind=8) :: nb_vari_redu, nbMapZone, nbCell, nb_vari_maxi, nbCellMesh, nb_elem_zone
    integer(kind=8) :: nt_vari, codret_dummy, nbCmpDyna
    integer(kind=8) :: posit, iret, affe_type, affe_indx, nume_elem
    integer(kind=8) :: jv_elga_cesd, jv_elga_cesl, jv_elgr_cesd, jv_elgr_cesl, jv_elga, jv_elgr
    integer(kind=8) :: ncmpvl, jicmp, poscmp, iCmp
    character(len=7) :: saux07
    character(len=8) :: saux08
    character(len=8), parameter :: base_name = '&&IRVARI'
    character(len=19) :: compor
    character(len=19), parameter :: vari_elga_s = '&&IRVARI.VARIELGA_S'
    character(len=19), parameter :: vari_elgr = '&&IRVARI.VARIELGR'
    character(len=19), parameter :: vari_elgr_s = '&&IRVARI.VARIELGR_S'
    character(len=19) :: vari_link
    character(len=19), parameter :: vari_redu = '&&IRVARI.VARIREDU'
    integer(kind=8), pointer :: v_vari_link(:) => null()
    character(len=16), pointer :: v_vari_redu(:) => null()
    character(len=19), parameter :: label_med = '&&IRVARI.LABELMED'
    character(len=19), parameter :: label_vxx = '&&IRVARI.LABELVXX'
    character(len=8), pointer :: v_label_vxx(:) => null()
    character(len=16), pointer :: v_label_med(:) => null()
    character(len=64) :: nomres
    real(kind=8), pointer :: v_elgr_cesv(:) => null()
    real(kind=8), pointer :: v_elga_cesv(:) => null()
    integer(kind=8), pointer :: v_compor_desc(:) => null()
    integer(kind=8), pointer :: v_compor_lima(:) => null()
    integer(kind=8), pointer :: v_compor_lima_lc(:) => null()
    character(len=19), parameter :: comporInfo = '&&IRVARI.INFO'
    integer(kind=8), pointer :: comporInfoInfo(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    character(len=16) :: field_type
    character(len=24)  :: indcmp
    character(len=32)  :: nomgd
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    field_type = 'VARI_ELGA'
    ASSERT(field_loca .eq. 'ELGA')
    codret = 0
!
! - Get name of <CARTE> COMPOR
!
    call rsexch('F', result, 'COMPORTEMENT', nume_store, compor, iret)
    if (iret .ne. 0) then
        call utmess('F', 'RESULT1_5')
    end if
    if (hasMFront(compor)) then
        call utmess('F', "COMPOR6_7")
    end if

! - Prepare informations about internal variables
    call comp_meca_pvar(ligrel, compor, comporInfo)

! - Access to informations
    call jeveuo(comporInfo(1:19)//'.INFO', 'L', vi=comporInfoInfo)
    nbCellMesh = comporInfoInfo(1)
    nbMapZone = comporInfoInfo(2)
    nb_vari_maxi = comporInfoInfo(3)
    nt_vari = comporInfoInfo(4)
!
    if (nt_vari .eq. 0) then
        codret = 300
        goto 999
    end if
    call jeveuo(comporInfo(1:19)//'.ZONE', 'L', vi=comporInfoZone)
!
! - Transform VARI_ELGA in VARI_ELGA_S
!
    call celces(vari_elga, 'V', vari_elga_s)
    call jeveuo(vari_elga_s//'.CESD', 'L', jv_elga_cesd)
    call jeveuo(vari_elga_s//'.CESL', 'L', jv_elga_cesl)
    call jeveuo(vari_elga_s//'.CESV', 'L', vr=v_elga_cesv)
!
! - Create list of internal variables and link to zone in <CARTE> COMPOR
!
    call comp_meca_uvar(comporInfo, base_name, vari_redu, nb_vari_redu, codret)

    call dismoi('NOM_GD', vari_elga_s, 'CHAMP', repk=nomgd)
    indcmp = '&&IRVARI.CMPINDIR'
!
! - Create component list in parallel context
!
    if (lfichUniq) then
        call impres_component_hpc(nomgd, vari_redu, ncmpvl, nb_vari_redu, indcmp)
    else
        if (nb_vari_redu .ne. 0) then
            call wkvect(indcmp, 'V V I', nb_vari_redu, jicmp)
            do iCmp = 1, nb_vari_redu
                zi(jicmp+iCmp-1) = iCmp
            end do
        end if
    end if

    call jeveuo(vari_redu, 'L', vk16=v_vari_redu)
! - Behaviours that cannot give name of internal state variables
    if (codret .eq. 200) then
        goto 999
    end if
! - Only elastic behaviours
    if (nb_vari_redu .eq. 0) then
        codret = 300
        goto 999
    end if
    call jeveuo(indcmp, 'L', jicmp)
!
! - Access to <CARTE> COMPOR
!
    call jeveuo(compor//'.DESC', 'L', vi=v_compor_desc)
    call jeveuo(jexnum(compor//'.LIMA', 1), 'L', vi=v_compor_lima)
    call jeveuo(jexatr(compor//'.LIMA', 'LONCUM'), 'L', vi=v_compor_lima_lc)
!
! - Prepare objects to reduced list of internal variables
!
    call wkvect(label_vxx, 'V V K8', nb_vari_redu, vk8=v_label_vxx)
    call wkvect(label_med, 'V V K16', 2*nb_vari_redu, vk16=v_label_med)
    do i_vari_redu = 1, nb_vari_redu
        call codent(i_vari_redu, 'G', saux07)
        v_label_vxx(i_vari_redu) = 'V'//saux07
        v_label_med(2*(i_vari_redu-1)+1) = 'V'//saux07
        v_label_med(2*(i_vari_redu-1)+2) = v_vari_redu(i_vari_redu)
    end do
!
! - Create VARI_ELGR_S on reduced list of internal variables
!
    call cescrm('V', vari_elgr_s, field_loca, 'VARI_R', nb_vari_redu, &
                v_label_vxx, vari_elga_s)
    call jeveuo(vari_elgr_s//'.CESD', 'L', jv_elgr_cesd)
    call jeveuo(vari_elgr_s//'.CESL', 'L', jv_elgr_cesl)
    call jeveuo(vari_elgr_s//'.CESV', 'L', vr=v_elgr_cesv)
!
! - Fill VARI_ELGR_S on reduced list of internal variables
!
    do iMapZone = 1, nbMapZone
        nb_elem_zone = comporInfoZone(iMapZone)
        if (nb_elem_zone .ne. 0) then
!
! --------- Get object to link zone to internal variables
!
            call codent(iMapZone, 'G', saux08)
            vari_link = base_name//saux08
            call jeveuo(vari_link, 'L', vi=v_vari_link)
!
! --------- Access to current zone in CARTE
!
            nbCell = 0
            affe_type = v_compor_desc(1+3+(iMapZone-1)*2)
            affe_indx = v_compor_desc(1+4+(iMapZone-1)*2)
            if (affe_type .eq. 3) then
                nbCell = v_compor_lima_lc(1+affe_indx)-v_compor_lima_lc(affe_indx)
                posit = v_compor_lima_lc(affe_indx)
            else if (affe_type .eq. 1) then
                nbCell = nbCellMesh
                posit = 0
            else
                ASSERT(.false.)
            end if
            call jelira(jexnum(comporInfo(1:19)//'.VARI', iMapZone), 'LONMAX', nb_vari_zone)
!
! --------- Loop on elements in zone of CARTE
!
            do iCell = 1, nbCell
                if (affe_type .eq. 3) then
                    nume_elem = v_compor_lima(posit+iCell-1)
                else if (affe_type .eq. 1) then
                    nume_elem = iCell
                else
                    ASSERT(.false.)
                end if
                nb_pt = zi(jv_elga_cesd-1+5+4*(nume_elem-1)+1)
                nb_spt = zi(jv_elga_cesd-1+5+4*(nume_elem-1)+2)
                nb_vari = zi(jv_elga_cesd-1+5+4*(nume_elem-1)+3)
                do i_pt = 1, nb_pt
                    do i_spt = 1, nb_spt
                        do i_vari = 1, nb_vari
                            call cesexi('C', jv_elga_cesd, jv_elga_cesl, nume_elem, i_pt, &
                                        i_spt, i_vari, jv_elga)
                            if (jv_elga .gt. 0 .and. i_vari .le. nb_vari_zone) then
                                i_vari_redu = v_vari_link(i_vari)
                                if (i_vari_redu .ne. 0) then
                                    poscmp = zi(jicmp+i_vari_redu-1)
                                    if (poscmp .eq. 0) cycle
                                    call cesexi('C', jv_elgr_cesd, jv_elgr_cesl, nume_elem, i_pt, &
                                                i_spt, poscmp, jv_elgr)
                                    ASSERT(jv_elgr .ne. 0)
                                    jv_elgr = abs(jv_elgr)
                                    v_elgr_cesv(jv_elgr) = v_elga_cesv(jv_elga)
                                    zl(jv_elgr_cesl-1+jv_elgr) = .true.
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end if
    end do
!
! - Transform VARI_ELGR_S in VARI_ELGR
!
    nomres = field_med(1:8)//'VARI_ELGA_NOMME'
    call cescel(vari_elgr_s, ligrel, ' ', ' ', 'OUI', &
                nume_elem, 'V', vari_elgr, 'F', codret_dummy)
!
! - Write in MED file
!
    nbCmpDyna = 0
    call irceme(ifi, nomres, vari_elgr, field_loca, ligrel, &
                nb_cmp_sele, cmp_name_sele, label_med, partie, numpt, &
                instan, nume_store, nbmaec, limaec, cara_elem, &
                carael, field_type, nbCmpDyna, lfichUniq, codret)

!
999 continue
!
! - Cleaning
!
    call detrsd('CHAM_ELEM_S', vari_elga_s)
    call detrsd('CHAM_ELEM_S', vari_elgr_s)
    call detrsd('CHAM_ELEM_S', vari_elgr)
    call jedetr(comporInfo(1:19)//'.ZONE')
    call jedetr(comporInfo(1:19)//'.INFO')
    call jedetr(comporInfo(1:19)//'.ELEM')
    call jedetr(comporInfo(1:19)//'.RELA')
    call jedetc('V', comporInfo(1:19)//'.VARI', 1)
    call jedetr(vari_redu)
    call jedetr(label_vxx)
    call jedetr(label_med)
    call jedetr(indcmp)
    do iMapZone = 1, nbMapZone
        call codent(iMapZone, 'G', saux08)
        vari_link = base_name//saux08
        call jedetr(vari_link)
    end do
!
    call jedema()
!
end subroutine
