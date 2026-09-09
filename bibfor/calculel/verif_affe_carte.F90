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
! aslint: disable=W0413
!
subroutine verif_affe_carte(modelLigrel, carte, comment, non_lin)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/exisdg.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/list_grma.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=19), intent(in) :: modelLigrel
    character(len=19), intent(in) :: carte
    character(len=*), intent(in) :: comment
    aster_logical, intent(in), optional ::  non_lin
!
!-----------------------------------------------------------------------
!   But :
!     Emettre des alarmes concernant les affectations douteuses d'une carte
!
!   Entrees:
!     modelLigrel     :  ligrel du modele
!     carte      :  sd_carte
!     comment    :  commentaire de description pour la sd_carte
!
!-----------------------------------------------------------------------
    character(len=3) :: tsca
    character(len=8) :: nomgd, cmpName, mesh, nommai
    character(len=24) :: lgrma(4)
    character(len=16) :: nomte
    character(len=80) :: valk(5)
    integer(kind=8) :: nbGrel, iGrel, iCmp, nbCmp, nbOption, nbte, k1, iexi, ient, iad1
    integer(kind=8) :: jnocmp, numgd, joptte, jligrmo, lielSize, iOption, optionNume, joptmod, jvale
    integer(kind=8) :: jmodeloc, nbin, kin, moloc, nbCell, cellNume, iret, typeNume, nbmapb, nbgrma
    integer(kind=8) :: nucalc, iLiel, kma, nec, nbma_verif, nbgdmx, code, decal, ico, kcmp2
    integer(kind=8), pointer :: a_un_sens(:) => null()
    integer(kind=8), pointer :: num_grel(:) => null()
    integer(kind=8), pointer :: numa_verif(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    integer(kind=8), pointer :: ptma(:) => null()
    integer(kind=8), pointer :: dg(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8)          :: list_ma_pb(5)
    aster_logical    :: verif_coef_drz
    aster_logical    :: exiq4_drz_nook, exiq4_coef_drz
    aster_logical    :: exiq3_coef_drz

!-----------------------------------------------------------------------
!
    call jemarq()

    verif_coef_drz = ASTER_FALSE
    exiq4_drz_nook = ASTER_FALSE
    exiq4_coef_drz = ASTER_FALSE
    exiq3_coef_drz = ASTER_FALSE

    call dismoi('NOM_GD', carte, 'CARTE', repk=nomgd)
    call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=nbCmp)
    call dismoi('NUM_GD', nomgd, 'GRANDEUR', repi=numgd)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jnocmp)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
    call dismoi('NB_GREL', modelLigrel, 'LIGREL', repi=nbGrel)
    call dismoi('NOM_MAILLA', modelLigrel, 'LIGREL', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCell)

!   -- 1. recuperation des objets des catalogues EF:
!   -------------------------------------------------
    call jelira('&CATA.OP.NOMOPT', 'NOMMAX', nbOption)
    call jelira('&CATA.TE.NOMTE', 'NOMMAX', nbte)
    call jeveuo('&CATA.TE.OPTTE', 'L', joptte)

!   -- 2. On calcule 2 tableaux :
!         A_UN_SENS(iGrel,iCmp)  : si 0=non , si 1=oui
!         NUM_GREL(2*(cellNume-1)+1)  : iGrel associe a la maille cellNume
!         NUM_GREL(2*(cellNume-1)+2)  : typeNume    associe a la maille cellNume
!   --------------------------------------------------------------------
    AS_ALLOCATE(vi=a_un_sens, size=nbGrel*nbCmp)
    a_un_sens = 0

    AS_ALLOCATE(vi=num_grel, size=2*nbCell)
    do iGrel = 1, nbGrel
        call jeveuo(jexnum(modelLigrel//'.LIEL', iGrel), 'L', jligrmo)
        call jelira(jexnum(modelLigrel//'.LIEL', iGrel), 'LONMAX', lielSize)
        typeNume = zi(jligrmo-1+lielSize)
        do iLiel = 1, lielSize-1
            cellNume = zi(jligrmo-1+iLiel)
            if (cellNume .gt. 0) then
                num_grel(2*(cellNume-1)+1) = iGrel
                num_grel(2*(cellNume-1)+2) = typeNume
            end if
        end do
        do iOption = 1, nbOption
            optionNume = zi(joptte-1+(typeNume-1)*nbOption+iOption)
            if (optionNume .gt. 0) then
                call jeveuo(jexnum('&CATA.TE.OPTMOD', optionNume), 'L', joptmod)
                nucalc = zi(joptmod-1+1)
                if (nucalc .ge. 0) then
                    nbin = zi(joptmod-1+2)
                    do kin = 1, nbin
                        moloc = zi(joptmod-1+3+kin)
                        call jeveuo(jexnum('&CATA.TE.MODELOC', moloc), 'L', jmodeloc)
                        if (zi(jmodeloc-1+2) .eq. numgd) then
                            ASSERT(zi(jmodeloc-1+4) .gt. 0)
                            do iCmp = 1, nbCmp
                                if (exisdg(zi(jmodeloc-1+5), iCmp)) then
                                    a_un_sens((iGrel-1)*nbCmp+iCmp) = 1
                                end if
                            end do
                        end if
                    end do
                end if
            end if
        end do
    end do
!
!   3. On parcourt les CMPS affectees volontairement dans la carte (sauf TOUT='OUI')
!
!   3.1 : on repère les mailles a vérifier (numa_verif(*)) :
    call jeveuo(carte//'.DESC', 'L', vi=desc)
    call jeveuo(carte//'.VALE', 'L', jvale)
    nbgdmx = desc(2)
!   si la carte est constante (TOUT='OUI'), on ne vérifie pas
    if (nbgdmx .eq. 1 .and. desc(3+1) .eq. 1) goto 999
!
    call etenca(carte, modelLigrel, iret)
    if (iret .gt. 0) goto 999
    call jeexin(carte//'.PTMA', iexi)
    if (iexi .eq. 0) goto 999
!
    nbma_verif = 0
    AS_ALLOCATE(vi=numa_verif, size=nbCell)
    call jeveuo(carte//'.PTMA', 'L', vi=ptma)
    do cellNume = 1, nbCell
        ient = ptma(cellNume)
        ! si la maille n'est pas affectee :
        if (ient .eq. 0) cycle
        ! si la maille est affectee par TOUT='OUI' :
        code = desc(3+2*(ient-1)+1)
        if (code .eq. 1) cycle
        ! on ne vérifie pas les mailles qui ne sont pas affectées dans le modèle
        ! (on ne saurait pas remplir le champ nomte du message)
        iGrel = num_grel(2*(cellNume-1)+1)
        if (iGrel .eq. 0) cycle
        !
        nbma_verif = nbma_verif+1
        numa_verif(nbma_verif) = cellNume
    end do
    !
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)

!
!   3.2 : on verifie les mailles a verifier (cmp par cmp) :
    do iCmp = 1, nbCmp
        cmpName = zk8(jnocmp-1+iCmp)
        verif_coef_drz = ASTER_FALSE
        ! Exceptions :
        ! E1) PESA_R / ROTA_R sont en général utilisés sans préciser les mailles
        if (nomgd .eq. 'PESA_R') cycle
        if (nomgd .eq. 'ROTA_R') cycle
        if (nomgd .eq. 'CAGNPO_R') cycle
        !
        ! E2) Valeurs fournies par le code d'AFFE_CHAR_MECA
        if (nomgd(1:5) .eq. 'FORC_' .and. cmpName .eq. 'REP') cycle
        if (nomgd(1:5) .eq. 'FORC_' .and. cmpName .eq. 'PLAN') cycle
        if (nomgd .eq. 'VENTCX_F' .and. cmpName .eq. 'FCXP') cycle
        !
        ! E3) Valeurs fournies en loucede par le code d'AFFE_CARA_ELEM
        if (nomgd .eq. 'CAMA_R' .and. cmpName .eq. 'C') cycle
        if (nomgd .eq. 'CACOQU_R' .and. cmpName .eq. 'KAPPA') cycle
        if (nomgd .eq. 'CACOQU_R' .and. cmpName .eq. 'CTOR') verif_coef_drz = ASTER_TRUE
        if (nomgd .eq. 'CAORIE_R' .and. cmpName .eq. 'ALPHA') cycle
        !
        if (nomgd .eq. 'CINFDI_R' .and. cmpName(1:3) .eq. 'REP') cycle
        if (nomgd .eq. 'CINFDI_R' .and. cmpName(1:3) .eq. 'SYM') cycle
        !
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'TSEC') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'HY1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'HZ1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'HY2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'HZ2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'EPY1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'EPY2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'EPZ1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'EPZ2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'EP1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'EP2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'R1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. cmpName .eq. 'R2') cycle
        !
        nbmapb = 0
        do kma = 1, nbma_verif
            cellNume = numa_verif(kma)
            iGrel = num_grel(2*(cellNume-1)+1)
            if ((a_un_sens((iGrel-1)*nbCmp+iCmp) .eq. 1) .and. (.not. verif_coef_drz)) cycle
            !
            ient = ptma(cellNume)
            decal = 3+2*nbgdmx+nec*(ient-1)
            dg => desc(decal+1:decal+nec)
            if (.not. exisdg(dg, iCmp)) cycle
            ! si la cmp est nulle, on n'alarme pas :
            !   comptage des cmps presentes pour pouvoir acceder a la valeur
            ico = 0
            do kcmp2 = 1, iCmp
                if (.not. (exisdg(dg, kcmp2))) cycle
                ico = ico+1
            end do
            iad1 = (ient-1)*nbCmp+ico
            !
            if (tsca .eq. 'R') then
                if (verif_coef_drz) then
                    if (zr(jvale-1+iad1) .lt. 0.d0) then
                        if (typmail(cellNume) .eq. MT_TRIA3) then
                            exiq3_coef_drz = exiq3_coef_drz .or. ASTER_TRUE
                        else if (typmail(cellNume) .eq. MT_QUAD4) then
                            exiq4_coef_drz = (exiq4_coef_drz .or. ASTER_TRUE)
                            exiq4_drz_nook = (exiq4_coef_drz) .and. &
                                             (zr(jvale-1+iad1) .gt. -1.0d12) .and. &
                                             (zr(jvale-1+iad1) .lt. -1.0d2)
                        end if
                    else
                        cycle
                    end if
                else
                    if (zr(jvale-1+iad1) .eq. 0.d0) cycle
                end if
            else if (tsca .eq. 'C') then
                if (abs(zc(jvale-1+iad1)) .eq. 0.d0) cycle
            else if (tsca .eq. 'I') then
                if (zi(jvale-1+iad1) .eq. 0) cycle
            else if (tsca .eq. 'L') then
                if (.not. zl(jvale-1+iad1)) cycle
            else if (tsca(1:2) .eq. 'K8') then
                if (zk8(jvale-1+iad1) .eq. ' ') cycle
                if (zk8(jvale-1+iad1) .eq. '&FOZERO') cycle
                if (zk8(jvale-1+iad1) .eq. 'GLOBAL') cycle
                if (zk8(jvale-1+iad1) .eq. 'LOCAL') cycle
                if (zk8(jvale-1+iad1) .eq. 'VENT') cycle
                if (zk8(jvale-1+iad1) .eq. 'LOCAL_PR') cycle
            else if (tsca(1:3) .eq. 'K16') then
                if (zk16(jvale-1+iad1) .eq. ' ') cycle
                if (zk16(jvale-1+iad1) .eq. '&FOZERO') cycle
            else if (tsca(1:3) .eq. 'K24') then
                if (zk24(jvale-1+iad1) .eq. ' ') cycle
                if (zk24(jvale-1+iad1) .eq. '&FOZERO') cycle
            else
                ASSERT(ASTER_FALSE)
            end if
            !
            ! s'il y a un problème
            nbmapb = nbmapb+1
            typeNume = num_grel(2*(cellNume-1)+2)
            if (nbmapb .eq. 1) call jenuno(jexnum('&CATA.TE.NOMTE', typeNume), nomte)
            if (nbmapb .le. 5) list_ma_pb(nbmapb) = cellNume
        end do
        !
        ! Message d'alarme en cas de probleme :
        if (nbmapb .gt. 0) then
            valk(1) = carte
            valk(2) = comment
            valk(3) = nomgd
            valk(4) = cmpName
            valk(5) = nomte
            if (present(non_lin)) then
                if (exiq3_coef_drz .or. exiq4_coef_drz) then
                    call utmess('F', 'CALCULEL_45', nk=5, valk=valk, si=nbmapb)
                end if
            end if
            if (exiq3_coef_drz .and. exiq4_coef_drz) then
                call utmess('A', 'CALCULEL_42', nk=5, valk=valk, si=nbmapb)
                cycle
            else if (exiq3_coef_drz .and. .not. exiq4_coef_drz) then
                call utmess('F', 'CALCULEL_43', nk=5, valk=valk, si=nbmapb)
            else if (.not. exiq3_coef_drz .and. exiq4_drz_nook) then
                call utmess('A', 'CALCULEL_44', nk=5, valk=valk, si=nbmapb)
                cycle
            else
                if (exiq3_coef_drz .or. exiq4_coef_drz) cycle
                if (a_un_sens((iGrel-1)*nbCmp+iCmp) .eq. 1) cycle

                if (cmpName .eq. 'C_METR') then
                    call utmess('F', 'MODELISA10_1')
                else
                    call utmess('A', 'CALCULEL_40', nk=5, valk=valk, si=nbmapb)
                end if
            end if
            do iLiel = 1, min(5, nbmapb)
                valk = ' '
                cellNume = list_ma_pb(iLiel)
                nommai = int_to_char8(cellNume)
                valk(1) = nommai

                call list_grma(mesh, cellNume, 4, lgrma, nbgrma)
                do k1 = 1, min(nbgrma, 3)
                    valk(1+k1) = lgrma(k1)
                end do
                if (nbgrma .gt. 3) valk(5) = '...'
                call utmess('I', 'CALCULEL_41', nk=5, valk=valk)
            end do
        end if
    end do
    AS_DEALLOCATE(vi=numa_verif)
    !
999 continue
    AS_DEALLOCATE(vi=a_un_sens)
    AS_DEALLOCATE(vi=num_grel)
    !
    call jedema()
end subroutine
