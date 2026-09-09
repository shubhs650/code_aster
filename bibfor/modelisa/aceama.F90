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

subroutine aceama(nomu, noma, lmax, nbocc)
    implicit none
#include "jeveux.h"
#include "asterfort/alcart.h"
#include "asterfort/dismoi.h"
#include "asterfort/eulnau.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/wkvect.h"
#include "asterfort/getvid.h"
#include "asterfort/copisd.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
#include "asterfort/chpchd.h"
    integer(kind=8) :: lmax, nbocc
    character(len=8) :: nomu, noma
!     AFFE_CARA_ELEM
!     AFFECTATION DES CARACTERISTIQUES POUR LES ELEMENTS MASSIF
! ----------------------------------------------------------------------
! IN  : NOMU   : NOM UTILISATEUR DE LA COMMANDE
! IN  : NOMA   : NOM DU MAILLAGE
! IN  : LMAX   : LONGUEUR
! IN  : NBOCC  : NOMBRE D'OCCURENCES DU MOT CLE MASSIF
! ----------------------------------------------------------------------
    real(kind=8) :: ang(3), orig(3), angeul(3)
    character(len=19) :: cartma, chorie
    character(len=24) :: tmpnma, tmpvma
!     ------------------------------------------------------------------
!
! --- CONSTRUCTION DES CARTES ET ALLOCATION
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ioc, jdcc, jdls, jdvc, naxe, neul
    integer(kind=8) :: ng, nm, norig, nrep, jdls2, ncham, ndim
    character(len=4) :: fieldDisc
!-----------------------------------------------------------------------
    call jemarq()
    cartma = nomu//'.CARMASSI'
!
! --- Si on a donné CHAM_ORIE alors on ne peut avoir qu'un seul mot-clé facteur MASSIF
!     Déjà vérifié par le superviseur.
    do ioc = 1, nbocc
        call getvid('MASSIF', 'CHAM_ORIE', iocc=ioc, scal=chorie, nbret=ncham)
        if (ncham .ne. 0) then
            ASSERT(nbocc .eq. 1)
        end if
    end do

    if (ncham .ne. 0) then
        call dismoi('TYPE_CHAMP', chorie, 'CHAMP', repk=fieldDisc)
        if (fieldDisc .eq. "ELEM") then
            call chpchd(chorie, "CART", " ", "NON", 'G', cartma, nomu)
        else
            call copisd('CHAMP_GD', 'G', chorie(1:19), cartma(1:19))
        end if

!

    else
!
! --- SINON ON CONSTRUIT LA CARTE AVEC LES INFORMATIONS (GROUP_MA, ANGLE_REP...) DONNEES SOUS MASSIF
        tmpnma = cartma//'.NCMP'
        tmpvma = cartma//'.VALV'
!
        call alcart('G', cartma, noma, 'CAMA_R')
        call jeveuo(tmpnma, 'E', jdcc)
        call jeveuo(tmpvma, 'E', jdvc)
!
        call wkvect('&&TMPMASSIF', 'V V K24', lmax, jdls)
        call wkvect('&&TMPMASSIF2', 'V V K8', lmax, jdls2)
!
!     STOCKAGE DE VALEURS NULLES SUR TOUT LE MAILLAGE
!
        zk8(jdcc) = 'C'
        zk8(jdcc+1) = 'ALPHA'
        zk8(jdcc+2) = 'BETA'
        zk8(jdcc+3) = 'KAPPA'
        zk8(jdcc+4) = 'X'
        zk8(jdcc+5) = 'Y'
        zk8(jdcc+6) = 'Z'
!
        zr(jdvc) = 1.d0
        zr(jdvc+1) = 0.d0
        zr(jdvc+2) = 0.d0
        zr(jdvc+3) = 0.d0
        zr(jdvc+4) = 0.d0
        zr(jdvc+5) = 0.d0
        zr(jdvc+6) = 0.d0
!
        call nocart(cartma, 1, 7)
!
        call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
!
! --- LECTURE DES VALEURS ET AFFECTATION DANS LA CARTE CARTMA
        do ioc = 1, nbocc
            ang(1) = 0.d0
            ang(2) = 0.d0
            ang(3) = 0.d0
            orig(1) = 0.d0
            orig(2) = 0.d0
            orig(3) = 0.d0
            call getvem(noma, 'GROUP_MA', 'MASSIF', 'GROUP_MA', ioc, &
                        lmax, zk24(jdls), ng)
            call getvem(noma, 'MAILLE', 'MASSIF', 'MAILLE', ioc, &
                        lmax, zk8(jdls2), nm)
            call getvr8('MASSIF', 'ANGL_REP', iocc=ioc, nbval=3, vect=ang(1), &
                        nbret=nrep)
            call getvr8('MASSIF', 'ANGL_EULER', iocc=ioc, nbval=3, vect=angeul(1), &
                        nbret=neul)
            call getvr8('MASSIF', 'ANGL_AXE', iocc=ioc, nbval=2, vect=ang(1), &
                        nbret=naxe)
            call getvr8('MASSIF', 'ORIG_AXE', iocc=ioc, nbval=3, vect=orig(1), &
                        nbret=norig)
!
            if (nrep .ne. 0) then
                zr(jdvc) = 1.d0
                zr(jdvc+1) = ang(1)
                zr(jdvc+2) = ang(2)
                zr(jdvc+3) = ang(3)
                zr(jdvc+4) = 0.d0
                zr(jdvc+5) = 0.d0
                zr(jdvc+6) = 0.d0
            else if (neul .ne. 0) then
                call eulnau(angeul, ang)
                zr(jdvc) = 1.d0
                zr(jdvc+1) = ang(1)
                zr(jdvc+2) = ang(2)
                zr(jdvc+3) = ang(3)
                zr(jdvc+4) = 0.d0
                zr(jdvc+5) = 0.d0
                zr(jdvc+6) = 0.d0
            else if (norig .ne. 0) then
                if (ndim .eq. 3 .and. naxe .eq. 0) then
                    call utmess('F', 'MODELISA10_2')
                elseif (ndim .eq. 2 .and. naxe .ne. 0) then
                    call utmess('F', 'MODELISA10_3')
                end if
                zr(jdvc) = -1.d0
                ! pas d'angle en 2D
                zr(jdvc+1) = ang(1)
                zr(jdvc+2) = ang(2)
                zr(jdvc+3) = 0.d0
                zr(jdvc+4) = orig(1)
                zr(jdvc+5) = orig(2)
                zr(jdvc+6) = orig(3)
            else
                ASSERT(.false.)
            end if
!
! ---    "GROUP_MA" = TOUTES LES MAILLES DE LA LISTE DE GROUPES MAILLES
            if (ng .gt. 0) then
                do i = 1, ng
                    call nocart(cartma, 2, 7, groupma=zk24(jdls+i-1))
                end do
            end if
!
! ---    "MAILLE" = TOUTES LES MAILLES DE LA LISTE DE MAILLES
!
            if (nm .gt. 0) then
                call nocart(cartma, 3, 7, mode='NOM', nma=nm, &
                            limano=zk8(jdls2))
            end if
!
        end do
!
        call jedetr('&&TMPMASSIF')
        call jedetr('&&TMPMASSIF2')
        call jedetr(tmpnma)
        call jedetr(tmpvma)
!
    end if
!
    call jedema()
end subroutine
