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
subroutine arlmai(mail, mailar, ndim, nom1, nom2, &
                  tabcor, nbma1, nbma2)
!
!
! ROUTINE ARLEQUIN
!
! CREATION DU PSEUDO-MAILLAGE
!
! ----------------------------------------------------------------------
!
!
! IN  MAIL   : NOM DU MAILLAGE
! IN  ULPMAI : UNITE LOGIQUE POUR L'IMPRESSION DU PSEUDO-MAILLAGE
! OUT MAILAR : NOM DU PSEUDO-MAILLAGE
! IN  NDIM   : DIMENSION DU PROBLEME
! IN  NOM1   : NOM DE LA SD DE STOCKAGE PREMIER GROUPE
! IN  NOM2   : NOM DE LA SD DE STOCKAGE SECOND GROUPE
! OUT TABCOR : TABLEAU DE CORRESPONDANCE
!            POUR CHAQUE NOUVEAU NUMERO ABSOLU DANS MAILAR
!             -> ANCIEN NUMERO ABSOLU DANS MAIL
!             -> SI NEGATIF, LA NOUVELLE MAILLE EST ISSUE D'UNE
!                DECOUPE DE LA MAILLE DE NUMERO ABSOLU ABS(NUM) DANS
!                MAIL
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/arlini.h"
#include "asterfort/arlmaf.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
!
!
!     ARGUMENTS:
!     ----------
    character(len=8) :: mail, mailar
    integer(kind=8) :: ndim, nbma1, nbma2
    character(len=10) :: nom1, nom2
    character(len=24) :: tabcor
!
!
    integer(kind=8) :: nbnomx
    parameter(nbnomx=MT_NNOMAX)
    integer(kind=8) :: nmain1, nmain2
    integer(kind=8) :: nnoin1, nnoin2, cxcumu
    character(len=24) :: k8bid
    integer(kind=8) :: icpl, i
    integer(kind=8) :: numma1, numma2
    integer(kind=8) :: imail
    integer(kind=8) :: nbno, nbmat, nctot
    integer(kind=8) :: iret
    character(len=24) :: maidim, cooval
    integer(kind=8) :: jdime, jcoor
    integer(kind=8) :: jcooro, jconxo, jcumuo
    character(len=19) :: ngrm1, ngrm2
    integer(kind=8) :: jtabco
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- NOM DES SDS TEMPORAIRES
!
    mailar = '&&ARL.MA'
    tabcor = '&&ARLMAI.TABCOR'
!
! --- DESTRUCTION DU PSEUDO-MAILLAGE S'IL EXISTE
!
    call exisd('MAILLAGE', mailar, iret)
    if (iret .ne. 0) then
        call detrsd('MAILLAGE', mailar)
    end if
!
! --- INITIALISATIONS
!
    ngrm1 = nom1(1:10)//'.GROUPEMA'
    ngrm2 = nom2(1:10)//'.GROUPEMA'
    call jelira(ngrm1, 'LONMAX', nbma1, k8bid)
    call jelira(ngrm2, 'LONMAX', nbma2, k8bid)
!
! --- MAILLAGE INITIAL
!
    call jeveuo(mail(1:8)//'.COORDO    .VALE', 'L', jcooro)
    call jeveuo(mail(1:8)//'.CONNEX         ', 'L', jconxo)
    call jeveuo(jexatr(mail(1:8)//'.CONNEX         ', 'LONCUM'), 'L', jcumuo)
    call dismoi('NB_NO_MAILLA', mail, 'MAILLAGE', nbno, k8bid, &
                'F', iret)
!
! --- ESTIMATION GROSSIERE DU NOMBRE D'ELEMENTS TOTAL
!
    nmain1 = nbma1
    nmain2 = nbma2
    nbmat = nmain1+nmain2
!
! --- ESTIMATION GROSSIERE DE LA LONGUEUR CUMULEE DES CONNECTIVITES
!
    nnoin1 = nbma1*nbnomx
    nnoin2 = nbma2*2
    nctot = nnoin1+nnoin2
!
! --- ESTIMATION DU NOMBRE DE NOEUDS TOTAL
! --- CREATIONS DES OBJETS DE BASE POUR LE PSEUDO MAILLAGE
!
    call arlini(mailar, 'V', ndim, nbno, nbmat, &
                nctot)
!
! --- ACCES AUX OBJETS DU PSEUDO MAILLAGE
!
    maidim = mailar(1:8)//'.DIME           '
    cooval = mailar(1:8)//'.COORDO    .VALE'
    call jeveuo(maidim, 'E', jdime)
    call jeveuo(cooval, 'E', jcoor)
!
! --- RECOPIE DES COORDONNEES DES ANCIENS NOEUDS
!
    do i = 1, 3*nbno
        zr(jcoor+i-1) = zr(jcooro+i-1)
    end do
!
! --- TABLEAU DE CORRESPONDANCE NOUV. MAILLES -> ANCIENNES MAILLES
!
    call wkvect(tabcor, 'V V I', nbmat, jtabco)
!
! --- INDICES POUR REMPLIR LES TABLEAUX
!
    imail = 0
    cxcumu = 0
!
! --- CREATION DE LA MAILLE DANS LE PSEUDO-MAILLAGE
!
    do icpl = 1, nmain1
!
! ----- MAILLAGE 1
!
        imail = imail+1
        call arlmaf(mail, mailar, ndim, ngrm1, icpl, &
                    zi(jconxo), zi(jcumuo), imail, numma1, cxcumu)
        zi(jtabco+imail-1) = numma1
    end do
    do icpl = 1, nmain2
!
! ----- MAILLAGE 2
!
        imail = imail+1
        call arlmaf(mail, mailar, ndim, ngrm2, icpl, &
                    zi(jconxo), zi(jcumuo), imail, numma2, cxcumu)
        zi(jtabco+imail-1) = numma2
    end do
!
! --- VRAI NOMBRE D'ELEMENTS
!
    zi(jdime-1+3) = imail
!
    call jedema()
!
end subroutine
