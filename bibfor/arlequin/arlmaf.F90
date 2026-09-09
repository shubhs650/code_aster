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
subroutine arlmaf(mail, mailar, dime, ngrma, ima, &
                  connex, loncum, imail, nummai, cxcumu)
!
!
! ROUTINE ARLEQUIN
!
! CREATION DU MAILLAGE INTERNE ARLEQUIN
! RECOPIE MAILLE
!
! ----------------------------------------------------------------------
!
!
! IN  MAIL   : NOM DU MAILLAGE
! IN  MAILAR : NOM DU PSEUDO-MAILLAGE ARLEQUIN
! IN  DIME   : DIMENSION DU PROBLEME
! IN  NGRMA  : NOM DU GROUPE ARLEQUIN DE LA MAILLE
! IN  IMA    : INDEX MAILLE COURANTE DANS NGRMA
! IN  CONNEX : CONNEXITE DES MAILLES
! IN  LONCUM : LONGUEUR CUMULEE DE CONNEX
! OUT NUMMAI : NUMERO ABSOLU DE LA MAILLE DANS LE MAILLAGE
! IN  IMAIL  : NUMERO MAILLE COURANTE DANS PSEUDO-MAILLAGE
! I/O CXCUMU : LONGUEUR CUMULEE DE LA CONNEXITE
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/arlgrm.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "MeshTypes_type.h"
!
!     ARGUMENTS:
!     ----------
    character(len=8) :: mail, mailar
    integer(kind=8) :: dime
    integer(kind=8) :: connex(*), loncum(*)
    character(len=19) :: ngrma
    integer(kind=8) :: imail
    integer(kind=8) :: ima, nummai, cxcumu
    integer(kind=8) :: itypma
    integer(kind=8) :: cxno(MT_NNOMAX), nbno, ino, cxmax
    integer(kind=8) :: jgcnx, jtypm
    character(len=8) :: nomel, nommai, k8bid
    character(len=24) :: mconn, mtypm
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- NOMS POUR ACCES AU PSEUDO MAILLAGE
!
    mconn = mailar(1:8)//'.CONNEX         '
    mtypm = mailar(1:8)//'.TYPMAIL        '
!
! --- GENERATION NOM DE LA MAILLE
!
    if (imail > 99999) then
        ASSERT(.false.)
    end if
    nomel(1:8) = 'm       '
    call codent(imail, 'D0', nomel(2:8))
!
! --- RECUPERATION INFOS
!
    call arlgrm(mail, ngrma, dime, ima, connex, &
                loncum, nummai, nommai, itypma, nbno, &
                cxno)
!
! --- RECOPIE DU TYPE
!
    call jeveuo(mtypm, 'E', jtypm)
    call jelira(mtypm, 'LONMAX', cxmax, k8bid)
    zi(jtypm+imail-1) = itypma
!
! --- RECOPIE DE LA CONNECTIVITE
!
    cxcumu = cxcumu+nbno
    call jelira(mconn, 'LONT', cxmax, k8bid)
    if (cxcumu > cxmax) then
        ASSERT(.false.)
    end if
    call jeecra(jexnum(mconn, imail), 'LONMAX', nbno, ' ')
    call jeveuo(jexnum(mconn, imail), 'E', jgcnx)
    do ino = 1, nbno
        zi(jgcnx+ino-1) = cxno(ino)
    end do
!
    call jedema()
!
end subroutine
