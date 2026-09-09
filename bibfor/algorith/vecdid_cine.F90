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

subroutine vecdid_cine(charge, infoch, numedd, cncine)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
!
!
    character(len=24), intent(in) :: charge, infoch, numedd
    character(len=19), intent(in) :: cncine
!
! --------------------------------------------------------------------------------------------------
!
! Mechanics - Load
!
! complete assembly vector for Dirichlet BC (DIDI) with cinematic loads
!
! --------------------------------------------------------------------------------------------------
!
! In  list_load        : name of datastructure for list of loads
! In  disp_didi        : displacement to compute DIDI loads
! In  vect_asse        : name of vect_asse for DIDI loads
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: neq, numgd, iddes, nec, idprno
    integer(kind=8) :: nchar, icha, iret
    integer(kind=8) :: jafci, nimp, imp, ino, iddl, ieq
    integer(kind=8) :: jvale, jdisp
    character(len=8) :: gd
    character(len=19) :: vect_asse_didi
    character(len=19) :: disp_didi
    integer(kind=8), pointer :: nequ(:) => null()
    integer(kind=8), pointer :: infc(:) => null()
    integer(kind=8), pointer :: vdidicine(:) => null()
    real(kind=8), pointer :: vdidi(:) => null()
    real(kind=8), pointer :: vasse(:) => null()
    character(len=24), pointer :: lcha(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    disp_didi = '&&CNCHAR.DIDI'
    call jeveuo(numedd(1:14)//'.NUME.NEQU', 'L', vi=nequ)
    neq = nequ(1)
    call dismoi('NOM_GD', numedd, 'NUME_DDL', repk=gd)
    call jenonu(jexnom('&CATA.GD.NOMGD', gd), numgd)
    call jeveuo(jexnum('&CATA.GD.DESCRIGD', numgd), 'L', iddes)
    nec = zi(iddes+2)
    call jeveuo(jexnum(numedd(1:14)//'.NUME.PRNO', 1), 'L', idprno)
!
! --- LISTE DES CHARGES
!
    call jelira(charge, 'LONMAX', nchar)
    call jeveuo(charge, 'L', vk24=lcha)
    call jeveuo(infoch, 'L', vi=infc)
!
    call wkvect("&&VEC.DIDI.CINE", 'V V I', neq, vi=vdidicine)
    do icha = 1, nchar
!
! --- VERIF SI CHARGE DE TYPE CINEMATIQUE DIRICHLET DIFFERENTIEL
!
        if (infc(icha+1) .ge. 0 .or. infc(1+3*nchar+2+icha) .eq. 0) then
            cycle
        end if
        call jeveuo(lcha(icha) (1:19)//'.AFCI', 'L', jafci)
        nimp = zi(jafci)
        do imp = 1, nimp
            ino = zi(jafci+3*(imp-1)+1)
            iddl = zi(jafci+3*(imp-1)+2)
            ieq = zi(idprno-1+(nec+2)*(ino-1)+1)+iddl-1
            vdidicine(ieq) = 1
        end do
!
    end do

    call jeveuo(disp_didi//'.VALE', 'L', vr=vdidi)
    call jeveuo(cncine//'.VALE', 'E', vr=vasse)
    do ieq = 1, neq
        vasse(ieq) = vasse(ieq)+vdidicine(ieq)*vdidi(ieq)
    end do

    call jedetr("&&VEC.DIDI.CINE")

    call jedema()
end subroutine
