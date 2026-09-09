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
subroutine op0175()
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterc/getres.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprsd.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/w175af.h"
#include "asterfort/w175ca.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
!
! CALC_FERRAILLAGE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv, n0, numeStoreInit
    integer(kind=8) :: iret, jvPara, ie, nbStore, iStore, numeStore
    character(len=8) :: resultIn, resultOut, model, caraElem
    character(len=16) :: crit, resultType, cmdName
    character(len=19) :: chfer2, chefge
    character(len=19), parameter :: chfer1 = '&&OP0175.CHFER1'
    real(kind=8) :: prec
    integer(kind=8), pointer :: listStore(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)

! - Get input/output results
    call getres(resultOut, resultType, cmdName)
    call getvid(' ', 'RESULTAT', scal=resultIn, nbret=n0)

! - Get all storing index
    call getvr8(' ', 'PRECISION', scal=prec, nbret=ie)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=ie)
    call rsutnu(resultIn, ' ', 0, '&&OP0175.NUME_ORDRE', nbStore, prec, crit, iret)
    ASSERT(iret .eq. 0)
    ASSERT(nbStore .gt. 0)
    call jeveuo('&&OP0175.NUME_ORDRE', 'L', vi=listStore)
!
!
!     -- ON PREND LE MODELE POUR LE 1ER INSTANT :
!     --------------------------------------------
    numeStoreInit = listStore(1)
!
    call rsadpa(resultIn, 'L', 1, 'MODELE', numeStoreInit, 0, sjv=jvPara)
    model = zk8(jvPara)
    ASSERT(model .ne. ' ')
    call getvtx(' ', 'CARA_ELEM', scal=caraElem, nbret=ie)
    ASSERT(caraElem .ne. ' ')
!
!     -- 1. ON CREE LE CHAMP DE DONNEES (CHFER1) :
!     ---------------------------------------------

    call w175af(model, chfer1)
    if (niv .gt. 1) then
        call imprsd('CARTE', chfer1, 6, 'CHFER1=')
    end if
!
!     -- 2. ON APPELLE L'OPTION FERRAILLAGE :
!     -------------------------------------------
    do iStore = 1, nbStore
        numeStore = listStore(iStore)
        call rsexch('F', resultIn, 'EFGE_ELNO', numeStore, chefge, iret)
        call rsexch(' ', resultIn, 'FERR_ELEM', numeStore, chfer2, iret)
        if (resultIn .eq. resultOut) then
            if (iret .eq. 0) then
                call utmess('A', 'CALCULEL_88', si=numeStore)
            end if
        end if
        call w175ca(model, caraElem, chfer1, chefge, chfer2)
        if (niv .gt. 1) then
            call imprsd('CHAMP', chfer2, 6, 'CHFER2=')
        end if
        call rsnoch(resultIn, 'FERR_ELEM', numeStore)
    end do
!
    call jedema()
end subroutine
