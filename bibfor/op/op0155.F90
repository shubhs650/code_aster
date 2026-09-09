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
subroutine op0155()
!
    use coorSyst_module, only: setOrieFields
    use result_module, only: rsCopyPara
    implicit none
!
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsutnu.h"
#include "asterfort/w155ce.h"
#include "asterfort/w155ex.h"
#include "asterfort/w155mx.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! POST_CHAMP
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv, nbStore, iret
    character(len=16) :: crit, resultType, cmdName
    character(len=8) :: resultIn, resultOut
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
    call getvid(' ', 'RESULTAT', scal=resultIn, nbret=iret)

! - Get all storing index
    prec = -1.d0
    crit = ' '
    call getvr8(' ', 'PRECISION', scal=prec, nbret=iret)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=iret)
    call rsutnu(resultIn, ' ', 0, '&&OP0155.NUME_ORDRE', nbStore, prec, crit, iret)
    ASSERT(iret .eq. 0)
    ASSERT(nbStore .gt. 0)
    call jeveuo('&&OP0155.NUME_ORDRE', 'L', vi=listStore)

! - Create output result
    call rscrsd('G', resultOut, resultType, nbStore)

! - EXTR_XXXX
    call w155ex(resultOut, resultIn, nbStore, listStore)

! - MIN_MAX_SP
    call w155mx(resultOut, resultIn, nbStore, listStore)

! - COQU_EXCENT
    call w155ce(resultOut, resultIn, nbStore, listStore)

! - Copy parameters
    call rsCopyPara(resultIn, resultOut, nbStore, listStore, copyField_=ASTER_TRUE)

! - RECOPIE DE L'OBJET .REFD
    call refdcp(resultIn, resultOut)
!
    call jedetr('&&OP0155.NUME_ORDRE')
    call jedema()
end subroutine
