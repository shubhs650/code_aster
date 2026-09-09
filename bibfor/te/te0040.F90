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
subroutine te0040(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cosiro.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/tecach.h"
#include "asterfort/elno_coq3d.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE 3D
!     OPTIONS : EPSI_ELNO
!               SIEF_ELNO
!               SIGM_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jvCompor, jvFieldPg
    integer(kind=8) :: jvFieldNo, jvGeom
    integer(kind=8) :: lzi, lzr, jvGano
    integer(kind=8) :: nbLayer, nso, nb1, nb2, npgsr, npgsn
    integer(kind=8) :: iret
    aster_logical :: lgreen
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    lgreen = .false.
    if (nomte .eq. 'MEC3QU9H') then
        nso = 4
    else if (nomte .eq. 'MEC3TR7H') then
        nso = 3
    end if

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    nbLayer = plateCara%nbLayer

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute global<=>local transformation
    call compCoorSystCO3D(nomte, jvGeom, &
                          plateCara, plateOrie)

! - Prepare fields
    if (option .eq. 'EPSI_ELNO') then
        call jevech('PDEFOPG', 'L', jvFieldPg)
        call jevech('PDEFONO', 'E', jvFieldNo)

    else if ((option .eq. 'SIEF_ELNO') .or. &
             (option .eq. 'SIGM_ELNO')) then
        call cosiro(plateCara, plateOrie, &
                    'PCONTRR', 'L', 'UI', 'G', &
                    jvFieldPg)
        call jevech('PSIEFNOR', 'E', jvFieldNo)
        call tecach('ONO', 'PCOMPOR', 'L', iret, iad=jvCompor)
        if (jvCompor .ne. 0) then
            if (zk16(jvCompor+2) .eq. 'GROT_GDEP') then
                lgreen = .true.
            end if
        end if
    else
        ASSERT(.false.)
    end if

! - Get GANO matrix
    call jevete('&INEL.'//nomte//'.B', ' ', jvGano)

! - Compute
    call elno_coq3d(plateCara, plateOrie, &
                    lgreen, option, nomte, &
                    nb2, npgsn, nso, nbLayer, &
                    zr(jvGano), zr(jvFieldPg), &
                    zr(jvFieldNo))
!
end subroutine
