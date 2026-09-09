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
subroutine ef0410(nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/vdefro.h"
#include "asterfort/vdxefgeElno.h"
#include "asterfort/vectan.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Compute EFGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npgt = 10
    integer(kind=8) :: jvEfge, jvGeom, lzi
    integer(kind=8) :: nbLayer, nb2, npgsr
    real(kind=8) :: efgeElno(8, 9), matevn(2, 2, npgt)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PEFFORR', 'E', jvEfge)

! - Access objects
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    ASSERT(npgsr .le. 10)

! - Get properties of shell
    call getCara(plateCara, plateOrie)
    nbLayer = plateCara%nbLayer

! - Compute global<=>local transformation
    call compCoorSystCO3D(nomte, jvGeom, &
                          plateCara, plateOrie)

! - Compute
    call vdxefgeElno(plateCara, plateOrie, &
                     nomte, zr(jvGeom), &
                     nbLayer, efgeElno, &
                     matevn)

! - From local to global
    call vdefro(nb2, matevn, efgeElno, zr(jvEfge))
!
end subroutine
