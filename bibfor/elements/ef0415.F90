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
subroutine ef0415(nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterfort/cosiro.h"
#include "asterfort/efcoq3d.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npge = 3
    integer(kind=8) :: jvSiefElga, jvEfgeElno, jvGeom
    integer(kind=8) :: lzi, lzr
    integer(kind=8) :: jvGano
    integer(kind=8) :: nb1, nb2, npgsr, npgsn, nso
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
    if (nomte .eq. 'MEC3QU9H') then
        nso = 4
    else if (nomte .eq. 'MEC3TR7H') then
        nso = 3
    end if

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute global<=>local transformation (intrinsec, for COQUE_3D)
    call compCoorSystCO3D(nomte, jvGeom, &
                          plateCara, plateOrie)

! - Get stress
    call cosiro(plateCara, plateOrie, &
                'PCONTRR', 'L', 'UI', 'G', &
                jvSiefElga)

! - Get GANO matrix
    call jevete('&INEL.'//nomte//'.B', ' ', jvGano)

! - Compute EFGE_ELNO
    call jevech('PEFFORR', 'E', jvEfgeElno)
    call efcoq3d(plateCara, plateOrie, &
                 nomte, nb1, nb2, &
                 npgsn, npgsr, npge, nso, &
                 zr(jvGeom), zr(lzr), &
                 zr(jvSiefElga), zr(jvGano), &
                 zr(jvEfgeElno))
!
end subroutine
