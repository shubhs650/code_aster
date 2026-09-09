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
subroutine te0410(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterfort/jevech.h"
#include "asterfort/postcoq3d.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_3D
!
! Options: DEGE_ELGA, DEGE_ELNO, EPSI_ELGA, SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)  :: jvGeom
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute global<=>local transformation
    call compCoorSystCO3D(nomte, jvGeom, &
                          plateCara, plateOrie)

! - Compute DEGE_ELGA, DEGE_ELNO, EPSI_ELGA, SIEF_ELGA
    call postcoq3d(plateCara, plateOrie, &
                   option, nomte, plateCara%nbLayer)
!
end subroutine
