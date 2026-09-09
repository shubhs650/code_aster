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
subroutine te0416(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterfort/cosiro.h"
#include "asterfort/forngr.h"
#include "asterfort/fornpd.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D - FORC_NODA / REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) ::  iret, icompo, jvGeom
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

! - Compute vector
    call tecach('ONO', 'PCOMPOR', 'L', iret, iad=icompo)
    if (icompo .eq. 0) then
        call fornpd(plateCara, plateOrie, &
                    option, nomte)
    else
        call jevech('PCOMPOR', 'L', vk16=compor)
        if (compor(DEFO) .eq. 'GROT_GDEP') then
            call forngr(plateCara, plateOrie, &
                        option, nomte)
        else if (compor(DEFO) (1:5) .eq. 'PETIT') then
            call fornpd(plateCara, plateOrie, &
                        option, nomte)
        else
            call utmess('F', 'ELEMENTS3_93', sk=compor(DEFO))
        end if
    end if
!
end subroutine
