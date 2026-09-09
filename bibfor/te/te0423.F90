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
subroutine te0423(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    implicit none
!
#include "asterfort/dxbsig.h"
#include "asterfort/dxefg2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKTG
!
! Options: CHAR_MECA_TEMP_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nno, i, jvGeom, jvVect
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: forcNoda(24), sigt(32)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of geometry
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! - CALCUL DES EFFORTS GENERALISES D'ORIGINE THERMIQUE AUX POINTS D'INTEGRATION
    call dxefg2(plateCara, plateOrie, &
                pgl, sigt)

! - CALCUL DES EFFORTS INTERNES D'ORIGINE THERMIQUE
    call dxbsig(plateCara, plateOrie, &
                nomte, 'FORC_NODA', &
                xyzl, pgl, sigt, &
                forcNoda)

! - AFFECTATION DU VECTEUR DES FORCES ELEMENTAIRES
    call jevech('PVECTUR', 'E', jvVect)
    do i = 1, nno*6
        zr(jvVect+i-1) = forcNoda(i)
    end do
!
end subroutine
