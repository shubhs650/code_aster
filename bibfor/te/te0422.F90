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
subroutine te0422(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/dxefgv.h"
#include "asterfort/dxefro.h"
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
! Elements: DKTG, Q4GG
!
! Options: SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nno, npg
    integer(kind=8) :: jvDisp, jvSief, jvGeom
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: depl(24), efge(32)
    character(len=8), parameter :: fami = 'RIGI'
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'SIEF_ELGA')
    call elrefe_info(fami=fami, nno=nno, npg=npg)
    efge = 0.d0

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Displacements
    call jevech('PDEPLAR', 'L', jvDisp)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of geometry and displacements
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)
    call utpvgl(nno, 6, pgl, zr(jvDisp), depl)

! - CALCUL DES EFFORTS GENERALISES AUX POINTS DE CALCUL
    call dxefgv(plateCara, plateOrie, &
                nomte, option, xyzl, pgl, depl, &
                efge)
!
    call jevech('PCONTRR', 'E', jvSief)
    call dxefro(npg, plateOrie%t2iu, efge, zr(jvSief))
!
end subroutine
