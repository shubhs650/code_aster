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
subroutine ef0033(nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/coqrep.h"
#include "asterfort/dxefgv.h"
#include "asterfort/dxefro.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DE EFGE_ELNO EN LINEAIRE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nno, npg
    integer(kind=8) :: jvDisp, jvEfge, jvGeom
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: depl(24), efgeElno(32)
    character(len=8), parameter :: fami = 'NOEU'
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, nno=nno, npg=npg)
    efgeElno = 0.d0

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

! - CALCUL DES EFFORTS GENERALISES VRAIS AUX POINTS DE CALCUL
    call dxefgv(plateCara, plateOrie, &
                nomte, 'EFGE_ELNO', xyzl, pgl, depl, efgeElno)

! - PASSAGE DES EFFORTS GENERALISES DU REPERE INTRINSEQUE A L'ELEMENT AU REPERE LOCAL DE LA COQUE
    call jevech('PEFFORR', 'E', jvEfge)
    call dxefro(nno, plateOrie%t2iu, efgeElno, zr(jvEfge))
!
end subroutine
