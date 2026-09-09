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
subroutine ef0031(nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/cosiro.h"
#include "asterfort/dxeffi.h"
#include "asterfort/dxefro.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/ppgan2.h"
#include "asterfort/tecach.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DE EFGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbEfgeNd = 8
    integer(kind=8) :: nno, npg, jgano
    integer(kind=8) :: jvEfge, jvGeom, jvSigm
    real(kind=8) :: pgl(3, 3), xyzl(3, 4), effgt(32)
    real(kind=8) :: effint(32)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno, npg=npg, jgano=jgano)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of displacements
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! - Get stress
    call cosiro(plateCara, plateOrie, &
                'PCONTRR', 'L', 'UI', 'G', &
                jvSigm)

! - Compute EFGE_ELGA
    call dxeffi(plateCara, plateOrie, &
                'EFGE_ELGA', nomte, zr(jvSigm), nbEfgeNd, &
                effint)

! - Change parametric => global
    call dxefro(npg, plateOrie%t2iu, effint, effgt)

! - Compute EFGE_ELNO
    call jevech('PEFFORR', 'E', jvEfge)
    call ppgan2(jgano, 1, nbEfgeNd, effgt, zr(jvEfge))
!
end subroutine
