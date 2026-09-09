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
subroutine bsthpl(plateCara, plateOrie, &
                  jvGeom, nomte, xyzl, &
                  bsigth)
!
    use plate_type
    use plateGeom_module, only: isPlateTria, isPlateQuad
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dxbsig.h"
#include "asterfort/dxefgt.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    type(plateCara_Para), intent(in) :: plateCara
    integer(kind=8), intent(in) :: jvGeom
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: xyzl(3, *)
    real(kind=8), intent(out) :: bsigth(24)
!
! --------------------------------------------------------------------------------------------------
!
!      CALCUL DU BSIGMA POUR LES CONTRAINTES THERMIQUES
!      (I.E. BT*D*ALPHA(T-TREF)) POUR LES ELEMENTS
!                                DE PLAQUE (DKT,DKQ,DST,DSQ,Q4G)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    real(kind=8) :: sigmTher(32)
    real(kind=8) :: pgl(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
    bsigth = zero
    ASSERT(ASTER_FALSE)
    ASSERT(jvGeom .ne. 0)

! - CALCUL DES EFFORTS GENERALISES D'ORIGNIE THERMIQUE AUX POINTS D'INTEGRATION
    call dxefgt(plateCara, plateOrie, &
                sigmTher)

! - CALCUL DE BT*SIGTH
    call dxbsig(plateCara, plateOrie, &
                nomte, 'FORC_NODA', &
                xyzl, pgl, sigmTher, &
                bsigth)
!
end subroutine
