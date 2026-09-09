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
subroutine te0506(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_rhs_module
    use FE_eval_module
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN THERMIQUE
!          CORRESPONDANT AU FLUX NON-LINEAIRE
!
!         OPTION : 'CHAR_THER_FLUTNL'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!----------------------------------------------------------------------
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/writeVector.h"
#include "asterfort/foderi.h"
#include "FE_module.h"
!
    type(FE_Skin) :: FESkin
    type(FE_Quadrature) :: FEQuad
    type(FE_Basis) :: FEBasis
!
    character(len=16) :: option, nomte
    real(kind=8) :: alpha, dalpha, tpg
    integer(kind=8) :: iflux, kp
    character(len=8) :: coef
    real(kind=8) :: rhs(MAX_BS_CG), valQP(MAX_QP)
    real(kind=8), pointer :: tempi(:) => null()
!
    call FESkin%init()
    call FEQuad%initFace(FESkin, "RIGI")
    call FEBasis%initFace(FESkin)
!
    call jevech('PTEMPEI', 'L', vr=tempi)
    call jevech('PFLUXNL', 'L', iflux)
!
    coef = zk8(iflux)
    if (coef(1:7) .eq. '&FOZERO') goto 999
!
    valQP = 0.d0
    do kp = 1, FEQuad%nbQuadPoints
        tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuad%points_param(1:3, kp))
!
        call foderi(coef, tpg, alpha, dalpha)
!
        valQP(kp) = alpha-dalpha*tpg
    end do
!
    call FeMakeRhsScal(FEQuad, FEBasis, valQP, rhs)
    call writeVector("PRESIDU", FEBasis%size, rhs)
!
999 continue
!
end subroutine
