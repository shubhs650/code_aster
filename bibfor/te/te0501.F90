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
subroutine te0501(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_eval_module
!
    implicit none
#include "jeveux.h"
#include "asterfort/ntfcma.h"
#include "asterfort/jevech.h"
#include "asterfort/rcfode.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:   OPTION : 'RIGI_THER_TRANS'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! ......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: kp, imate, ifon(6)
    real(kind=8) :: tpg, alpha, dalpha
    real(kind=8) :: rigi(MAX_BS_CG, MAX_BS_CG)
    real(kind=8) ::  valQPK(3, 3, MAX_QP)
    real(kind=8), pointer :: tempi(:) => null()
    aster_logical :: aniso
! ----------------------------------------------------------------------
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPEI', 'L', vr=tempi)
!
    aniso = ASTER_FALSE
    call ntfcma(' ', zi(imate), aniso, ifon)
!
    valQPK = 0.d0
    do kp = 1, FEQuadCell%nbQuadPoints
        tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuadCell%points_param(1:3, kp))
        call rcfode(ifon(2), tpg, alpha, dalpha)
        valQPK(1, 1, kp) = alpha
        valQPK(2, 2, kp) = alpha
        valQPK(3, 3, kp) = alpha
    end do
!
    call FEStiffJacoScal(FEQuadCell, FEBasis, valQPK, rigi)
    call writeMatrix("PMATTTR", FEBasis%size, FEBasis%size, ASTER_TRUE, rigi)
!
end subroutine
