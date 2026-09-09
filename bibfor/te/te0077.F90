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
subroutine te0077(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_mass_module
!
    implicit none
#include "FE_module.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/writeMatrix.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'MASS_THER'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: icodre(1)
    character(len=16) :: phenom
    real(kind=8) :: valQP(MAX_QP), cp(1)
    real(kind=8) :: mass(MAX_BS_CG, MAX_BS_CG)
    integer(kind=8) ::  imate, itemps, kp
    character(len=8), parameter :: famiM = "MASS"
!
!-----------------------------------------------------------------------
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, famiM)
    call FEBasis%initCell(FECell)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTR', 'L', itemps)
!
    call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
    if (phenom .eq. 'THER') then
        do kp = 1, FEQuadCell%nbQuadPoints
            call rcvalb(famiM, kp, 1, '+', zi(imate), &
                        ' ', phenom, 1, 'INST', [zr(itemps)], &
                        1, 'RHO_CP', cp, icodre(1), 1)
            valQP(kp) = cp(1)
        end do
    else if (phenom .eq. 'THER_ORTH') then
        do kp = 1, FEQuadCell%nbQuadPoints
            call rcvalb(famiM, kp, 1, '+', zi(imate), &
                        ' ', phenom, 1, 'INST', [zr(itemps)], &
                        1, 'RHO_CP', cp, icodre(1), 1)
            valQP(kp) = cp(1)
        end do
    else
        call utmess('F', 'ELEMENTS2_63')
    end if
!
    call FEMassMatScal(FEQuadCell, FEBasis, mass, valQP)
!
    call writeMatrix("PMATTTR", FEBasis%size, FEBasis%size, ASTER_TRUE, mass)
!
end subroutine
