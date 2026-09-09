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
subroutine te0309(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_mass_module
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
!
    character(len=16) :: nomte, option
!.......................................................................
!
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES DE FLUX FLUIDE EN MECANIQUE
!          ELEMENTS ISOPARAMETRIQUES 1D
!
!          OPTION : 'FLUX_FLUI_X 'OU 'FLUX_FLUI_Y 'OU 'FLUX_FLUI_Z '
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    type(FE_Skin) :: FESkin
    type(FE_Quadrature) :: FEQuad
    type(FE_basis) :: FEBasis
!
    real(kind=8) :: normal(3)
    integer(kind=8) :: index, kp
    real(kind=8) :: mass(MAX_BS_CG, MAX_BS_CG)
    real(kind=8) :: valQP(MAX_QP)
!-----------------------------------------------------------------------
    call FESkin%init()
    call FEQuad%initFace(FESkin, "RIGI")
    call FEBasis%initFace(FESkin)
!
    if (option(11:11) .eq. 'X') then
        index = 1
    elseif (option(11:11) .eq. 'Y') then
        index = 2
    elseif (option(11:11) .eq. 'Z') then
        index = 3
    else
        ASSERT(ASTER_FALSE)
    end if
!
    do kp = 1, FEQuad%nbQuadPoints
        normal = FESkin%normal(FEQuad%points_param(1:2, kp))
        valQP(kp) = normal(index)
    end do
!
    call FEMassMatScal(FEQuad, FEBasis, mass, valQP)
!
    call writeMatrix("PMATTTR", FEBasis%size, FEBasis%size, ASTER_TRUE, mass)
!
end subroutine
