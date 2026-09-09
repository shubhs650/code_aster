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
subroutine te0494(nomopt, nomte)
!
    use HHO_type
    use HHO_size_module, only: hhoTherDofs
    use HHO_init_module, only: hhoInfoInitCellAndFace
    use HHO_basis_module
    use FE_algebra_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/HHO_basis_module.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO - Generic
!  Option: HHO_PRECALC_BS
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: nomte, nomopt
!
! -- Local variables
!
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_basis_cell) :: hhoBasisCell
    type(HHO_basis_face) :: hhoBasisFace
    real(kind=8) :: basis(6*MAX_FACE_COEF+MAX_CELL_COEF)
    integer(kind=8) :: dec, iFace, nbCoeffFace, nbCoeffCell
    integer(kind=8) :: max_deg_cell, max_deg_face, face_size, cell_size
!
    ASSERT(nomopt .eq. 'HHO_PRECALC_BS')
!
! --- Retrieve HHO informations
!
    call hhoInfoInitCellAndFace(hhoCell, hhoData)
!
    call hhoGetMaxDegree(max_deg_cell, max_deg_face)
!
    dec = 1
    do iFace = 1, hhoCell%nbfaces
        call hhoBasisFace%initialize(hhoCell%faces(iFace))
        face_size = hhoBasisFace%BSSize(0, max_deg_face)
        nbCoeffFace = face_size*(face_size+1)/2
        call dcopy_1(nbCoeffFace, hhoBasisFace%coeff_mono, basis(dec))
        dec = dec+nbCoeffFace
    end do
!
    call hhoBasisCell%initialize(hhoCell)
    cell_size = hhoBasisCell%BSSize(0, max_deg_cell)
    nbCoeffCell = cell_size*(cell_size+1)/2
    call dcopy_1(nbCoeffCell, hhoBasisCell%coeff_mono, basis(dec))
    dec = dec+nbCoeffCell
!
! -- Save - the name is not PCHHOBS because reading this field in basis
!
    call writeVector('PCHHOBO', dec-1, basis)
!
end subroutine
