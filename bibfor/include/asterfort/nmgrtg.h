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
#include "asterf_types.h"
#include "FE_module.h"
!
interface
    subroutine nmgrtg(FEBasis, coorpg, weight, BGSEval, &
        lVect, lMatr, lMatrPred, &
        fPrev, fCurr, dsidep, sigmPrev, &
        sigmCurr, matsym, matuu, vectu)
use FE_basis_module

        type(FE_basis), intent(in) :: FEBasis
        real(kind=8), intent(in) :: dsidep(6, 6), weight, coorpg(3), BGSEval(3, MAX_BS_CG)
        real(kind=8), intent(in) :: sigmCurr(6), sigmPrev(6), fPrev(3, 3), fCurr(3, 3)
        real(kind=8), intent(inout) :: matuu(*), vectu(*)
        aster_logical, intent(in) :: matsym, lVect, lMatr, lMatrPred
    end subroutine nmgrtg
end interface
