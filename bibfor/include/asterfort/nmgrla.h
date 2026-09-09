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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine nmgrla(FECell, FEBasis, FEQuad, &
                      nno, npg, ndim, &
                      typmod, option, &
                      compor, carcri, multComp, &
                      BEHInteg, &
                      instam, instap, &
                      dispPrev, dispIncr, &
                      lgpg, sigmPrev, vim, &
                      sigmCurr, vip, &
                      matsym, matuu, vectu, &
                      codret)
        use FE_topo_module
        use FE_quadrature_module
        use FE_basis_module
        use Behaviour_type
        type(FE_Cell), intent(in) :: FECell
        type(FE_Quadrature), intent(in) :: FEQuad
        type(FE_basis), intent(in) :: FEBasis
        integer(kind=8), intent(in) :: nno, npg, ndim
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option
        character(len=16), intent(in) :: compor(COMPOR_SIZE), multComp
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        type(Behaviour_Integ), intent(inout) :: BEHInteg
        real(kind=8), intent(in) :: instam, instap
        real(kind=8), intent(inout) :: dispPrev(ndim*nno), dispIncr(ndim*nno)
        integer(kind=8), intent(in) :: lgpg
        real(kind=8), intent(inout) :: sigmPrev(2*ndim, npg), vim(lgpg, npg)
        real(kind=8), intent(inout) :: sigmCurr(2*ndim, npg), vip(lgpg, npg)
        aster_logical, intent(in) :: matsym
        real(kind=8), intent(inout) :: matuu(*), vectu(ndim*nno)
        integer(kind=8), intent(inout) :: codret
    end subroutine nmgrla
end interface
