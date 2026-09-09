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
    subroutine nmdlog(FECell, FEBasis, FEQuad, &
                      nno, npg, ndim, &
                      typmod, option, &
                      compor, carcri, multComp, &
                      BEHInteg, &
                      instm, instp, &
                      dispPrev, dispIncr, &
                      lgpg, sigmPrev, vim, &
                      sigmCurr, vip, &
                      matsym, fint, matuu, &
                      codret)
        use FE_topo_module
        use FE_quadrature_module
        use FE_basis_module
        use Behaviour_type
        type(FE_Cell), intent(in) :: FECell
        type(FE_Quadrature), intent(in) :: FEQuad
        type(FE_basis), intent(in) :: FEBasis
        integer(kind=8), intent(in) :: ndim, nno, npg, lgpg
        character(len=8) :: typmod(2)
        character(len=16) :: option
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        character(len=16), intent(in) :: multComp
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        real(kind=8) :: instm
        real(kind=8) :: instp
        aster_logical :: matsym
        real(kind=8) :: dispPrev(*)
        real(kind=8) :: dispIncr(*)
        real(kind=8) :: sigmPrev(2*ndim, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: sigmCurr(2*ndim, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: fint(ndim*nno)
        real(kind=8) :: matuu(*)
        integer(kind=8) :: codret
    end subroutine nmdlog
end interface
