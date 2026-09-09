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
subroutine sigimc(materPara, &
                  nbsig, npg, time, &
                  epsini, sigma)
!
    use MaterialPara_type
    use MaterialPara_module
    implicit none
!
#include "asterfort/dmatmc.h"
!
    type(Material_Para), intent(inout) :: materPara
    integer(kind=8), intent(in) :: nbsig, npg
    real(kind=8), intent(in) :: time
    real(kind=8), intent(in) :: epsini(nbsig*npg)
    real(kind=8), intent(out) :: sigma(nbsig*npg)
!
! --------------------------------------------------------------------------------------------------
!
!      SIGIMC   -- CALCUL DES  CONTRAINTES INITIALES
!                  AUX POINTS D'INTEGRATION
!                  POUR LES ELEMENTS ISOPARAMETRIQUES
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: d(36)
    integer(kind=8), parameter :: ksp = 1
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: i, kpg, j
!
! --------------------------------------------------------------------------------------------------
!
    sigma(1:nbsig*npg) = zero

! - Loop on Gauss points
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Compute elasticity matrix
        call dmatmc(materPara, "+", time, &
                    nbsig, d)

! ----- Compute stress
        do i = 1, nbsig
            do j = 1, nbsig
                sigma(i+nbsig*(kpg-1)) = sigma(i+nbsig*(kpg-1))+ &
                                         d(j+(i-1)*nbsig)*epsini(j+nbsig*(kpg-1))
            end do
        end do
    end do
end subroutine
