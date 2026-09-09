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
subroutine vdesga(plateCara, plateOrie, &
                  kwgt, nb1, &
                  disp, btild, vectBaseKpg, &
                  hasTemp_, alpha_, tempKpg_, &
                  siefKpg_, epsiKpg_)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/matrc.h"
#include "asterfort/moytem.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: kwgt, nb1
    real(kind=8), intent(in) :: disp(42), btild(5, 42), vectBaseKpg(3, 3)
    aster_logical, optional, intent(in) :: hasTemp_
    real(kind=8), optional, intent(in) :: alpha_, tempKpg_
    real(kind=8), optional, intent(out) :: siefKpg_(6, *), epsiKpg_(6, *)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8) :: matrElas(5, 5), tempMoye
    real(kind=8) :: epsi(5), sigm(5)
    real(kind=8) :: kappa
    integer(kind=8) :: i, k, nbLayer, npg, iret
    aster_logical :: hasTemp
!
! --------------------------------------------------------------------------------------------------
!
    hasTemp = ASTER_FALSE
    if (present(hasTemp_)) then
        hasTemp = hasTemp_
    end if
    call elrefe_info(fami=fami, npg=npg)

! - Get plate parameters
    nbLayer = plateCara%nbLayer
    kappa = plateCara%shearCoef

! - Compute tensor of strains
    epsi = 0.d0
    do i = 1, 5
        do k = 1, 5*nb1+2
            epsi(i) = epsi(i)+btild(i, k)*disp(k)
        end do
    end do

! - Mechanical strains (only to compute stress !)
    if (hasTemp) then
        ASSERT(.not. present(epsiKpg_))
        epsi(1) = epsi(1)-alpha_*tempKpg_
        epsi(2) = epsi(2)-alpha_*tempKpg_
        epsi(3) = epsi(3)
        epsi(4) = epsi(4)
        epsi(5) = epsi(5)
    end if

    if (present(siefKpg_)) then
! ----- Compute mean temperature (on all point and "sous-point" gauss)
        call moytem(fami, npg, 3*nbLayer, '+', tempMoye, iret)

! ----- Compute elastic matrix
        call matrc(plateOrie, vectBaseKpg, tempMoye, kappa, matrElas)

! ----- Compute stress
        sigm = 0.d0
        do i = 1, 5
            do k = 1, 5
                sigm(i) = sigm(i)+matrElas(i, k)*epsi(k)
            end do
        end do
        siefKpg_(1, kwgt) = sigm(1)
        siefKpg_(2, kwgt) = sigm(2)
        siefKpg_(3, kwgt) = 0.d0
        siefKpg_(4, kwgt) = sigm(3)
        siefKpg_(5, kwgt) = sigm(4)
        siefKpg_(6, kwgt) = sigm(5)
    end if
!
    if (present(epsiKpg_)) then
        epsiKpg_(1, kwgt) = epsi(1)
        epsiKpg_(2, kwgt) = epsi(2)
        epsiKpg_(3, kwgt) = 0.d0
        epsiKpg_(4, kwgt) = epsi(3)/2.d0
        epsiKpg_(5, kwgt) = epsi(4)/2.d0
        epsiKpg_(6, kwgt) = epsi(5)/2.d0
    end if
!
end subroutine
