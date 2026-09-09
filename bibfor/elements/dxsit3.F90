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
subroutine dxsit3(plateCara, plateOrie, &
                  jvMaterCode, sigma)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dxmate.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/plate_type.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: jvMaterCode
    real(kind=8), intent(out) :: sigma(*)
!
! --------------------------------------------------------------------------------------------------
!
!       CALCUL DES CONTRAINTES VRAIES
!        (==SIGMA_MECA - SIGMA_THER).
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbcmp = 6
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: npg
    integer(kind=8) :: i, j, iLayer, icpg, igauh, kpg, ipgh, nbLayer
    integer(kind=8) :: npgh
    integer(kind=8) :: multic
    real(kind=8) :: zero, epsth(2)
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: h(3, 3), d(4, 4)
    real(kind=8) ::  epais
    aster_logical :: coupmf
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg)
    zero = 0.d0

! - Get parameters of shell
    epais = plateCara%thick

! - Get layers of shell
    nbLayer = plateCara%nbLayer
    if (plateCara%type .eq. PLATE_DKTG) then
        ASSERT(nbLayer .eq. 1)
        npgh = 1
    else
        npgh = 3
    end if
    ASSERT(nbLayer .ge. 1)

! - CARACTERISTIQUES DES MATERIAUX
    call dxmate(plateCara, plateOrie, &
                fami, df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)

! - CALCUL DE LA MATRICE DE HOOKE EN MEMBRANE
    ASSERT(multic .eq. 0)
    do i = 1, 3
        do j = 1, 3
            h(i, j) = dm(i, j)/epais
        end do
    end do

! - passage a la matrice de hooke complete
    d = 0.d0
    d(1:2, 1:2) = h(1:2, 1:2)
    d(1, 4) = h(1, 3)
    d(2, 4) = h(2, 3)
    d(4, 4) = h(3, 3)
    d(4, 1) = h(3, 1)
    d(4, 2) = h(3, 2)
!
    do kpg = 1, npg
        do iLayer = 1, nbLayer
            do igauh = 1, npgh
                icpg = nbcmp*npgh*nbLayer*(kpg-1)+ &
                       nbcmp*npgh*(iLayer-1)+ &
                       nbcmp*(igauh-1)

!         -- INTERPOLATION DE ALPHA EN FONCTION DE LA TEMPERATURE
                ipgh = npgh*(iLayer-1)+igauh
                call verift('RIGI', kpg, ipgh, '+', jvMaterCode, &
                            epsth_=epsth(1))
                epsth(2) = epsth(1)

!           -- CALCUL DES CONTRAINTES VRAIES (==SIGMA_MECA - SIGMA_THER)
                do i = 1, 4
                    do j = 1, 2
                        sigma(icpg+i) = sigma(icpg+i)-epsth(j)*d(i, j)
                    end do
                end do
            end do
        end do
    end do
!
end subroutine
