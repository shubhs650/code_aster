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
subroutine te0232(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getDensity.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_AXIS
! Option: CHAR_MECA_ROTA_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    real(kind=8) :: dfdx(3), nx, ny, poids, cour, rx, ry
    integer(kind=8) :: nno, kpg, npg, i
    integer(kind=8) :: ipoids, ivf, idfdk
    real(kind=8) :: rho
    integer(kind=8) :: jvGeom, jvRota, jvVect, jvMaterc
    real(kind=8) :: rota_speed, rota_axis(3), rota_cent(3)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'CHAR_MECA_ROTA_R')

! - Finite element parameters
    call elrefe_info(fami='RIGI', nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Load
    call jevech('PROTATR', 'L', jvRota)
    rota_speed = zr(jvRota-1+1)
    rota_axis(1) = zr(jvRota-1+2)
    rota_axis(2) = zr(jvRota-1+3)
    rota_axis(3) = zr(jvRota-1+4)
    rota_cent(1) = zr(jvRota-1+5)
    rota_cent(2) = zr(jvRota-1+6)
    rota_cent(3) = zr(jvRota-1+7)

! - Checs
! AXE=Oy et CENTRE=ORIGINE
    if (abs(rota_axis(1)) .gt. r8miem() .or. abs(rota_axis(3)) .gt. r8miem()) then
        call utmess('F', 'CHARGES2_65')
    end if
    if (abs(rota_axis(2)) .le. r8miem()) then
        call utmess('F', 'CHARGES2_65')
    end if
    if (abs(rota_cent(1)) .gt. r8miem() .or. abs(rota_cent(2)) .gt. r8miem() .or. &
        abs(rota_cent(3)) .gt. r8miem()) then
        call utmess('F', 'CHARGES2_66')
    end if

! - Get density
    call jevech('PMATERC', 'L', jvMaterc)
    call getDensity(zi(jvMaterc), rho)

! - OUT fields
    call jevech('PVECTUR', 'E', jvVect)

! - Computation
    do kpg = 1, npg
        call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                    cour, poids, nx, ny)
        poids = poids*rho*rota_speed**2*plateCara%thick
        rx = zero
        ry = zero
        do i = 1, nno
            rx = rx+zr(jvGeom+2*i-2)*zr(ivf+(kpg-1)*nno+i-1)
            ry = ry+zr(jvGeom+2*i-1)*zr(ivf+(kpg-1)*nno+i-1)
        end do
        poids = poids*rx
        do i = 1, nno
            zr(jvVect+3*i-3) = zr(jvVect+3*i-3)+ &
                               poids*rota_axis(2)**2*rx*zr(ivf+(kpg-1)*nno+i-1)
        end do
    end do
end subroutine
