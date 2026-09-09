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
subroutine te0227(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterc/r8depi.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getDensity.h"
#include "asterfort/jevech.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_AXIS
! Option: MASS_INER
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: dfdx(3), r, rm, poids, cour, nx, ny, yg
    real(kind=8) :: rho, x(3), y(3), xxi, xyi, yyi
    real(kind=8) :: matine(6), volume, depi
    integer(kind=8) :: nno, ipoids, ivf, idfdk, jvGeom, jvMaterc
    integer(kind=8) :: kpg, npg, i, j, lcastr
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    depi = r8depi()
    call elrefe_info(fami='RIGI', nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)
    do i = 1, nno
        x(i) = zr(jvGeom-2+2*i)
        y(i) = zr(jvGeom-1+2*i)
    end do

! - Get density
    call jevech('PMATERC', 'L', jvMaterc)
    call getDensity(zi(jvMaterc), rho)
    rm = rho*plateCara%thick
!
    call jevech('PMASSINE', 'E', lcastr)
!
    volume = 0.d0
    matine = 0.d0

    do kpg = 1, npg
        call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                    cour, poids, nx, ny)
        r = 0.d0
        do i = 1, nno
            r = r+zr(jvGeom+2*(i-1))*zr(ivf+(kpg-1)*nno+i-1)
        end do
        poids = poids*r
        volume = volume+poids
!
        do i = 1, nno
!           --- CDG ---
            zr(lcastr+1) = zr(lcastr+1)+poids*x(i)*zr(ivf+(kpg-1)*nno+i-1)
            zr(lcastr+2) = zr(lcastr+2)+poids*y(i)*zr(ivf+(kpg-1)*nno+i-1)
!           --- INERTIE ---
            xxi = 0.d0
            xyi = 0.d0
            yyi = 0.d0
            do j = 1, nno
                xxi = xxi+x(i)*zr(ivf+(kpg-1)*nno+i-1)*x(j)*zr(ivf+(kpg-1)*nno+j-1)
                xyi = xyi+x(i)*zr(ivf+(kpg-1)*nno+i-1)*y(j)*zr(ivf+(kpg-1)*nno+j-1)
                yyi = yyi+y(i)*zr(ivf+(kpg-1)*nno+i-1)*y(j)*zr(ivf+(kpg-1)*nno+j-1)
            end do
            matine(1) = matine(1)+poids*yyi
            matine(2) = matine(2)+poids*xyi
            matine(3) = matine(3)+poids*xxi
        end do
    end do
!
    yg = zr(lcastr+2)/volume
    zr(lcastr) = depi*volume*rm
    zr(lcastr+3) = yg
    zr(lcastr+1) = 0.d0
    zr(lcastr+2) = 0.d0
!
!    --- ON DONNE LES INERTIES AU CDG ---
    matine(6) = matine(3)*rm*depi
    matine(1) = matine(1)*rm*depi+matine(6)/2.d0-zr(lcastr)*yg*yg
    matine(2) = 0.d0
    matine(3) = matine(1)
    zr(lcastr+4) = matine(1)
    zr(lcastr+5) = matine(3)
    zr(lcastr+6) = matine(6)
    zr(lcastr+7) = matine(2)
    zr(lcastr+8) = matine(4)
    zr(lcastr+9) = matine(5)
!
end subroutine
