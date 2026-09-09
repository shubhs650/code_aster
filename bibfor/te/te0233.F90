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
subroutine te0233(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/getDensity.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_AXIS
! Option: CHAR_MECA_PESA_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: kpg
    real(kind=8) :: dfdx(3), nx, ny, poids, cour, rx
    integer(kind=8) :: nno, npg, i, ivectu, jvPesa
    integer(kind=8) :: ipoids, ivf, idfdk, jvGeom, jvMaterc
    real(kind=8) :: rho
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'CHAR_MECA_PESA_R')

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
    call jevech('PPESANR', 'L', jvPesa)

! - Get density
    call jevech('PMATERC', 'L', jvMaterc)
    call getDensity(zi(jvMaterc), rho)
!
    call jevech('PVECTUR', 'E', ivectu)
    do kpg = 1, npg
        call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                    cour, poids, nx, ny)
        poids = poids*rho*zr(jvPesa)*plateCara%thick
        rx = 0.d0
        do i = 1, nno
            rx = rx+zr(jvGeom+2*i-2)*zr(ivf+(kpg-1)*nno+i-1)
        end do
        poids = poids*rx
        do i = 1, nno
            zr(ivectu+3*i-2) = zr(ivectu+3*i-2)+poids*zr(jvPesa+2)*zr(ivf+(kpg-1)*nno+i-1)
        end do
    end do
end subroutine
