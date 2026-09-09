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
subroutine te0226(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getDensity.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/vecma.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_AXIS
! Option: MASS_MECA / M_GAMMA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: kpg, ii, jj, i, j, kd1, kd2, kd3, ij1, ij2, ij3
    real(kind=8) :: dfdx(3), r, rm, rf, rmf, poids, cour, nx, ny, vfi, vfj
    real(kind=8) :: matp(9, 9), matrMass(45), rho
    integer(kind=8) :: ipoids, ivf, idfdk
    integer(kind=8) :: jvGeom, jvMaterc, jvMatr, jvVect, jvAcce
    integer(kind=8) :: nno, npg, nddl, nvec
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk)
    nddl = 3*nno
    nvec = nddl*(nddl+1)/2
    matrMass = 0.d0

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get density
    call jevech('PMATERC', 'L', jvMaterc)
    call getDensity(zi(jvMaterc), rho)
    rm = rho*plateCara%thick
    rf = rho*plateCara%thick**3/12.d0
!
    do kpg = 1, npg
        call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), &
                    dfdx, cour, poids, nx, ny)
        r = 0.d0
        do i = 1, nno
            r = r+zr(jvGeom+2*(i-1))*zr(ivf+(kpg-1)*nno+i-1)
        end do
        poids = poids*r
        rmf = rf*(cour+nx/r)
!
        kd1 = 5
        kd2 = 3
        kd3 = 2
        do i = 1, 3*nno, 3
            kd1 = kd1+3*i-6
            kd2 = kd2+3*i-3
            kd3 = kd3+3*i
            ii = (i+2)/3
            do j = 1, i, 3
                jj = (j+2)/3
                ij1 = kd1+j-2
                ij2 = kd2+j-2
                ij3 = kd3+j-2
                vfi = zr(ivf+(kpg-1)*nno+ii-1)
                vfj = zr(ivf+(kpg-1)*nno+jj-1)
                matrMass(ij1) = matrMass(ij1)+vfi*vfj*poids*rm
                matrMass(ij2) = 0.0d0
                matrMass(ij2+1) = matrMass(ij1)
                matrMass(ij3) = matrMass(ij3)+vfi*vfj*poids*rmf*ny
                matrMass(ij3+1) = matrMass(ij3+1)-vfi*vfj*poids*rmf*nx
                matrMass(ij3+2) = matrMass(ij3+2)+vfi*vfj*poids*rf
            end do
            do j = 1, i-3, 3
                jj = (j+2)/3
                ij1 = kd1+j-2
                ij2 = kd2+j-2
                ij3 = kd3+j-2
                matrMass(ij1+1) = matrMass(ij2)
                matrMass(ij1+2) = matrMass(ij3)
                matrMass(ij2+2) = matrMass(ij3+1)
            end do
        end do
    end do
!
    if (option .eq. 'MASS_MECA') then
        call jevech('PMATUUR', 'E', jvMatr)
        do i = 1, nvec
            zr(jvMatr+i-1) = matrMass(i)
        end do

    else if (option .eq. 'M_GAMMA') then
        call jevech('PACCELR', 'L', jvAcce)
        call jevech('PVECTUR', 'E', jvVect)
        call vecma(matrMass, nvec, matp, nddl)
        call pmavec('ZERO', nddl, matp, zr(jvAcce), zr(jvVect))

    else
        ASSERT(.false.)
    end if
!
end subroutine
