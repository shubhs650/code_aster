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
subroutine te0006(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
!
    character(len=16), intent(in) :: option, nomte
!
    integer(kind=8) :: ndim, nno, nnos, npg, ino
    integer(kind=8) :: ipoids, ivf, idfde, jgano, nbsig
    integer(kind=8) :: iconti, iconto, ncmp, icmp, ipg, jv_geom
!
    real(kind=8) :: volume, moyenne, poids(100)
    real(kind=8) :: somme(6), xx
    aster_logical :: laxi
!
    character(len=8) :: fami
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    nbsig = nbsigm()
!
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    if (ndim .eq. 3) then
        ncmp = 6
    else if (ndim .eq. 2) then
        ncmp = 4
    end if
    ASSERT(nbsig .eq. ncmp)
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PCONTRR', 'L', iconti)
    call jevech('PSIEFNOR', 'E', iconto)
!
    ASSERT(npg .le. 100)
    ASSERT(ncmp .le. 6)
!
    volume = 0.
    somme(1:6) = 0.
    do ipg = 1, npg
        xx = 0.d0
        do ino = 1, nno
            xx = xx+zr(jv_geom+2*(ino-1)+0)*zr(ivf+(ipg-1)*nno+ino-1)
        end do
        if (ndim .eq. 3) then
            call dfdm3d(nno, ipg, ipoids, idfde, zr(jv_geom), poids(ipg))
        else if (ndim .eq. 2) then
            call dfdm2d(nno, ipg, ipoids, idfde, zr(jv_geom), poids(ipg))
        else
            ASSERT(ASTER_FALSE)
        end if
        if (laxi) poids(ipg) = poids(ipg)*xx
        volume = volume+poids(ipg)
        do icmp = 1, ncmp
            somme(icmp) = somme(icmp)+zr(iconti+(ipg-1)*ncmp+icmp-1)*poids(ipg)
        end do
    end do
!
    do icmp = 1, ncmp
        moyenne = somme(icmp)/volume
        do ipg = 1, npg
            zr(iconto+(ipg-1)*ncmp+icmp-1) = moyenne
        end do
    end do
!
end subroutine
