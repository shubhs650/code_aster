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
! aslint: disable=W0413
!
subroutine te0397(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_AXIS
! Option: CHAR_MECA_PRES_*
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPara = 4
    character(len=8), parameter :: paraName(nbPara) = (/"X   ", "Y   ", "Z   ", "INST"/)
    real(kind=8) :: paraVale(nbPara)
    integer(kind=8) :: nno, nddl, kpg, npg, ipoids, ivf, idfdk, jvGeom
    integer(kind=8) :: ivectu, i, l, ipres, ier, iadzi, iazk24, itemps
    real(kind=8) ::  poids, r, fx, fy, f3, nx, ny, cour, dfdx(3), pr
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk)
!
!
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PVECTUR', 'E', ivectu)
    nddl = 3
!
    if (option .eq. 'CHAR_MECA_PRES_R') then
        call jevech('PPRESSR', 'L', ipres)
        do kpg = 1, npg
            call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                        cour, poids, nx, ny)
            r = 0.d0
            fx = 0.d0
            fy = 0.d0
            do i = 1, nno
                l = (kpg-1)*nno+i
!-----------------------------------------------------
!              LE SIGNE MOINS CORRESPOND A LA CONVENTION :
!                 UNE PRESSION POSITIVE PROVOQUE UN GONFLEMENT
!-----------------------------------------------------
                f3 = -zr(ipres+i-1)
                fx = fx+nx*f3*zr(ivf+l-1)
                fy = fy+ny*f3*zr(ivf+l-1)
                r = r+zr(jvGeom+2*(i-1))*zr(ivf+l-1)
            end do
            poids = poids*r
            do i = 1, nno
                l = (kpg-1)*nno+i
                zr(ivectu+nddl*(i-1)) = zr(ivectu+nddl*(i-1))+fx*zr(ivf+l-1)*poids
                zr(ivectu+nddl*(i-1)+1) = zr(ivectu+nddl*(i-1)+1)+fy*zr(ivf+l-1)*poids
            end do
        end do
    else if (option .eq. 'CHAR_MECA_PRES_F') then
        call jevech('PPRESSF', 'L', ipres)
        call jevech('PINSTR', 'L', itemps)
        paraVale(4) = zr(itemps)
        do i = 0, nno-1
            paraVale(1) = zr(jvGeom+3*i)
            paraVale(2) = zr(jvGeom+3*i+1)
            paraVale(3) = zr(jvGeom+3*i+2)
            call fointe('FM', zk8(ipres), nbPara, paraName, paraVale, &
                        pr, ier)
            if (pr .ne. 0.d0) then
                call tecael(iadzi, iazk24)
                call utmess('F', 'ELEMENTS4_92', si=zi(iadzi-1+1))
            end if
        end do

    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
