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
subroutine te0225(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
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
! Option: CHAR_MECA_TEMP_R
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: nbProp = 3
    character(len=16), parameter :: propName(nbProp) = (/'E    ', 'NU   ', 'ALPHA'/)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0, deux = 2.d0
    integer(kind=8) :: i, ip, kpg, jvGeom, ivectt, jvMaterc
    integer(kind=8) :: ivf, idfdk, nno, npg, jcoopg, j
    integer(kind=8) :: ipoids, iret1, iret2, iret3, iret4
    real(kind=8) :: tempRefe
    character(len=32) :: elasKeyword
    real(kind=8) :: dfdx(3), r, cour, jac, cosa, sina
    real(kind=8) :: tpg1, tpg2, tpg3, tpg, x3
    real(kind=8) :: h, epsthe, nu, coef, axis
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, nno=nno, npg=npg, &
                     jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk)

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    h = plateCara%thick

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

    call jevech('PVECTUR', 'E', ivectt)

! - TEMPERATURE DE REFERENCE
    call rcvarc(' ', 'TEMP', 'REF', fami, 1, &
                1, tempRefe, iret1)
! - RECUPERATION DE LA NATURE DU MATERIAU DANS PHENOM
    call jevech('PMATERC', 'L', jvMaterc)
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)
!
    if (elasKeyword .eq. 'ELAS') then
        axis = un
        do kpg = 1, npg
            call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                        cour, jac, cosa, sina)
            r = zero
            tpg = zero
            call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                        1, tpg2, iret2)
            call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                        2, tpg1, iret3)
            call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                        3, tpg3, iret4)
            do i = 1, nno
                r = r+zr(jvGeom+2*i-2)*zr(ivf+(kpg-1)*nno+i-1)
            end do
            jac = jac*r

!---- UTILISATION DE 4 POINTS DE GAUSS DANS L'EPAISSEUR
!---- COMME POUR LA LONGUEUR
            do ip = 1, npg
                x3 = zr(jcoopg+ip-1)
                tpg = tpg1*(un-x3**2)+x3*(tpg3*(un+x3)-tpg2*(un-x3))/deux
                call rcvalb('RIGI', 1, 1, '+', zi(jvMaterc), &
                            ' ', 'ELAS', 1, 'TEMP', [tpg], &
                            2, propName, propVale, propCode, 1)
                call rcvalb('RIGI', 1, 1, '+', zi(jvMaterc), &
                            ' ', 'ELAS', 1, 'TEMP', [tpg], &
                            1, propName(3), propVale(3), propCode(3), 0)
                if (((iret1+iret2+iret3+iret4) .ge. 1) .and. (propCode(3) .eq. 0)) then
                    call utmess('F', 'CALCULEL_15')
                else if (propCode(3) .ne. 0) then
                    epsthe = 0.d0
                else
                    epsthe = (tpg-tempRefe)*propVale(3)
                end if
                nu = propVale(2)
                coef = propVale(1)*jac*epsthe*zr(ipoids+ip-1)*(h/deux)
                coef = coef/(un-nu)
                do i = 1, nno
                    j = 3*(i-1)
                    zr(ivectt+j) = zr(ivectt+j)+ &
                                   coef*(axis*zr(ivf+(kpg-1)*nno+i-1)/r-dfdx(i)*sina)
                    zr(ivectt+j+1) = zr(ivectt+j+1)+ &
                                     coef*dfdx(i)*cosa
                    zr(ivectt+j+2) = zr(ivectt+j+2)- &
                                     coef*x3*h/deux*(axis*zr(ivf+(kpg-1)*nno+i-1)*sina/r-dfdx(i))
                end do
            end do
        end do
    else
        call utmess('F', 'ELEMENTS3_49')
    end if
end subroutine
