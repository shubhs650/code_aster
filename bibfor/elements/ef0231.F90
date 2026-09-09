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
subroutine ef0231()
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/ppgan2.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_AXIS
! Option: EFGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 3
    character(len=16), parameter :: propName(nbProp) = (/'E    ', 'NU   ', 'ALPHA'/)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName(nbPara) = (/'TEMP'/)
    real(kind=8) :: paraVale(nbPara)
    real(kind=8) :: e, nu, alpha
    real(kind=8) :: tpg, tgmoy, tgsup, tginf, tempRefe
    real(kind=8) :: x3, epsiKpg(5), c, h, siefTher
    real(kind=8) :: e11, e22, k11, k22, ep11, ep22
    real(kind=8) :: dfdx(3), efgeElno(24)
    real(kind=8) :: jac, r, cosa, sina, cour
    integer(kind=8) :: i, kpg, jvGeom, jvMaterc, jvDisp
    integer(kind=8) :: nno, npg, idfdk, ivf, iret, iret2, iret1, iret3, iret4
    integer(kind=8) :: jcoopg, ip, correc, nbLayer
    integer(kind=8) :: jvEfgeElno, ipoids, jgano
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno, &
                     npg=npg, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, jgano=jgano)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get displacements
    call jevech('PDEPLAR', 'L', jvDisp)

! - Get materials parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Get reference temperature
    call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, &
                1, tempRefe, iret)

! - Get plate properties
    nbLayer = plateCara%nbLayer
    if (nbLayer .le. 0) then
        call utmess('F', 'ELEMENTS_12')
    end if
    if (nbLayer .gt. 30) then
        call utmess('F', 'ELEMENTS3_50')
    end if
    h = plateCara%thick
    correc = nint(plateCara%metric)
!
    efgeElno = 0.d0
    do kpg = 1, npg
        call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                    cour, jac, cosa, sina)
        r = 0.d0
        do i = 1, nno
            r = r+zr(ivf+(kpg-1)*nno+i-1)*zr(jvGeom+2*i-2)
        end do

! ----- Compute strains
        epsiKpg = 0.d0
        do i = 1, nno
            epsiKpg(1) = epsiKpg(1)+dfdx(i)*zr(jvDisp+3*i-3)
            epsiKpg(2) = epsiKpg(2)+dfdx(i)*zr(jvDisp+3*i-2)
            epsiKpg(3) = epsiKpg(3)+dfdx(i)*zr(jvDisp+3*i-1)
            epsiKpg(4) = epsiKpg(4)+zr(ivf+(kpg-1)*nno+i-1)*zr(jvDisp+3*i-3)
            epsiKpg(5) = epsiKpg(5)+zr(ivf+(kpg-1)*nno+i-1)*zr(jvDisp+3*i-1)
        end do
        e11 = epsiKpg(2)*cosa-epsiKpg(1)*sina
        k11 = epsiKpg(3)
        e22 = epsiKpg(4)/r
        k22 = -epsiKpg(5)*sina/r

! ----- Get temperatures
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kpg, 1, tginf, iret1)
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kpg, 2, tgmoy, iret2)
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kpg, 3, tgsup, iret3)
        iret4 = iret1+iret2+iret3
        ASSERT(iret4 .eq. 0 .or. iret4 .eq. 3)

! ----- UTILISATION DE 4 POINTS DE GAUSS DANS L'EPAISSEUR COMME POUR LA LONGUEUR
        do ip = 1, npg
! --------- Mean temperature
            x3 = zr(jcoopg-1+ip)
            if (iret4 .eq. 0) then
                tpg = tgmoy*(1.d0-(x3)**2)+ &
                      tgsup*x3*(1.d0+x3)/2.d0- &
                      tginf*x3*(1.d0-x3)/2.d0
            else
                tpg = r8nnem()
            end if
            x3 = x3*h/2.d0

! --------- Get material parameters
            paraVale(1) = tpg
            call rcvalb('RIGI', 1, 1, '+', &
                        zi(jvMaterc), ' ', 'ELAS', &
                        nbPara, paraName, [paraVale], &
                        2, propName, propVale, &
                        propCode, 1)
            call rcvalb('RIGI', 1, 1, '+', &
                        zi(jvMaterc), ' ', 'ELAS', &
                        nbPara, paraName, [paraVale], &
                        1, propName(3), propVale(3), &
                        propCode(3), 0)
            e = propVale(1)
            nu = propVale(2)
            alpha = propVale(3)
            c = e/(1.d0-nu*nu)
            if (iret4 .eq. 0) then
                if ((propCode(3) .ne. 0) .or. (iret .eq. 1)) then
                    call utmess('F', 'CALCULEL_15')
                else
                    siefTher = (tpg-tempRefe)*alpha*e/(1.d0-nu)
                end if
            else
                siefTher = 0.d0
            end if

! --------- Compute force
            ep11 = (e11+x3*k11)/(1.d0+(correc*x3*cour))
            ep22 = (e22+x3*k22)/(1.d0+(correc*cosa*x3/r))
            efgeElno(6*(kpg-1)+1) = efgeElno(6*(kpg-1)+1)+ &
                                    zr(ipoids-1+ip)*(h/2.d0)*(c*(ep11+nu*ep22)-siefTher)
            efgeElno(6*(kpg-1)+2) = efgeElno(6*(kpg-1)+2)+ &
                                    zr(ipoids-1+ip)*(h/2.d0)*(c*(nu*ep11+ep22)-siefTher)
            efgeElno(6*(kpg-1)+4) = efgeElno(6*(kpg-1)+4)+ &
                                    zr(ipoids-1+ip)*x3*(h/2.d0)*(c*(ep11+nu*ep22)-siefTher)
            efgeElno(6*(kpg-1)+5) = efgeElno(6*(kpg-1)+5)+ &
                                    zr(ipoids-1+ip)*x3*(h/2.d0)*(c*(nu*ep11+ep22)-siefTher)
        end do
        efgeElno(6*(kpg-1)+3) = 0.d0
        efgeElno(6*(kpg-1)+6) = 0.d0
    end do
!
    call jevech('PEFFORR', 'E', jvEfgeElno)
    call ppgan2(jgano, 1, 6, efgeElno, zr(jvEfgeElno))
!
end subroutine
