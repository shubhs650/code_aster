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
subroutine te0237(option, nomte)
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
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
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
! Option: SIEF_ELGA / EPSI_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npge = 3
    integer(kind=8), parameter :: nbProp = 3
    character(len=16), parameter :: propName(nbProp) = (/'E    ', 'NU   ', 'ALPHA'/)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName(nbPara) = (/'TEMP'/)
    real(kind=8) :: paraVale(nbPara)
    real(kind=8) :: e, nu, alpha
    real(kind=8) :: tpg, tpgmoy, tpginf, tpgsup, tempRefe
    real(kind=8) :: x3, epsiKpg(5), c1, c2, h, siefTher, niv
    real(kind=8) :: e11, e22, k11, k22, ep11, ep22, ep12, esx3
    real(kind=8) :: dfdx(3)
    real(kind=8) :: jac, r, cosa, sina, cour, correc, zmin, hLayer
    integer(kind=8) :: i, kpg, jvGeom, jvMaterc, jvDisp, jvSief, jvEpsi
    integer(kind=8) :: itab(7)
    integer(kind=8) :: nno, npg, idfdk, ivf, iret, iret1, iret2, iret3, idec, inte
    integer(kind=8) :: iLayer, ipoids, iret4, isp, nbcmp
    integer(kind=8) :: nbLayer
    real(kind=8) :: si11, si12, si22, zic
    real(kind=8), parameter :: ki(3) = (/-1.d0, 0.d0, +1.d0/)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get displacements
    call jevech('PDEPLAR', 'L', jvDisp)

! - Get plate properties
    nbLayer = plateCara%nbLayer
    if (nbLayer .le. 0) then
        call utmess('F', 'ELEMENTS_12')
    end if
    if (nbLayer .gt. 30) then
        call utmess('F', 'ELEMENTS3_50')
    end if
    h = plateCara%thick
    correc = plateCara%metric
    zmin = -h/2.d0
    hLayer = h/nbLayer

! - Get output fields
    if (option .eq. 'EPSI_ELGA') then
        call tecach('OOO', 'PDEFOPG', 'E', iret, nval=7, itab=itab)
        jvEpsi = itab(1)
    else if (option .eq. 'SIEF_ELGA') then
        call jevech('PMATERC', 'L', jvMaterc)
        call tecach('OOO', 'PCONTRR', 'E', iret, nval=7, itab=itab)
        jvSief = itab(1)
        call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, &
                    1, tempRefe, iret)
    else
        ASSERT(.false.)
    end if
!
    nbcmp = itab(2)/itab(3)
    ASSERT(nbcmp .gt. 0)

!
    do iLayer = 1, nbLayer
        do inte = 1, npge
            niv = ki(inte)
            if (inte .eq. 1) then
                zic = zmin+(iLayer-1)*hLayer
            else if (inte .eq. 2) then
                zic = zmin+hLayer/2.d0+(iLayer-1)*hLayer
            else
                zic = zmin+hLayer+(iLayer-1)*hLayer
            end if
            x3 = zic
            do kpg = 1, npg
                idec = nbcmp*(kpg-1)*npge*nbLayer+ &
                       nbcmp*(iLayer-1)*npge+ &
                       nbcmp*(inte-1)
                call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                            cour, jac, cosa, sina)
                r = 0.d0
                do i = 1, nno
                    r = r+zr(ivf+(kpg-1)*nno+i-1)*zr(jvGeom+2*i-2)
                end do

! ------------- Compute strains
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
                esx3 = epsiKpg(5)+epsiKpg(1)*cosa+epsiKpg(2)*sina
                e22 = epsiKpg(4)/r
                k22 = -epsiKpg(5)*sina/r
                ep22 = (e22+x3*k22)/(1.d0+(correc*x3*cosa/r))
                ep11 = (e11+x3*k11)/(1.d0+(correc*x3*cour))
                ep12 = esx3/(1.d0+(correc*x3*cour))

                if (option .eq. 'EPSI_ELGA') then
                    zr(jvEpsi+idec-1+1) = ep11
                    zr(jvEpsi+idec-1+2) = ep22
                    zr(jvEpsi+idec-1+3) = ep12

                else if (option .eq. 'SIEF_ELGA') then
! ---------------- Get temperatures
                    isp = 3*(iLayer-1)
                    call rcvarc(' ', 'TEMP', '+', 'RIGI', kpg, &
                                isp+1, tpginf, iret1)
                    call rcvarc(' ', 'TEMP', '+', 'RIGI', kpg, &
                                isp+2, tpgmoy, iret2)
                    call rcvarc(' ', 'TEMP', '+', 'RIGI', kpg, &
                                isp+3, tpgsup, iret3)
                    iret4 = iret1+iret2+iret3
                    ASSERT(iret4 .eq. 0 .or. iret4 .eq. 3)
                    if (iret4 .eq. 0) then
                        tpg = tpgsup*niv*(1.d0+niv)/2.d0+ &
                              tpgmoy*(1.d0-(niv)**2)- &
                              tpginf*niv*(1.d0-niv)/2.d0
                    else
                        tpg = r8nnem()
                    end if

! ----------------- Get material parameters
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
                    c1 = e/(1.d0+nu)
                    c2 = c1/(1.d0-nu)

! ----------------- Thermal stress
                    if (iret4 .eq. 0) then
                        if ((propCode(3) .ne. 0) .or. (iret .eq. 1)) then
                            call utmess('F', 'CALCULEL_15')
                        else
                            siefTher = (tpg-tempRefe)*alpha*e/(1.d0-nu)
                        end if
                    else
                        siefTher = 0.d0
                    end if

! ----------------- Mechanical stress
                    si11 = c2*(ep11+nu*ep22)-siefTher
                    si22 = c2*(ep22+nu*ep11)-siefTher
                    si12 = c1*ep12
                    zr(jvSief+idec-1+1) = si11
                    zr(jvSief+idec-1+2) = si22
                    zr(jvSief+idec-1+4) = si12
                else
                    ASSERT(.false.)
                end if
            end do
        end do
    end do
!
end subroutine
