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
! => real zero (affect here)
!
subroutine te0234(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterf_types.h"
#include "asterfort/defgen.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/effi.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/moytpg.h"
#include "asterfort/rcvala.h"
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
! Option: FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: eps = 1.d-3
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0, deux = 2.d0
    integer(kind=8), parameter :: npge = 3
    integer(kind=8), parameter :: nbProp = 2
    character(len=16), parameter :: propName(nbProp) = (/'E ', 'NU'/)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName(nbPara) = (/'TEMP'/)
    real(kind=8) ::  paraVale(nbPara)
    integer(kind=8) :: nno, npg, nbsp, itab(7), nbLayer
    integer(kind=8) :: jvGeom, jvMaterc, jvSief, jvDisp, jvVect
    integer(kind=8) :: iLayer, inte, kpki, k1, i, iret, kpg
    real(kind=8) :: cisail, zic, coef, rhos, rhot, epsx3, gsx3, sgmsx3
    real(kind=8) :: dfdx(3)
    real(kind=8) :: test, test2, nu, h, cosa, sina, cour, r
    real(kind=8) :: jacp, kappa, correc
    real(kind=8) :: eps2d(4), sigtdi(5), sigmtd(5)
    real(kind=8) :: x3
    integer(kind=8) :: ipoids, ivf, idfdk
    aster_logical :: testl1, testl2
    real(kind=8) :: zmin, hLayer
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdk)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get plate properties
    nbLayer = plateCara%nbLayer
    if (nbLayer .le. 0) then
        call utmess('F', 'ELEMENTS_12')
    end if
    if (nbLayer .gt. 30) then
        call utmess('F', 'ELEMENTS3_50')
    end if
    h = plateCara%thick
    kappa = plateCara%shearCoef
    correc = plateCara%metric
    zmin = -h/2.d0
    hLayer = h/nbLayer

! - Material properties
    call jevech('PMATERC', 'L', jvMaterc)

! - Get stress
    call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, itab=itab)
    jvSief = itab(1)
    nbsp = itab(7)
    if (nbsp .ne. npge*nbLayer) then
        call utmess('F', 'ELEMENTS_4')
    end if

! - Get displacements
    call jevech('PDEPLAR', 'L', jvDisp)

! - INITIALISATION DU VECTEUR FORCE INTERNE
    call jevech('PVECTUR', 'E', jvVect)
    do i = 1, 3*nno
        zr(jvVect+i-1) = 0.d0
    end do
!
    kpki = 0
    do kpg = 1, npg
        call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                    cour, jacp, cosa, sina)
        r = zero
        do i = 1, nno
            r = r+zr(jvGeom+2*i-2)*zr(ivf+(kpg-1)*nno+i-1)
        end do
        jacp = jacp*r

! ----- Get temperature
        call moytpg('RIGI', kpg, npge, '+', paraVale(1), iret)

! ----- Get elasticity parameters
        call rcvala(zi(jvMaterc), ' ', 'ELAS', &
                    nbPara, paraName, paraVale, &
                    nbProp, propName, propVale, &
                    propCode, 1)
        nu = propVale(2)
        cisail = propVale(1)/(un+nu)

! ----- MODI_METRIQUE ?
        test = abs(h*cour/deux)
        if (test .ge. un) correc = zero
        test2 = abs(h*cosa/(deux*r))
        if (test2 .ge. un) correc = zero
!
        testl1 = (test .le. eps .or. correc .eq. zero)
        testl2 = (test2 .le. eps .or. correc .eq. zero .or. &
                  abs(cosa) .le. eps .or. abs(cour*r) .le. eps .or. &
                  abs(cosa-cour*r) .le. eps)
!
        sigmtd = zero
        do iLayer = 1, nbLayer
            do inte = 1, npge
                if (inte .eq. 1) then
                    zic = zmin+(iLayer-1)*hLayer
                    coef = 1.d0/3.d0
                else if (inte .eq. 2) then
                    zic = zmin+hLayer/2.d0+(iLayer-1)*hLayer
                    coef = 4.d0/3.d0
                else
                    zic = zmin+hLayer+(iLayer-1)*hLayer
                    coef = 1.d0/3.d0
                end if
                x3 = zic

! ------------- Apply MODI_METRIQUE (or not !)
                if (testl1) then
                    rhos = 1.d0
                else
                    rhos = 1.d0+x3*cour
                end if
                if (testl2) then
                    rhot = 1.d0
                else
                    rhot = 1.d0+x3*cosa/r
                end if

! ------------- CALCUL DES COMPOSANTES DE DEFORMATIONS TRIDIMENSIONNELLES EPSSS, EPSTT, EPSSX3
                call defgen(testl1, testl2, nno, r, x3, &
                            sina, cosa, cour, zr(ivf+(kpg-1)*nno), dfdx, &
                            zr(jvDisp), eps2d, epsx3)

! ------------- CONSTRUCTION DE LA DEFORMATION GSX3 ET DE LA CONTRAINTE SGMSX3
                gsx3 = 2.d0*epsx3
                sgmsx3 = cisail*kappa*gsx3/2.d0

! ------------- CALCUL DES CONTRAINTES TILDE
                kpki = kpki+1
                k1 = 4*(kpki-1)
                sigtdi(1) = zr(jvSief-1+k1+1)/rhos
                sigtdi(2) = x3*zr(jvSief-1+k1+1)/rhos
                sigtdi(3) = zr(jvSief-1+k1+2)/rhot
                sigtdi(4) = x3*zr(jvSief-1+k1+2)/rhot
                sigtdi(5) = sgmsx3/rhos
                do i = 1, 5
                    sigmtd(i) = sigmtd(i)+sigtdi(i)*0.5d0*hLayer*coef
                end do
            end do
        end do

! ----- CALCUL DES EFFORTS INTERIEURS
        call effi(nomte, sigmtd, zr(ivf+(kpg-1)*nno), dfdx, jacp, &
                  sina, cosa, r, zr(jvVect))
    end do
!
end subroutine
