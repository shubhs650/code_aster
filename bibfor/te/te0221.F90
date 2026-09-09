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
subroutine te0221(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/moytpg.h"
#include "asterfort/rcvala.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_AXIS
! Option: RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 2
    character(len=16), parameter :: propName(nbProp) = (/'E ', 'NU'/)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName(nbPara) = (/'TEMP'/)
    real(kind=8) ::  paraVale(nbPara)
    real(kind=8), parameter :: eps = 1.d-3
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0, deux = 2.d0, trois = 3.d0, douze = 12.d0
    integer(kind=8) :: ndim, nnos
    real(kind=8) :: gss
    real(kind=8) :: dfdx(3)
    real(kind=8) :: test, test2, nu, h, cosa, sina, cour, r
    real(kind=8) :: coefxx, coefyy, coefxy, coeff1, coeff2
    real(kind=8) :: css, ctt, cts, dss, dts, dtt, bss, btt, bts, vfi, vfj
    real(kind=8) :: c1, c2, c3, cons, cons2, jacp, kappa, correc
    integer(kind=8) :: nno, kpg, npg, jvMatr
    integer(kind=8) :: ii, jj, i, j, ij1, ij2, ij3, kd1, kd2, kd3
    integer(kind=8) :: ipoids, ivf, idfdk, jvGeom, jvMaterc
    integer(kind=8) :: iret
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk)

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    h = plateCara%thick
    kappa = plateCara%shearCoef
    correc = plateCara%metric

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PMATUUR', 'E', jvMatr)
!
    do kpg = 1, npg
        call dfdm1d(nno, zr(ipoids+kpg-1), zr(idfdk+(kpg-1)*nno), zr(jvGeom), dfdx, &
                    cour, jacp, cosa, sina)
        r = zero
        do i = 1, nno
            r = r+zr(jvGeom+2*i-2)*zr(ivf+(kpg-1)*nno+i-1)
        end do

! ----- Get temperature
        call moytpg('RIGI', kpg, 3, '+', paraVale(1), iret)
        test = abs(h*cour/deux)
        if (test .ge. un) correc = zero

! ----- Get elasticity parameters
        call rcvala(zi(jvMaterc), ' ', 'ELAS', &
                    nbPara, paraName, paraVale, &
                    nbProp, propName, propVale, &
                    propCode, 1)
        nu = propVale(2)
        coeff1 = propVale(1)/(un-(nu*nu))
        coeff2 = propVale(1)/(un+nu)
!
!     CALCUL DES COEFFICIENTS RIGIDITE COQUE AXI
!     COUR : COURBURE RS
!     H    : EPAISSEUR
!     COSA : COSINUS ANGLE ALPHA: NORMALE/ HORIZONTALE
!     R    : RAYON VECTEUR
        c1 = -cour*h/(un-(cour*h/deux)**2)
        c2 = -cosa*h/(r**2-(cosa*h/deux)**2)
        if (test .le. eps .or. correc .eq. zero) then
            css = coeff1*h
            bss = zero
            dss = coeff1*h**3/douze
            gss = coeff2*kappa*h/deux
        else
            cons = log((un+(h*cour/deux))/(un-(h*cour/deux)))
            css = ((c1+cons)*cosa/(r*cour**2)+cons/cour)
            bss = -((c1+deux*cons-h*cour)*cosa/(r*cour**3)+cons/(cour**2)-h/cour)*coeff1
            dss = ((cons-h*cour)/cour**3+cosa*(c1+cons*trois-deux*h*cour)/(cour**4*r) &
                   )*coeff1
            gss = css*coeff2*kappa/deux
            css = css*coeff1
        end if
        test2 = abs(h*cosa/(deux*r))
        if (test2 .ge. un) correc = zero
        if (test2 .le. eps .or. correc .eq. zero) then
            ctt = coeff1*h
            btt = zero
            dtt = coeff1*h**3/douze
        else
            cons2 = log((r+(h*cosa/deux))/(r-(h*cosa/deux)))
            c3 = r/cosa
            ctt = (cons2*c3+cour*c3*c3*(r*c2+cons2))*coeff1
            btt = -(c3**3*c2*r*cour+cour*c3*c3*(deux*cons2*c3-h)+c3*(c3*cons2-h))*coeff1
            dtt = (cour*c3**4*r*c2- &
                   h*c3*c3*(un+deux*cour*c3)+ &
                   cons2*c3**3*(un+trois*c3*cour))*coeff1
        end if
        if (abs(cosa) .le. eps .or. abs(cour*r) .le. eps .or. &
            abs(cosa-cour*r) .le. eps .or. correc .eq. zero) then
            cts = coeff1*h
            bts = zero
            dts = coeff1*h**3/douze
        else
            c3 = r/cosa
            cts = ((cour*r**2*cons2-cosa**2*cons/cour)/(cosa*(cour*r-cosa)))*coeff1
            bts = -(-h*(un/cour+c3)+ &
                    cosa*cons/(cour**2*(cosa-cour*r))+ &
                    cons2*cour*cosa*c3**3/(r*cour-cosa))*coeff1
            dts = (-h*(un+c3*cour+c3*c3*cour*cour)/cour**2+ &
                   cons2*cour*c3**3*r/(r*cour-cosa)- &
                   cons/(cour**3*(cour*c3-un)))*coeff1
        end if
        jacp = jacp*r
!
        coefxx = sina*sina
        coefyy = cosa*cosa
        coefxy = -cosa*sina
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
                ij1 = jvMatr+kd1+j-3
                ij2 = jvMatr+kd2+j-3
                ij3 = jvMatr+kd3+j-3
                vfi = zr(ivf+(kpg-1)*nno+ii-1)
                vfj = zr(ivf+(kpg-1)*nno+jj-1)
                zr(ij1) = zr(ij1)+dfdx(ii)*dfdx(jj)*jacp*(coefxx*css+coefyy*gss)
                zr(ij2) = zr(ij2)+dfdx(ii)*dfdx(jj)*jacp*coefxy*(css-gss)
                zr(ij2+1) = zr(ij2+1)+dfdx(ii)*dfdx(jj)*jacp*(coefyy*css+coefxx*gss)
                zr(ij3) = zr(ij3)+dfdx(jj)*jacp*(cosa*gss*vfi+sina*bss*dfdx(ii))
                zr(ij3+1) = zr(ij3+1)-dfdx(jj)*jacp*(cosa*bss*dfdx(ii)-sina*gss*vfi)
                zr(ij3+2) = zr(ij3+2)+jacp*(dss*dfdx(ii)*dfdx(jj)+gss*vfi*vfj)
            end do
            do j = 1, i-3, 3
                jj = (j+2)/3
                ij1 = jvMatr+kd1+j-3
                ij2 = jvMatr+kd2+j-3
                vfi = zr(ivf+(kpg-1)*nno+ii-1)
                vfj = zr(ivf+(kpg-1)*nno+jj-1)
                zr(ij1+1) = zr(ij1+1)+dfdx(ii)*dfdx(jj)*jacp*coefxy*(css-gss)
                zr(ij1+2) = zr(ij1+2)+dfdx(ii)*jacp*(cosa*gss*vfj+sina*bss*dfdx(jj))
                zr(ij2+2) = zr(ij2+2)+dfdx(ii)*jacp*(sina*gss*vfj-cosa*bss*dfdx(jj))
            end do
        end do
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
                ij1 = jvMatr+kd1+j-3
                ij2 = jvMatr+kd2+j-3
                ij3 = jvMatr+kd3+j-3
                vfi = zr(ivf+(kpg-1)*nno+ii-1)
                vfj = zr(ivf+(kpg-1)*nno+jj-1)
                zr(ij1) = zr(ij1)+ &
                          jacp*(ctt*vfi*vfj/(r*r)- &
                                nu*cts*sina*(dfdx(ii)*vfj+dfdx(jj)*vfi)/r)
                zr(ij2) = zr(ij2)+ &
                          jacp*nu*cts*cosa*dfdx(ii)*vfj/r
                zr(ij3) = zr(ij3)+ &
                          jacp*nu*bts*deux*(coefxx*vfi*dfdx(jj)-dfdx(ii)*vfj)/r- &
                          jacp*btt*sina*vfi*vfj/r
                zr(ij3+1) = zr(ij3+1)+ &
                            jacp*nu*bts*deux*coefxy*dfdx(jj)*vfi/r
                zr(ij3+2) = zr(ij3+2)+ &
                            jacp*sina*(dtt*sina*vfi*vfj/(r*r)+ &
                                       nu*dts*(vfi*dfdx(jj)+dfdx(ii)*vfj)/r)
            end do
            do j = 1, i-3, 3
                jj = (j+2)/3
                ij1 = jvMatr+kd1+j-3
                ij2 = jvMatr+kd2+j-3
                vfi = zr(ivf+(kpg-1)*nno+ii-1)
                vfj = zr(ivf+(kpg-1)*nno+jj-1)
                zr(ij1+1) = zr(ij1+1)+ &
                            jacp*nu*cts*cosa*dfdx(jj)*vfi/r
                zr(ij1+2) = zr(ij1+2)+ &
                            jacp*nu*bts*deux*(coefxx*vfj*dfdx(ii)-dfdx(jj)*vfi)/r- &
                            jacp*btt*sina*vfj*vfi/r
                zr(ij2+2) = zr(ij2+2)+ &
                            jacp*nu*bts*deux*coefxy*dfdx(ii)*vfj/r
            end do
        end do
    end do
end subroutine
