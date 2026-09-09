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
subroutine te0433(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystGrid
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cargri.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getDensity.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmgrib.h"
#include "asterfort/rcvalb.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: GRILLE_MEMBRANE / GRILLE_EXCENTRE
!
! Elements: EPOT_ELEM
!           EPSI_ELGA
!           SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: nbProp = 1
    integer(kind=8) :: propCode(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'E'/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: nddl, nno, npg, i, j, n, kpg
    integer(kind=8) :: ipoids, ivf, idfde, jvGeom, jvMaterc
    integer(kind=8) :: icontp, imass, idepl, idefo, inr
    real(kind=8) :: dff(2, 8), vff(8), b(6, 8), p(3, 6), jac
    real(kind=8) :: dir11(3), densit, pgl(3, 3), distn
    real(kind=8) :: epsm, epsg(9), epsthe, sig, sigg(9), rho, epot
    real(kind=8) :: x(8), y(8), z(8), volume, cdg(3), ppg, xxi, yyi, zzi
    real(kind=8) :: matine(6), vro
    aster_logical :: lexc
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    lexc = (lteatt('MODELI', 'GRC'))

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    if (lexc) then
        call compCoorSystPara(plateCara, zr(jvGeom), pgl)
        nddl = 6
    else
        nddl = 3
    end if
    call compCoorSystGrid(pgl, plateCara, plateOrie)

! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elrefe_info(fami=fami, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)

! - Input fields
    if ((option .eq. 'SIEF_ELGA') .or. (option .eq. 'EPOT_ELEM')) then
        call jevech('PDEPLAR', 'L', idepl)
        call jevech('PMATERC', 'L', jvMaterc)
    else if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEPLAR', 'L', idepl)
        epsg = 0.d0
    else if (option .eq. 'MASS_INER') then
        call jevech('PMATERC', 'L', jvMaterc)
    end if

! - Output fields
    if (option .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', icontp)
        sigg = 0.d0
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PENERDR', 'E', inr)
        epot = 0.d0
    else if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', idefo)
    else if (option .eq. 'MASS_INER') then
        call jevech('PMASSINE', 'E', imass)
        call getDensity(zi(jvMaterc), rho)
    end if

! - LECTURE DES CARACTERISTIQUES DE GRILLE ET CALCUL DE LA DIRECTION D'ARMATURE
    call cargri(plateCara, plateOrie, &
                densit, distn, dir11)

! - COORDONNEES PHYSIQUES DES NOEUDS
    if (option .eq. 'MASS_INER') then
        do i = 1, nno
            x(i) = zr(jvGeom+3*(i-1))
            y(i) = zr(jvGeom+3*i-2)
            z(i) = zr(jvGeom+3*i-1)
        end do
        if (lexc) then
            x(i) = x(i)+plateOrie%gridNorm(1)
            y(i) = y(i)+plateOrie%gridNorm(2)
            z(i) = z(i)+plateOrie%gridNorm(3)
        end if
        cdg = 0.d0
        matine = 0.d0
    end if
!
    volume = 0.d0

    do kpg = 1, npg
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do

! ----- CALCUL DE LA MATRICE "B" : DEPL NODAL --> EPS11 ET DU JACOBIEN
        call nmgrib(nno, zr(jvGeom), dff, dir11, lexc, &
                    plateOrie%gridNorm, b, jac, p)

! --- SIEF_ELGA, EPOT_ELEM : ON CALCULE LA CONTRAINTE AU PG
        if ((option .eq. 'SIEF_ELGA') .or. (option .eq. 'EPOT_ELEM')) then
! --------- Strains
            epsm = 0.d0
            do i = 1, nno
                do j = 1, nddl
                    epsm = epsm+b(j, i)*zr(idepl+(i-1)*nddl+j-1)
                end do
            end do
            call verift(fami, kpg, 1, '+', zi(jvMaterc), &
                        epsth_=epsthe)
            epsm = epsm-epsthe

! --------- Stress
            call rcvalb(fami, kpg, 1, '+', zi(jvMaterc), &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        1, propName, propVale, propCode, 0)
            sig = propVale(1)*epsm
!
            if (option .eq. 'EPOT_ELEM') then
                epot = epot+(sig*epsm*zr(ipoids+kpg-1)*jac*densit)/2
            else
                sigg(kpg) = sig
            end if

        else if (option .eq. 'EPSI_ELGA') then
            do i = 1, nno
                do j = 1, nddl
                    epsg(kpg) = epsg(kpg)+b(j, i)*zr(idepl+(i-1)*nddl+j-1)
                end do
            end do

        else if (option .eq. 'MASS_INER') then
            volume = volume+zr(ipoids+kpg-1)*densit*jac
            ppg = zr(ipoids+kpg-1)*jac*densit
            do i = 1, nno
                cdg(1) = cdg(1)+ppg*vff(i)*x(i)
                cdg(2) = cdg(2)+ppg*vff(i)*y(i)
                cdg(3) = cdg(3)+ppg*vff(i)*z(i)
                xxi = 0.d0
                yyi = 0.d0
                zzi = 0.d0
                do j = 1, nno
                    xxi = xxi+x(i)*vff(i)*vff(j)*x(j)
                    yyi = yyi+y(i)*vff(i)*vff(j)*y(j)
                    zzi = zzi+z(i)*vff(i)*vff(j)*z(j)
                    matine(2) = matine(2)+x(i)*vff(i)*vff(j)*y(j)*ppg
                    matine(4) = matine(4)+x(i)*vff(i)*vff(j)*z(j)*ppg
                    matine(5) = matine(5)+y(i)*vff(i)*vff(j)*z(j)*ppg
                end do
                matine(1) = matine(1)+ppg*(yyi+zzi)
                matine(3) = matine(3)+ppg*(xxi+zzi)
                matine(6) = matine(6)+ppg*(xxi+yyi)
            end do
        end if
    end do
!
    if (option .eq. 'SIEF_ELGA') then
        do kpg = 1, npg
            zr(icontp+kpg-1) = sigg(kpg)
        end do
    else if (option .eq. 'EPOT_ELEM') then
        zr(inr) = epot
    else if (option .eq. 'EPSI_ELGA') then
        do kpg = 1, npg
            zr(idefo+kpg-1) = epsg(kpg)
        end do
    else if (option .eq. 'MASS_INER') then
        vro = rho/volume
        zr(imass) = rho*volume
        zr(imass+1) = cdg(1)/volume
        zr(imass+2) = cdg(2)/volume
        zr(imass+3) = cdg(3)/volume
        zr(imass+4) = matine(1)*rho-vro*(cdg(2)*cdg(2)+cdg(3)*cdg(3))
        zr(imass+5) = matine(3)*rho-vro*(cdg(1)*cdg(1)+cdg(3)*cdg(3))
        zr(imass+6) = matine(6)*rho-vro*(cdg(1)*cdg(1)+cdg(2)*cdg(2))
        zr(imass+7) = matine(2)*rho-vro*(cdg(1)*cdg(2))
        zr(imass+8) = matine(4)*rho-vro*(cdg(1)*cdg(3))
        zr(imass+9) = matine(5)*rho-vro*(cdg(2)*cdg(3))
    end if
!
end subroutine
