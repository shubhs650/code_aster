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
subroutine te0436(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystMemb
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getDensity.h"
#include "asterfort/jevech.h"
#include "asterfort/mbcine.h"
#include "asterfort/mbrigi.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: MEMBRANE
!
! Options: EPOT_ELEM, SIEF_ELGA, EPSI_ELGA, EFGE_ELGA, MASS_INER
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: nddl = 3, ncomp = 3
    integer(kind=8) :: nno, npg
    integer(kind=8) :: i, j, n, c, cc, kpg, jvSief, jvEfge, k
    integer(kind=8) :: ipoids, ivf, idfde, iret
    integer(kind=8) :: jvGeom, jvMaterc, jvDisp, jvEner, jvEpsi, jvMassIner, jcCompor
    real(kind=8) :: dff(2, 9), vff(9), b(3, 3, 9), jac
    real(kind=8) :: epot
    real(kind=8) :: epsm(3), epsg(3, 9), epsthe, sigmMemb(3), sigg(3, 9), matrRigi(3, 3)
    real(kind=8) :: rho, rhog
    real(kind=8) :: x(9), y(9), z(9), surfac, cdg(3), ppg, xxi, yyi, zzi
    real(kind=8) :: massIner(6)
    real(kind=8) :: vro
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Compute global<=>local transformation
    call compCoorSystMemb(plateOrie)

! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elrefe_info(fami=fami, nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde)

    if (option .eq. 'EFGE_ELGA') then
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PEFGER', 'E', jvEfge)
        do k = 1, 3*npg
            zr(jvEfge-1+k) = zr(jvSief-1+k)
        end do
        goto 999
    end if

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Input fields
    if ((option .eq. 'SIEF_ELGA') .or. (option .eq. 'EPOT_ELEM')) then
        call jevech('PDEPLAR', 'L', jvDisp)
        call jevech('PMATERC', 'L', jvMaterc)
    else if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEPLAR', 'L', jvDisp)
        epsg = 0.d0
    else if (option .eq. 'MASS_INER') then
        call jevech('PMATERC', 'L', jvMaterc)
    end if

! - Output field
    if (option .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', jvSief)
        sigg = 0.d0
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PENERDR', 'E', jvEner)
        epot = 0.d0
    else if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', jvEpsi)
    else if (option .eq. 'MASS_INER') then
        call jevech('PMASSINE', 'E', jvMassIner)
    end if

! - ON INTERDIT CERTAINES OPTIONS POUR LES GRANDES DEFORMATIONS
    call tecach('NNO', 'PCOMPOR', 'L', iret, iad=jcCompor)
    if (((option .eq. 'EPSI_ELGA') .or. (option .eq. 'EPOT_ELEM')) .and. &
        (iret .eq. 0) .and. (zk16(jcCompor+2) (1:9) .eq. 'GROT_GDEP')) then
        call utmess('F', 'MEMBRANE_8', sk=option)
    end if

! - COORDONNEES PHYSIQUES DES NOEUDS
    if (option .eq. 'MASS_INER') then
        do i = 1, nno
            x(i) = zr(jvGeom+3*(i-1))
            y(i) = zr(jvGeom+3*i-2)
            z(i) = zr(jvGeom+3*i-1)
        end do
        cdg = 0.d0
        massIner = 0.d0
        surfac = 0.d0
        rhog = 0.d0
    end if

    do kpg = 1, npg
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do

! ----- CALCUL DE LA MATRICE "B" :
        call mbcine(plateOrie, &
                    nno, zr(jvGeom), dff, &
                    b, jac)
!
! --- SIEF_ELGA, EPOT_ELEM : ON CALCULE LA CONTRAINTE AU PG
!
        if ((option .eq. 'SIEF_ELGA') .or. (option .eq. 'EPOT_ELEM')) then
! ------    CALCUL DE LA DEFORMATION MEMBRANAIRE DANS LE REPERE LOCAL
            epsm = 0.d0
            do n = 1, nno
                do i = 1, nddl
                    do c = 1, ncomp
                        epsm(c) = epsm(c)+b(c, i, n)*zr(jvDisp+(n-1)*nddl+i-1)
                    end do
                end do
            end do

! ------    RETRAIT DE LA DEFORMATION THERMIQUE
            call verift(fami, kpg, 1, '+', zi(jvMaterc), &
                        epsth_=epsthe)
            epsm(1) = epsm(1)-epsthe
            epsm(2) = epsm(2)-epsthe

! ------    CALCUL DE LA CONTRAINTE AU PG
            call mbrigi(fami, kpg, jvMaterc, matrRigi)
            sigmMemb = 0.d0
            do c = 1, ncomp
                do cc = 1, ncomp
                    sigmMemb(c) = sigmMemb(c)+epsm(cc)*matrRigi(cc, c)
                end do
            end do
            if (option .eq. 'EPOT_ELEM') then
                do c = 1, ncomp
                    epot = epot+(sigmMemb(c)*epsm(c)*zr(ipoids+kpg-1)*jac)/2
                end do
            else
                do c = 1, ncomp
                    sigg(c, kpg) = sigmMemb(c)
                end do
            end if

        else if (option .eq. 'EPSI_ELGA') then
            do n = 1, nno
                do i = 1, nddl
                    do c = 1, ncomp
                        epsg(c, kpg) = epsg(c, kpg)+ &
                                       b(c, i, n)*zr(jvDisp+(n-1)*nddl+i-1)
                    end do
                end do
            end do

        else if (option .eq. 'MASS_INER') then
            call getDensity(zi(jvMaterc), rho, 'ELAS_MEMBRANE')
            surfac = surfac+zr(ipoids+kpg-1)*jac
            rhog = rhog+rho*zr(ipoids+kpg-1)*jac
            ppg = zr(ipoids+kpg-1)*jac
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
                    massIner(2) = massIner(2)+x(i)*vff(i)*vff(j)*y(j)*ppg
                    massIner(4) = massIner(4)+x(i)*vff(i)*vff(j)*z(j)*ppg
                    massIner(5) = massIner(5)+y(i)*vff(i)*vff(j)*z(j)*ppg
                end do
                massIner(1) = massIner(1)+ppg*(yyi+zzi)
                massIner(3) = massIner(3)+ppg*(xxi+zzi)
                massIner(6) = massIner(6)+ppg*(xxi+yyi)
            end do
        end if
    end do

    if (option .eq. 'SIEF_ELGA') then
        do kpg = 1, npg
            do c = 1, ncomp
                zr(jvSief+(kpg-1)*ncomp+c-1) = sigg(c, kpg)
            end do
        end do

    else if (option .eq. 'EPOT_ELEM') then
        zr(jvEner) = epot

    else if (option .eq. 'EPSI_ELGA') then
        do kpg = 1, npg
            do c = 1, ncomp
                zr(jvEpsi+(kpg-1)*ncomp+c-1) = epsg(c, kpg)
            end do
        end do

    else if (option .eq. 'MASS_INER') then
        vro = rhog/surfac
        zr(jvMassIner) = rhog*surfac
        zr(jvMassIner+1) = cdg(1)/surfac
        zr(jvMassIner+2) = cdg(2)/surfac
        zr(jvMassIner+3) = cdg(3)/surfac
        zr(jvMassIner+4) = massIner(1)*rhog-vro*(cdg(2)*cdg(2)+cdg(3)*cdg(3))
        zr(jvMassIner+5) = massIner(3)*rhog-vro*(cdg(1)*cdg(1)+cdg(3)*cdg(3))
        zr(jvMassIner+6) = massIner(6)*rhog-vro*(cdg(1)*cdg(1)+cdg(2)*cdg(2))
        zr(jvMassIner+7) = massIner(2)*rhog-vro*(cdg(1)*cdg(2))
        zr(jvMassIner+8) = massIner(4)*rhog-vro*(cdg(1)*cdg(3))
        zr(jvMassIner+9) = massIner(5)*rhog-vro*(cdg(2)*cdg(3))
    end if

999 continue

end subroutine
