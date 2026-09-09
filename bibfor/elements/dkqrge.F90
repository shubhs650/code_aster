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
subroutine dkqrge(plateCara, plateOrie, &
                  xyzl, rigiGeom)
!
    use plate_type
    use plateGeom_module, only: isPlateDKT, isPlateDKTG
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/cosiro.h"
#include "asterfort/dkqbnl.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxqloc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gquad4.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/prmama.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8), intent(in) :: xyzl(3, *)
    real(kind=8), intent(out) :: rigiGeom(*)
!
! --------------------------------------------------------------------------------------------------
!
!     matrice de rigiGeomidite geometrique de l'element de plaque dkg
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbsig = 6, nbcon = 8
    integer(kind=8) :: ndim, npg, ipoids, icoopg
    integer(kind=8) :: iLayer, nbLayer
    integer(kind=8) :: iret, ier, idec
    integer(kind=8) :: jtab(7), jsigm, nbsp, npgh
    integer(kind=8) :: i, j, kpg
    real(kind=8) :: poids, caraq4(25)
    real(kind=8) :: hLayer, h, hb, hm, hh, cb, cm, ch
    real(kind=8) :: bnl(2, 12), bnli(12, 2)
    real(kind=8) :: jacob(5), qsi, eta
    real(kind=8) :: flex(12, 12), memb(64), mefl(96)
    real(kind=8) :: flexi(12, 12)
    real(kind=8) :: effint(32), effgt(32)
    real(kind=8) :: ctor
    real(kind=8) :: sixxb, siyyb, sixyb
    real(kind=8) :: sixxm, siyym, sixym
    real(kind=8) :: sixxh, siyyh, sixyh
    real(kind=8) :: nxx, nyy, nxy, normal(2, 2)
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, jpoids=ipoids, jcoopg=icoopg)
!
    flex = 0.d0
    memb = 0.d0
    mefl = 0.d0

! - Get parameters
    h = plateCara%thick
    ctor = plateCara%coefRigiDRZ
    nbLayer = plateCara%nbLayer

! - Contraintes dans les couches
    call tecach('OOO', 'PCONTRR', 'L', iret, nval=7, itab=jtab)
    jsigm = jtab(1)
    npg = jtab(3)
    nbsp = jtab(7)
    npgh = 3
!
    if (isPlateDKT(plateCara)) then
        ASSERT(nbsp .eq. nbLayer*npgh)
        ASSERT(jtab(2) .eq. nbsig*npg)
    end if

! - Passage des contraintes/efforts généralisés dans le repere paramétrique
    if (isPlateDKT(plateCara)) then
        call cosiro(plateCara, plateOrie, &
                    'PCONTRR', 'L', 'UI', 'G', &
                    jsigm)
    else if (isPlateDKTG(plateCara)) then
        ASSERT(plateCara%nbLayer .eq. 1)
        call tecach('OOO', 'PCONTRR', 'L', iret, nval=7, itab=jtab)
        jsigm = jtab(1)
        do i = 1, nbcon*npg
            effgt(i) = zr(jsigm-1+i)
        end do
        call dxefro(npg, plateOrie%t2iu, effgt, effint)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Calcul des grandeurs geometriques
    call gquad4(xyzl, caraq4)

    do kpg = 1, npg
        qsi = zr(icoopg-1+ndim*(kpg-1)+1)
        eta = zr(icoopg-1+ndim*(kpg-1)+2)

! ----- calcul du jacobien sur le quadrangle
        call jquad4(xyzl, qsi, eta, jacob)
        poids = zr(ipoids+kpg-1)*jacob(1)

! ----- calcul des efforts de membrane
        nxx = 0.d0
        nyy = 0.d0
        nxy = 0.d0
        if (isPlateDKT(plateCara)) then
            hb = -h/2
            do iLayer = 1, nbLayer
                idec = ((kpg-1)*nbLayer+(iLayer-1))*npgh*nbsig
                hLayer = h/nbLayer
                hm = hb+hLayer/2.d0
                hh = hm+hLayer/2.d0
!         -- sixxb, siyyb, ... : contraintes au bas de la couche
                sixxb = zr(jsigm-1+idec+1)
                siyyb = zr(jsigm-1+idec+2)
                sixyb = zr(jsigm-1+idec+4)
!         -- sixxm, siyym, ... : contraintes au milieu de la couche
                sixxm = zr(jsigm-1+idec+1+nbsig)
                siyym = zr(jsigm-1+idec+2+nbsig)
                sixym = zr(jsigm-1+idec+4+nbsig)
!         -- sixxh, siyyh, ... : contraintes en haut de la couche
                sixxh = zr(jsigm-1+idec+1+2*nbsig)
                siyyh = zr(jsigm-1+idec+2+2*nbsig)
                sixyh = zr(jsigm-1+idec+4+2*nbsig)
!         -- on integre dans l'epaisseur de chaque couche
!            avec une forrmule de newton-cotes a 3 points
!            les coefficients sont 1/6, 4/6 et 1/6
                cb = hLayer/6
                cm = 4.d0*hLayer/6
                ch = hLayer/6
!         -- nxx, nyy, nxy = somme de sixx, siyy, sixy :
                nxx = nxx+cb*sixxb+cm*sixxm+ch*sixxh
                nyy = nyy+cb*siyyb+cm*siyym+ch*siyyh
                nxy = nxy+cb*sixyb+cm*sixym+ch*sixyh
!         -- mise a jour de hb pour la couche suivante :
                hb = hb+hLayer
            end do
        else if (isPlateDKTG(plateCara)) then
            nxx = effint((kpg-1)*nbcon+1)
            nyy = effint((kpg-1)*nbcon+2)
            nxy = effint((kpg-1)*nbcon+3)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        normal(1, 1) = nxx*poids
        normal(2, 2) = nyy*poids
        normal(1, 2) = nxy*poids
        normal(2, 1) = normal(1, 2)

! ----- Matrice bnl des deformations non-lineaires de membrane
        call dkqbnl(qsi, eta, jacob(2), bnl)

! ----- Bending part
        call prmama(3, bnl, 2, 2, 12, &
                    normal, 2, 2, 2, bnli, &
                    12, 12, 2, ier)
        flexi = matmul(bnli, bnl)
        do i = 1, 12
            do j = 1, 12
                flex(i, j) = flex(i, j)+flexi(i, j)
            end do
        end do
    end do

! - Assemblying sub-matrix
    call dxqloc(flex, memb, mefl, ctor, rigiGeom)
!
end
