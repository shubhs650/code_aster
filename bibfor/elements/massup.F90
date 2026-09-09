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
! aslint: disable=W1306
!
subroutine massup(jvMaterCode, &
                  option, ndim, dlns, nno, nnos, &
                  npg, ipoids, idfde, &
                  geom, vff1, imatuu, icodre, igeom, &
                  ivf)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2j.h"
#include "asterfort/dfdm3j.h"
#include "asterfort/lteatt.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "jeveux.h"
!
    character(len=16) :: option
    integer(kind=8) :: jvMaterCode, ndim, nno, nnos, npg
!
! --------------------------------------------------------------------------------------------------
!
!   CALCUL DE LA MATRICE DE MASSE
!                          POUR ELEMENTS DONT LES NOEUDS SOMMETS
!                          ONT + DE DDL QUE LES DEPLACEMENTS
!
! --------------------------------------------------------------------------------------------------
!
!        DONNEES:   NDIM   -->  DIMENSION DU PROBLEME
!                   DLNS   -->  DEGRES DE LIBERTE AU NOEUD SOMMET
!                   NNO    -->  NOMBRE DE NOEUD
!                   NNOS   -->  NOMBRE DE NOEUD SOMMET
!                   NPG    -->  NOMBRE DE POIDS DE GAUSS
!                   IPOIDS -->  POSITION DES POIDS DE GAUSS DANS ZR
!                   IDFDE  -->
!                   GEOM   -->  COORDONNEES DE L ELEMENT
!                   VFF1   -->  VALEUR DES FONCTIONS DE FORME AUX PG
!                   IMATUU -->  POSITION DE LA MATRICE DE MASSE DANS ZR
!        RESULTATS: ICODRE -->  CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "FPG1"
    character(len=1), parameter :: poum = "+"
    integer(kind=8), parameter :: kpgCst = 1, ksp = 1
    integer(kind=8) :: i, j, k, l, ik, ijkl, dlns, kpg
    integer(kind=8) :: ipoids, idfde, imatuu
    integer(kind=8) :: n1, n2, j2, k2, idiag
    integer(kind=8) :: igeom, ivf, i2, idec
    real(kind=8) :: vff1(nno, npg), geom(ndim, nno), rho(1), r
    real(kind=8) :: a(ndim, ndim, nno, nno), matv(ndim*nno*(ndim*nno+1)/2)
    real(kind=8) :: poids, wgt, trace, alpha
    integer(kind=8) :: icodre(1)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    idec = dlns-ndim

! - Initializations of material parameters on current cell
    call initParaCell(fami, jvMaterCode, materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Initialization for current point
    call initParaPoin(kpgCst, ksp, materPara)

    call r8inir(ndim*ndim*nno*nno, 0.d0, a, 1)
    call r8inir(ndim*nno*(ndim*nno+1)/2, 0.d0, matv, 1)

! - Get parameters
    call rcvalb(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                poum, materPara%jvMaterCode, &
                ' ', materPara%elasKeyword, &
                0, ' ', [0.d0], &
                1, 'RHO', rho, icodre, 1)
!
    if (ndim .eq. 2) then
        do kpg = 1, npg
            k = (kpg-1)*nno
            call dfdm2j(nno, kpg, idfde, geom, poids)
            poids = abs(poids)*zr(ipoids+kpg-1)
            if (lteatt('AXIS', 'OUI')) then
                r = 0.0d0
                do i = 1, nno
                    r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
                end do
                poids = poids*r
            end if
            do i = 1, nno
                do j = 1, i
                    a(1, 1, i, j) = a(1, 1, i, j)+rho(1)*poids*vff1(i, kpg)*vff1(j, kpg)
                    a(2, 2, i, j) = a(1, 1, i, j)
                end do
            end do
        end do

    else if (ndim .eq. 3) then
        do kpg = 1, npg
            call dfdm3j(nno, kpg, idfde, geom, poids)
            poids = abs(poids)*zr(ipoids+kpg-1)
            do i = 1, nno
                do j = 1, i
                    a(1, 1, i, j) = a(1, 1, i, j)+rho(1)*poids*vff1(i, kpg)*vff1(j, kpg)
                    a(2, 2, i, j) = a(1, 1, i, j)
                    a(3, 3, i, j) = a(1, 1, i, j)
                end do
            end do
        end do
    else
        ASSERT(ASTER_FALSE)
    end if

! - PASSAGE DU STOCKAGE RECTANGULAIRE (A) AU STOCKAGE TRIANGULAIRE (ZR)
    do k = 1, ndim
        do l = 1, ndim
            do i = 1, nno
                ik = ((ndim*i+k-ndim-1)*(ndim*i+k-ndim))/2
                do j = 1, i
                    ijkl = ik+ndim*(j-1)+l
                    matv(ijkl) = a(k, l, i, j)
                end do
            end do
        end do
    end do
!
    if (option .eq. 'MASS_MECA') then
        do k = 1, nno
            do n1 = 1, ndim
                i = ndim*k+n1-ndim
                if (k .le. nnos) then
                    i2 = i+idec*(k-1)
                else
                    i2 = i+idec*nnos
                end if
                do l = 1, nno
                    do n2 = 1, ndim
                        j = ndim*l+n2-ndim
                        if (j .gt. i) goto 405
                        if (l .le. nnos) then
                            j2 = j+idec*(l-1)
                        else
                            j2 = j+idec*nnos
                        end if
                        zr(imatuu+i2*(i2-1)/2+j2-1) = matv(i*(i-1)/2+j)
                    end do
                end do
405             continue
            end do
        end do
    elseif (option .eq. 'MASS_MECA_DIAG' .or. &
            option .eq. 'MASS_MECA_EXPLI') then
! - CALCUL DE LA MASSE DE L'ELEMENT
        wgt = a(1, 1, 1, 1)
        do i = 2, nno
            do j = 1, i-1
                wgt = wgt+2*a(1, 1, i, j)
            end do
            wgt = wgt+a(1, 1, i, i)
        end do
!
! - CALCUL DE LA TRACE EN TRANSLATION SUIVANT X
        trace = 0.d0
        do i = 1, nno
            trace = trace+a(1, 1, i, i)
        end do
!
! - CALCUL DU FACTEUR DE DIAGONALISATION
        alpha = wgt/trace
!
! - PASSAGE DU STOCKAGE RECTANGULAIRE (A) AU STOCKAGE TRIANGULAIRE (ZR)
        k = 0
        do j = 1, nno
            do i = 1, 3
                k = k+1
                if (idec .eq. 0) then
                    idiag = k*(k+1)/2
                else
                    if (j .le. nnos) then
                        k2 = k+idec*(j-1)
                    else
                        k2 = k+idec*nnos
                    end if
                    idiag = k2*(k2+1)/2
                end if
                zr(imatuu+idiag-1) = a(i, i, j, j)*alpha
            end do
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
