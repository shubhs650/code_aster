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
subroutine te0392(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/caatdb.h"
#include "asterfort/cast3d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/dmatmc.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elraga.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/invjac.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_SI
! Option: RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: idfde2, kpg, jvMaterc, imatuu, ipoid2
    integer(kind=8) :: nbsig, nno, npg
    real(kind=8) :: jacgau, time
    integer(kind=8) :: jvGeom, ipoids, ivf, idfde
    aster_logical :: calbn
    integer(kind=8) :: i, ino, j, k, proj, nbpg2
    integer(kind=8) :: ndim, nnos, kp
    real(kind=8) :: d(6, 6), s
    real(kind=8) :: poipg2(8), b(6, 81), b0(6, 3, 8)
    real(kind=8) :: jac, invja(3, 3), bi(3, 8), hx(3, 4)
    real(kind=8) :: gam(4, 8), coopg2(24), h(8, 4), dh(4, 24)
    real(kind=8) :: bn(6, 3, 8)
    real(kind=8) :: dfdx(8), dfdy(8), dfdz(8)
    real(kind=8) :: nu, nub, nu12
    data h/1.d0, 1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 1.d0, 1.d0,&
     &        1.d0, -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0,&
     &        1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0,&
     &       -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0, 1.d0, -1.d0/
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!

! - Finite element informations
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)

! - Initializations
    b = 0.d0
    nbsig = nbsigm()

! - Get current time
    time = r8vide()

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Set output
    call jevech('PMATUUR', 'E', imatuu)
    do i = 1, 300
        zr(imatuu-1+i) = 0.0d0
    end do

! - Compute [Bi] (mean value for derivate of shape functions)
    do kpg = 1, npg
        call dfdm3d(nno, kpg, ipoids, idfde, zr(jvGeom), &
                    jac, dfdx, dfdy, dfdz)
        do ino = 1, nno
            bi(1, ino) = dfdx(ino)
            bi(2, ino) = dfdy(ino)
            bi(3, ino) = dfdz(ino)
        end do
    end do
!
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Compute matrix [B]: displacement -> strain (first order)
        call dfdm3d(nno, kpg, ipoids, idfde, zr(jvGeom), &
                    jacgau, dfdx, dfdy, dfdz)

! ----- Modify matrix [B] for underintegrated elements
        do i = 1, 8
            j = 3*(i-1)+1
            b(1, j) = bi(1, i)
            b(2, j+1) = bi(2, i)
            b(3, j+2) = bi(3, i)
            b(4, j) = bi(2, i)
            b(4, j+1) = bi(1, i)
            b(5, j) = bi(3, i)
            b(5, j+2) = bi(1, i)
            b(6, j+1) = bi(3, i)
            b(6, j+2) = bi(2, i)
        end do
        do i = 1, nno
            do j = 1, 3
                do k = 1, 6
                    b0(k, j, i) = b(k, (i-1)*3+j)
                end do
            end do
        end do

! ----- Compute Hooke matrix [D]
        call dmatmc(materPara, '+', time, &
                    nbsig, d)

! ----- Compute "center" rigidity matrix [KC]
        call caatdb(nno, b0, d, b0, jacgau, &
                    zr(imatuu))

    end do

! - Gamma ratio
    do i = 1, 4
        do k = 1, 3
            hx(k, i) = 0.d0
            do j = 1, nno
                hx(k, i) = hx(k, i)+h(j, i)*zr(jvGeom-1+3*(j-1)+k)
            end do
        end do
    end do
!
    do i = 1, 4
        do j = 1, nno
            s = 0.d0
            do k = 1, 3
                s = s+hx(k, i)*bi(k, j)
            end do
            gam(i, j) = 0.125d0*(h(j, i)-s)
        end do
    end do

! - Poisson ration for ASQBI
    kpg = 1
    call get_elas_para(fami, zi(jvMaterc), '+', kpg, ksp, &
                       materPara%elasID, materPara%elasKeyword, &
                       nu_=nu, nu12_=nu12)
    if (materPara%elasID .eq. ELAS_ISOT) then
        nub = nu/(1.d0-nu)
    else
        nub = nu12/(1.d0-nu12)
    end if

! - Projection type
!           0 AUCUNE
!           1 ADS
!           2 ASBQI
    proj = 2
    calbn = ASTER_FALSE

! - Finite element informations for underintegrated element
    call elraga('HE8', 'FPG8    ', ndim, nbpg2, coopg2, poipg2)
    call elrefe_info(elrefe='HE8', fami='MASS', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=nbpg2, jpoids=ipoid2, jdfde=idfde2)

! - Compute corrected stabilization matrix [K_STAB]
    do kpg = 1, nbpg2
        kp = 3*(kpg-1)
        call invjac(nno, kpg, ipoid2, idfde2, zr(jvGeom), &
                    invja, jac)
        do i = 1, 3
            dh(1, kp+i) = coopg2(3*kpg-1)*invja(i, 3)+coopg2(3*kpg)*invja(i, 2)
        end do
        do i = 1, 3
            dh(2, kp+i) = coopg2(3*kpg-2)*invja(i, 3)+coopg2(3*kpg)*invja(i, 1)
        end do
        do i = 1, 3
            dh(3, kp+i) = coopg2(3*kpg-2)*invja(i, 2)+coopg2(3*kpg-1)*invja(i, 1)
        end do
        do i = 1, 3
            dh(4, kp+i) = coopg2(3*kpg-2)*coopg2(3*kpg-1)*invja(i, 3)+coopg2(3*kpg-1)*coop&
                         &g2(3*kpg)*invja(i, 1)+coopg2(3*kpg-2)*coopg2(3*kpg)*invja(i, 2)
        end do
        call cast3d(proj, gam, dh, b0, nno, &
                    kpg, nub, nu, d, calbn, &
                    bn, jac, zr(imatuu))
    end do
!
end subroutine
