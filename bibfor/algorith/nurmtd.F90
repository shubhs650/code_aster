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
subroutine nurmtd(ndim, nno1, nno2, npg, iw, &
                  vff1, vff2, idff1, vu, &
                  vp, typmod, jvGeom, jvMaterCode, mini, &
                  matr)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calkbb.h"
#include "asterfort/calkbp.h"
#include "asterfort/calkce.h"
#include "asterfort/dfdmip.h"
#include "asterfort/r8inir.h"
#include "asterfort/tanbul.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    aster_logical :: mini
    integer(kind=8) :: ndim, nno1, nno2, npg, iw, idff1
    integer(kind=8) :: vu(3, MT_NNOMAX), vp(MT_NNOMAX)
    integer(kind=8) :: jvGeom, jvMaterCode
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg)
    character(len=8) :: typmod(*)
    real(kind=8) :: matr(*)
!
! --------------------------------------------------------------------------------------------------
!
!          CALCUL DE LA RIGIDITE MECANIQUE POUR LES ELEMENTS
!          INCOMPRESSIBLES POUR LES PETITES DEFORMATIONS
!          3D/D_PLAN/AXIS
!          ROUTINE APPELEE PAR TE0597
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO1    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AUX DEPLACEMENTS
! IN  NNO2    : NOMBRE DE NOEUDS DE L'ELEMENT LIES A LA PRESSION
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS
! IN  VFF1    : VALEUR  DES FONCTIONS DE FORME LIES AUX DEPLACEMENTS
! IN  VFF2    : VALEUR  DES FONCTIONS DE FORME LIES A LA PRESSION
! IN  IDFF1   : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  VU      : TABLEAU DES INDICES DES DDL DE DEPLACEMENTS
! IN  VP      : TABLEAU DES INDICES DES DDL DE PRESSION
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IGEOM   : POINTEUR SUR LES COORDONEES DES NOEUDS
! IN  MATE    : MATERIAU CODE
! IN  MINI    : STABILISATION BULLE - MINI ELEMENT
! OUT MATR    : MATRICE DE RIGIDITE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    character(len=8), parameter :: fami = "RIGI"
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    character(len=16), parameter :: relaComp = 'ELAS'
    aster_logical, parameter :: miniForBubble = ASTER_TRUE
    aster_logical :: axi
    integer(kind=8) :: kpg
    integer(kind=8) :: ia, na, sa, ib, nb, sb, ja, jb
    integer(kind=8) :: os, kk
    integer(kind=8) :: vuiana, vpsa
    integer(kind=8) :: idim
    real(kind=8) :: r, w, dff1(nno1, ndim)
    real(kind=8) :: dsidep(2*ndim, 2*ndim)
    real(kind=8) :: def(2*ndim, nno1, ndim), deftr(nno1, ndim)
    real(kind=8) :: ddev(2*ndim, 2*ndim), devd(2*ndim, 2*ndim)
    real(kind=8) :: dddev(2*ndim, 2*ndim)
    real(kind=8) :: coorBary(3)
    real(kind=8) :: t1
    real(kind=8) :: idev(6, 6), idev2(4, 4)
    real(kind=8) :: alpha
    real(kind=8) :: presm(nno2), presd(nno2)
    real(kind=8) :: kbb(ndim, ndim), kbp(ndim, nno2)
    real(kind=8) :: kce(nno2, nno2), rce(nno2)
    real(kind=8) :: fm(3, 3)
!
    data fm/1.d0, 0.d0, 0.d0,&
     &                  0.d0, 1.d0, 0.d0,&
     &                  0.d0, 0.d0, 1.d0/
    data idev2/2.d0, -1.d0, -1.d0, 0.d0,&
     &                 -1.d0, 2.d0, -1.d0, 0.d0,&
     &                 -1.d0, -1.d0, 2.d0, 0.d0,&
     &                  0.d0, 0.d0, 0.d0, 3.d0/
    data idev/2.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0,&
     &                 -1.d0, 2.d0, -1.d0, 0.d0, 0.d0, 0.d0,&
     &                 -1.d0, -1.d0, 2.d0, 0.d0, 0.d0, 0.d0,&
     &                  0.d0, 0.d0, 0.d0, 3.d0, 0.d0, 0.d0,&
     &                  0.d0, 0.d0, 0.d0, 0.d0, 3.d0, 0.d0,&
     &                  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 3.d0/
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    axi = typmod(1) .eq. 'AXIS'
    call r8inir(nno2, 0.d0, presm, 1)
    call r8inir(nno2, 0.d0, presd, 1)

! - RECUPERATION  DES DONNEEES RELATIVES AU REPERE D'ORTHOTROPIE
! - COORDONNEES DU BARYCENTRE ( POUR LE REPRE CYLINDRIQUE )
    coorBary = 0.d0
    do ia = 1, nno1
        do idim = 1, ndim
            coorBary(idim) = coorBary(idim)+zr(jvGeom+idim+ndim*(ia-1)-1)/nno1
        end do
    end do

! - Initializations of material parameters on current cell
    call initParaCell(fami, jvMaterCode, materPara)

! - Set local coordinate system from user
    call getUserLCSWithBaryCenter(ndim, coorBary, materPara%lcsPara)

! - CALCUL POUR CHAQUE POINT DE GAUSS
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! - CALCUL DES ELEMENTS GEOMETRIQUES
! - CALCUL DE DFDI,F,EPS,R(EN AXI) ET POIDS
        call dfdmip(ndim, nno1, axi, zr(jvGeom), kpg, &
                    iw, vff1(1, kpg), idff1, r, w, &
                    dff1)
!
! - CALCUL DE LA MATRICE B EPS_ij=B_ijkl U_kl
        if (ndim .eq. 2) then
            do na = 1, nno1
                do ia = 1, ndim
                    def(1, na, ia) = fm(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = fm(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = 0.d0
                    def(4, na, ia) = (fm(ia, 1)*dff1(na, 2)+fm(ia, 2)*dff1(na, 1))/rac2
                end do
            end do
!
! - TERME DE CORRECTION (3,3) AXI QUI PORTE EN FAIT SUR LE DDL 1
            if (axi) then
                do na = 1, nno1
                    def(3, na, 1) = fm(3, 3)*vff1(na, kpg)/r
                end do
            end if
        else
            do na = 1, nno1
                do ia = 1, ndim
                    def(1, na, ia) = fm(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = fm(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = fm(ia, 3)*dff1(na, 3)
                    def(4, na, ia) = (fm(ia, 1)*dff1(na, 2)+fm(ia, 2)*dff1(na, 1))/rac2
                    def(5, na, ia) = (fm(ia, 1)*dff1(na, 3)+fm(ia, 3)*dff1(na, 1))/rac2
                    def(6, na, ia) = (fm(ia, 2)*dff1(na, 3)+fm(ia, 3)*dff1(na, 2))/rac2
                end do
            end do
        end if
!
! - CALCUL DE TRACE(B)
        do na = 1, nno1
            do ia = 1, ndim
                deftr(na, ia) = def(1, na, ia)+def(2, na, ia)+def(3, na, ia)
            end do
        end do

! ----- Compute "bubble" matrix
        call tanbul(materPara, relaComp, &
                    ndim, miniForBubble, &
                    alpha, dsidep)

        dsidep(4, 4) = dsidep(4, 4)/2.d0
        if (ndim .eq. 3) then
            dsidep(5, 5) = dsidep(5, 5)/2.d0
            dsidep(6, 6) = dsidep(6, 6)/2.d0
        end if
!
! - CALCUL DE LA MATRICE DE CONDENSATION STATIQUE
        if (mini) then
            call calkbb(nno1, ndim, w, def, dsidep, &
                        kbb)
            call calkbp(nno2, ndim, w, dff1, kbp)
            call calkce(nno1, ndim, kbp, kbb, presm, &
                        presd, kce, rce)
        else
            call r8inir(nno2*nno2, 0.d0, kce, 1)
        end if
        devd(:, :) = 0.0d0
        ddev(:, :) = 0.0d0
        dddev(:, :) = 0.0d0
!
        if (ndim .eq. 3) then
            devd(1:6, 1:6) = matmul(idev/3.d0, dsidep(1:6, 1:6))
            ddev(1:6, 1:6) = matmul(dsidep(1:6, 1:6), idev/3.d0)
            dddev(1:6, 1:6) = matmul(devd(1:6, 1:6), idev/3.d0)
        else
            devd(1:4, 1:4) = matmul(idev2/3.d0, dsidep(1:4, 1:4))
            ddev(1:4, 1:4) = matmul(dsidep(1:4, 1:4), idev2/3.d0)
            dddev(1:4, 1:4) = matmul(devd(1:4, 1:4), idev2/3.d0)
        end if
!
! - CALCUL DE LA MATRICE DE RIGIDITE
! - TERME K:UX
        do na = 1, nno1
            do ia = 1, ndim
                vuiana = vu(ia, na)
                os = (vuiana-1)*vuiana/2
!
! - TERME K:UU      KUU(NDIM,NNO1,NDIM,NNO1)
                do nb = 1, nno1
                    do ib = 1, ndim
                        if (vu(ib, nb) .le. vuiana) then
                            kk = os+vu(ib, nb)
                            t1 = 0.d0
                            do ja = 1, 2*ndim
                                do jb = 1, 2*ndim
                                    t1 = t1+def(ja, na, ia)*dddev(ja, jb)*def(jb, nb, ib)
                                end do
                            end do
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
                end do
!
! - TERME K:UP      KUP(NDIM,NNO1,NNO2)
                do sb = 1, nno2
                    if (vp(sb) .lt. vuiana) then
                        kk = os+vp(sb)
                        t1 = deftr(na, ia)*vff2(sb, kpg)
                        matr(kk) = matr(kk)+w*t1
                    end if
                end do
            end do
        end do
!
! - TERME K:PX
        do sa = 1, nno2
            vpsa = vp(sa)
            os = (vpsa-1)*vpsa/2
!
! - TERME K:PU      KPU(NDIM,NNO2,NNO1)
            do nb = 1, nno1
                do ib = 1, ndim
                    if (vu(ib, nb) .lt. vpsa) then
                        kk = os+vu(ib, nb)
                        t1 = vff2(sa, kpg)*deftr(nb, ib)
                        matr(kk) = matr(kk)+w*t1
                    end if
                end do
            end do
!
! - TERME K:PP      KPP(NNO2,NNO2)
            do sb = 1, nno2
                if (vp(sb) .le. vpsa) then
                    kk = os+vp(sb)
                    t1 = -vff2(sa, kpg)*vff2(sb, kpg)*alpha
                    matr(kk) = matr(kk)+w*t1-kce(sa, sb)
                end if
            end do
        end do
    end do
end subroutine
