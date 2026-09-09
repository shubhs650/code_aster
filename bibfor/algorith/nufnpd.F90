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
subroutine nufnpd(ndim, nno1, nno2, npg, &
                  iw, vff1, vff2, idff1, &
                  vu, vp, &
                  typmod, relaComp, geomi, sig, &
                  ddl, mini, vect, &
                  materPara)
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
#include "asterfort/nmepsi.h"
#include "asterfort/r8inir.h"
#include "asterfort/tanbul.h"
#include "blas/ddot.h"
#include "MeshTypes_type.h"
!
    aster_logical :: mini
    integer(kind=8) :: ndim, nno1, nno2, npg, iw, idff1
    integer(kind=8) :: vu(3, MT_NNOMAX), vp(MT_NNOMAX)
    real(kind=8) :: geomi(ndim, nno1)
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg)
    real(kind=8) :: sig(2*ndim+1, npg), ddl(*), vect(*)
    character(len=8) :: typmod(2)
    character(len=16) :: relaComp
    type(Material_Para), intent(inout) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
!          CALCUL DES FORCES NODALES POUR LES ELEMENTS
!          INCOMPRESSIBLES POUR LES PETITES DEFORMATIONS
!          3D/D_PLAN/AXIS
!          ROUTINE APPELEE PAR TE0591
!
! --------------------------------------------------------------------------------------------------
!
! IN  MINI    : STABILISATION BULLE - MINI ELEMENT
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
! IN  GEOMI   : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  MATE    : MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  DDL     : DEGRES DE LIBERTE A L'INSTANT PRECEDENT
! IN  SIG     : CONTRAINTES A L'INSTANT PRECEDENT
! OUT VECT    : FORCES INTERNES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    aster_logical, parameter :: grand = ASTER_FALSE
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    aster_logical :: axi
    integer(kind=8) :: nddl, kpg
    integer(kind=8) :: sa, na, ia, kk
    real(kind=8) :: deplm(3*MT_NNOMAX)
    real(kind=8) :: presm(MT_NNOMAX), presd(MT_NNOMAX), pm
    real(kind=8) :: dff1(nno1, ndim)
    real(kind=8) :: fm(3, 3)
    real(kind=8) :: r, w
    real(kind=8) :: def(2*ndim, nno1, ndim)
    real(kind=8) :: epsm(6), sigma(6)
    real(kind=8) :: divum
    real(kind=8) :: t1, t2
    real(kind=8) :: alpha
    real(kind=8) :: dsbdep(2*ndim, 2*ndim), kbb(ndim, ndim), kbp(ndim, nno2)
    real(kind=8) :: kce(nno2, nno2), rce(nno2)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    axi = typmod(1) .eq. 'AXIS'
    nddl = nno1*ndim+nno2
    call r8inir(nddl, 0.d0, vect, 1)

! - EXTRACTION DES CHAMPS
    do na = 1, nno1
        do ia = 1, ndim
            deplm(ia+ndim*(na-1)) = ddl(vu(ia, na))
        end do
    end do
!
    do sa = 1, nno2
        presm(sa) = ddl(vp(sa))
        presd(sa) = 0.d0
    end do

! - CALCUL POUR CHAQUE POINT DE GAUSS
    do kpg = 1, npg

! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)
!
! - CALCUL DES ELEMENTS GEOMETRIQUES
        call r8inir(6, 0.d0, epsm, 1)
        call dfdmip(ndim, nno1, axi, geomi, kpg, &
                    iw, vff1(1, kpg), idff1, r, w, &
                    dff1)
        call nmepsi(ndim, nno1, axi, grand, vff1(1, kpg), &
                    r, dff1, deplm, fm, epsm)
!
        divum = epsm(1)+epsm(2)+epsm(3)
!
! - CALCUL DE LA MATRICE B EPS_ij=B_ijkl U_kl
! - DEF (XX,YY,ZZ,2/RAC(2)XY,2/RAC(2)XZ,2/RAC(2)YZ)
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
! - CALCUL DE LA PRESSION
        b_n = to_blas_int(nno2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pm = ddot(b_n, vff2(1, kpg), b_incx, presm, b_incy)
!
! - CALCUL DES CONTRAINTES MECANIQUES A L'EQUILIBRE
        do ia = 1, 3
            sigma(ia) = sig(ia, kpg)
        end do
        do ia = 4, 2*ndim
            sigma(ia) = sig(ia, kpg)*rac2
        end do

! ----- CALCUL DE LA MATRICE D'ELASTICITE BULLE
        call tanbul(materPara, relaComp, &
                    ndim, mini, &
                    alpha, dsbdep)

! ----- CALCUL DE LA MATRICE DE CONDENSATION STATIQUE
        if (mini) then
            call calkbb(nno1, ndim, w, def, dsbdep, &
                        kbb)
            call calkbp(nno2, ndim, w, dff1, kbp)
            call calkce(nno1, ndim, kbp, kbb, presm, &
                        presd, kce, rce)
        else
            call r8inir(nno2, 0.d0, rce, 1)
        end if

! - VECTEUR FINT:U
        do na = 1, nno1
            do ia = 1, ndim
                kk = vu(ia, na)
                b_n = to_blas_int(2*ndim)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                t1 = ddot(b_n, sigma, b_incx, def(1, na, ia), b_incy)
                vect(kk) = vect(kk)+w*t1
            end do
        end do

! - VECTEUR FINT:P
        t2 = (divum-pm*alpha)
        do sa = 1, nno2
            kk = vp(sa)
            t1 = vff2(sa, kpg)*t2
            vect(kk) = vect(kk)+w*t1-rce(sa)
        end do
    end do
!
end subroutine
