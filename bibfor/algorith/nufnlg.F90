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
subroutine nufnlg(ndim, nno1, nno2, npg, &
                  iw, vff1, vff2, idff1, &
                  vu, vp, &
                  typmod, relaComp, geomi, sig, &
                  ddl, vect, &
                  materPara)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmepsi.h"
#include "asterfort/nmmalu.h"
#include "asterfort/r8inir.h"
#include "asterfort/tanbul.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
#include "MeshTypes_type.h"
!
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
!          INCOMPRESSIBLES POUR LES GRANDES DEFORMATIONS
!          3D/D_PLAN/AXIS
!          ROUTINE APPELEE PAR TE0596
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
    aster_logical, parameter :: grand = ASTER_TRUE, mini = ASTER_FALSE
    aster_logical :: axi
    integer(kind=8) :: vij(3, 3), lij(3, 3)
    integer(kind=8) :: nddl, ndu, kpg
    integer(kind=8) :: kl, sa, na, ia, ja, kk
    real(kind=8) :: geomm(3*MT_NNOMAX), jm, wm, epsm(6)
    real(kind=8) :: deplm(3*MT_NNOMAX), presm(MT_NNOMAX), pm
    real(kind=8) :: dff1(nno1, 4)
    real(kind=8) :: fm(3, 3)
    real(kind=8) :: r, w
    real(kind=8) :: tau(6), taudv(6), tauhy
    real(kind=8) :: t1, t2
    real(kind=8) :: kr(6), id(3, 3)
    real(kind=8) :: alpha
    real(kind=8) :: dsbdep(2*ndim, 2*ndim)
    blas_int :: b_incx, b_incy, b_n
    data vij/1, 4, 5,&
     &                  4, 2, 6,&
     &                  5, 6, 3/
    data kr/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
    data id/1.d0, 0.d0, 0.d0,&
     &                  0.d0, 1.d0, 0.d0,&
     &                  0.d0, 0.d0, 1.d0/
!
! --------------------------------------------------------------------------------------------------
!
    axi = typmod(1) .eq. 'AXIS'
    nddl = nno1*ndim+nno2
    ndu = ndim
    if (axi) ndu = 3
!
    call r8inir(nddl, 0.d0, vect, 1)
    call r8inir(6, 0.d0, tau, 1)

! - REACTUALISATION DE LA GEOMETRIE ET EXTRACTION DES CHAMPS
    do na = 1, nno1
        do ia = 1, ndim
            geomm(ia+ndim*(na-1)) = geomi(ia, na)+ddl(vu(ia, na))
            deplm(ia+ndim*(na-1)) = ddl(vu(ia, na))
        end do
    end do
!
    do sa = 1, nno2
        presm(sa) = ddl(vp(sa))
    end do

! - CALCUL POUR CHAQUE POINT DE GAUSS
    do kpg = 1, npg

! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)
!
! - CALCUL DES ELEMENTS GEOMETRIQUES
        call dfdmip(ndim, nno1, axi, geomi, kpg, &
                    iw, vff1(1, kpg), idff1, r, w, &
                    dff1)
        call nmepsi(ndim, nno1, axi, grand, vff1(1, kpg), &
                    r, dff1, deplm, fm, epsm)
        call dfdmip(ndim, nno1, axi, geomm, kpg, &
                    iw, vff1(1, kpg), idff1, r, wm, &
                    dff1)
        call nmmalu(nno1, axi, r, vff1(1, kpg), dff1, &
                    lij)
!
        jm = fm(1, 1)*(fm(2, 2)*fm(3, 3)-fm(2, 3)*fm(3, 2))- &
             fm(2, 1)*(fm(1, 2)*fm(3, 3)-fm(1, 3)*fm(3, 2))+ &
             fm(3, 1)*(fm(1, 2)*fm(2, 3)-fm(1, 3)*fm(2, 2))
!
! - CALCUL DE LA PRESSION
        b_n = to_blas_int(nno2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pm = ddot(b_n, vff2(1, kpg), b_incx, presm, b_incy)
!
! - CONTRAINTE DE KIRCHHOFF
        b_n = to_blas_int(2*ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sig(1, kpg), b_incx, tau, b_incy)
        b_n = to_blas_int(2*ndim)
        b_incx = to_blas_int(1)
        call dscal(b_n, jm, tau, b_incx)
        tauhy = (tau(1)+tau(2)+tau(3))/3.d0
        do kl = 1, 6
            taudv(kl) = tau(kl)-tauhy*kr(kl)
        end do

! ----- CALCUL DE LA MATRICE D'ELASTICITE BULLE
        call tanbul(materPara, relaComp, &
                    ndim, mini, &
                    alpha, dsbdep)

! - VECTEUR FINT:U
        do na = 1, nno1
            do ia = 1, ndu
                kk = vu(ia, na)
                t1 = 0.d0
                do ja = 1, ndu
                    t2 = taudv(vij(ia, ja))+pm*id(ia, ja)
                    t1 = t1+t2*dff1(na, lij(ia, ja))
                end do
                vect(kk) = vect(kk)+w*t1
            end do
        end do

! - VECTEUR FINT:P
        t2 = log(jm)-pm*alpha
        do sa = 1, nno2
            kk = vp(sa)
            t1 = vff2(sa, kpg)*t2
            vect(kk) = vect(kk)+w*t1
        end do
    end do
!
end subroutine
