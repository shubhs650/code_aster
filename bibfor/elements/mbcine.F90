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
subroutine mbcine(plateOrie, &
                  nno, geom, dff, &
                  b, jac)
!
    use plate_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/subaco.h"
#include "asterfort/subacv.h"
#include "asterfort/sumetr.h"
#include "asterfort/utmess.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: nno
    real(kind=8), intent(in) :: geom(3, nno), dff(2, nno)
    real(kind=8), intent(out) :: b(3, 3, nno), jac
!
! --------------------------------------------------------------------------------------------------
!
!      CALCUL DE LA MATRICE B ET DU JACOBIEN POUR LES MEMBRANES
!
! --------------------------------------------------------------------------------------------------
!
! IN  NNO          NOMBRE DE NOEUDS
! IN  GEOM         COORDONNEES DES NOEUDS
! IN  DFF          DERIVEE DES F. DE FORME
! IN  ALPHA,BETA   ANGLES NAUTIQUES ORIENTANT LE COMPORTEMENT
!                             ORTHOTROPE DE LA MEMBRANE (EN RADIAN)
! OUT B            MATRICE DE PASSAGE DEPL. NODAL --> DEF. MEMBRANAIRES
! OUT JAC          JACOBIEN DE LA TRANSFORMATION
!
! --------------------------------------------------------------------------------------------------
!
!     LES TROIS INDICES DE B CORRESPONDENT RESPECTIVEMENT :
!               - A LA COMPOSANTE DE EPSILON MEMBRANAIRE PARMI
!                     (EPS11, EPS22, SQRT(2)EPS12)
!                     CALCULEES DANS LA BASE LOCALE
!               - A LA COMPOSANTE DU VECTEUR DEPLACEMENT
!               - AU NUMERO DU NOEUD
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, n, gamma
    real(kind=8) :: vdirec(3), vortho(3), alpha, beta
    real(kind=8) :: cova(3, 3), metr(2, 2), cnva(3, 2), a(2, 2)
    real(kind=8) :: projn, denomi
    real(kind=8) :: factor, dicnva(2), orcnva(2)
!
! --------------------------------------------------------------------------------------------------
!
    b = 0.d0
    jac = 0.d0
    ASSERT(plateOrie%lUpdate)
    alpha = plateOrie%alpha
    beta = plateOrie%beta

! - CALCUL DES COORDONNEES COVARIANTES ET CONTRAVARIANTES DE SURFACE
    call subaco(nno, dff, geom, cova)
    call sumetr(cova, metr, jac)
    call subacv(cova, metr, jac, cnva, a)

! - CALCUL ET PROJECTION DU VECTEUR DIRECTION SUR LA SURFACE
    vdirec(1) = cos(beta)*cos(alpha)
    vdirec(2) = cos(beta)*sin(alpha)
    vdirec(3) = -sin(beta)
!
    projn = 0.d0
    do i = 1, 3
        projn = projn+vdirec(i)*cova(i, 3)
    end do
!
    if (abs(1.d0-projn*projn) .le. r8prem()) then
        call utmess('F', 'ELEMENTS_3')
    end if
!
    denomi = sqrt(1.d0-projn*projn)
    do i = 1, 3
        vdirec(i) = (vdirec(i)-projn*cova(i, 3))/denomi
    end do
!
! - CALCUL DU VECTEUR TANGENT ORTHOGONAL AU VECTEUR DIRECTION
!
    vortho(1) = cova(2, 3)*vdirec(3)-cova(3, 3)*vdirec(2)
    vortho(2) = cova(3, 3)*vdirec(1)-cova(1, 3)*vdirec(3)
    vortho(3) = cova(1, 3)*vdirec(2)-cova(2, 3)*vdirec(1)

! - LE TERME DE CISAILLEMENT EST SYMETRISE ET MULTIPLIE PAR SQRT(2)
    factor = 1.d0/sqrt(2.d0)
!
! - ON PRECALCULE CERTAINS PRODUITS SCALAIRES
    dicnva = 0.d0
    orcnva = 0.d0
    do gamma = 1, 2
        do i = 1, 3
            dicnva(gamma) = dicnva(gamma)+vdirec(i)*cnva(i, gamma)
            orcnva(gamma) = orcnva(gamma)+vortho(i)*cnva(i, gamma)
        end do
    end do
!
!
! - ON BOUCLE SUR LES DEGRES DE LIBERTE
    do n = 1, nno
        do i = 1, 3
            do gamma = 1, 2
                b(1, i, n) = b(1, i, n)+dff(gamma, n)*dicnva(gamma)*vdirec(i)
                b(2, i, n) = b(2, i, n)+dff(gamma, n)*orcnva(gamma)*vortho(i)
                b(3, i, n) = b(3, i, n)+factor*dff(gamma, n)*(dicnva(gamma)*vortho(i)+ &
                                                              orcnva(gamma)*vdirec(i))
            end do
        end do
    end do
!
end subroutine
