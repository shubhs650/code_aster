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
subroutine xcinem(axi, igeom, nnop, nnos, idepl, &
                  ndim, he, nfiss, nfh, &
                  nfe, ddls, ddlm, fk, dkdgl, &
                  ff, dfdi, f, eps, grad, &
                  heavn)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/indent.h"
#include "asterfort/xcalc_code.h"
#include "MeshTypes_type.h"
#include "asterfort/xcalc_heav.h"
!
    aster_logical, intent(in) :: axi
    integer(kind=8), intent(in) :: igeom
    integer(kind=8), intent(in) :: nnop
    integer(kind=8), intent(in) :: nnos
    integer(kind=8), intent(in) :: idepl
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: nfiss
    real(kind=8), intent(in) :: he(nfiss)
    integer(kind=8), intent(in) :: nfh
    integer(kind=8), intent(in) :: nfe
    integer(kind=8), intent(in) :: ddls
    integer(kind=8), intent(in) :: ddlm
    integer(kind=8), intent(in) :: heavn(nnop, 5)
    real(kind=8), intent(in) :: fk(MT_NNOMAX3D, 3, 3)
    real(kind=8), intent(in) :: dkdgl(MT_NNOMAX3D, 3, 3, 3)
    real(kind=8), intent(in) :: ff(nnop)
    real(kind=8), intent(in) :: dfdi(nnop, ndim)
    real(kind=8), intent(out) :: f(3, 3)
    real(kind=8), intent(out) :: eps(6)
    real(kind=8), intent(out) :: grad(ndim, ndim)
!
! ----------------------------------------------------------------------
!
! X-FEM : CALCUL DES ELEMENTS CINEMATIQUES F, EPS ET GRAD CONNAISSANT
!         FF ET DFDI (VALEURS DES FF CLASSIQUES ET DE LEUR DERIVEES DANS
!         LE REPERE DE REFERENCE PREALABLEMENT CALCULEES AVEC REEREF)
!
! ----------------------------------------------------------------------
!
!
! IN   AXI   : INDIQUER POUR MODEL AXIS
! IN  NNOP   : NOMBRE DE NOEUDS DE L'ELT DE RÉF PARENT
!   L'ORDRE DES DDLS DOIT ETRE 'DC' 'H1' 'E1' 'E2' 'E3' 'E4' 'LAGC'
! IN  DEPL   : DEPLACEMENT RÉEL À PARTIR DE LA CONF DE REF
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  HE     : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ÉLT
! IN  R      : RADIUS POUR CALCULER EPSILON_33 POUR AXI
! IN UR      : DEPLACEMNET RADIAL POUR CALCULER EPSILON_33 POUR AXI
! IN  NFH    : NOMBRE DE FONCTIONS HEAVYSIDE (PAR NOEUD)
! IN  NFE    : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  DDLT   : NOMBRE DE DDLS TOTAL PAR NOEUD
! IN  DKDGL  : DÉRIVÉES DES FONCTIONS D'ENRICHISSEMENT
! IN  FF     : FONCTIONS DE FORMES EN XE
! IN  DFDI   : DÉRIVÉES DES FONCTIONS DE FORMES EN XE
! OUT F      : GRADIENT DE LA TRANSFORMATION
! OUT EPS    : DÉFORMATIONS
! OUT GRAD   : GRADIENT DES DÉPLACEMENTS
! IN  HEAVN  : DEFINITION DES DOMAINES DES FONCTIONS HEAVISIDES
!
    real(kind=8) :: zero, un, rac2, r, ur
    integer(kind=8) :: i, j, n, p, ig, cpt, nn, hea_se, alp
    real(kind=8) :: kron(3, 3), tmp, epstab(3, 3)
    aster_logical :: ldec
!
! ----------------------------------------------------------------------
!
!
! --- INITIALISATIONS
!
    zero = 0.d0
    un = 1.d0
    rac2 = sqrt(2.d0)
    hea_se = xcalc_code(nfiss, he_real=[he])
    ASSERT(nfe .le. 1)
    r = 0.d0
    ur = 0.d0
!
! --- MATRICE IDENTITE
!
    kron(:, :) = zero
    do p = 1, 3
        kron(p, p) = un
    end do
!
! --- CALCUL DES GRADIENTS : GRAD(U) ET F
!
    do j = 1, 3
        do i = 1, 3
            f(i, j) = kron(i, j)
        end do
    end do
!
    do j = 1, ndim
        do i = 1, ndim
            grad(i, j) = zero
        end do
    end do
!
    ldec = .false.
    if (ddlm .eq. 0 .or. ddlm .eq. -1 .or. ddlm .eq. ddls) ldec = .true.
!
! --- L'ORDRE DES DDLS DOIT ETRE 'DC' 'H1' 'E1' 'E2' 'E3' 'E4' 'LAGC'
!
    do n = 1, nnop
        if (ldec) then
! --- DDLM=-1 PERMET D'EVITER D'AVOIR A FOURNIR DDLM DANS CHAQUE CAS
            nn = ddls*(n-1)
        else
            call indent(n, ddls, ddlm, nnos, nn)
        end if
!
        cpt = 0
!
! -- DDLS CLASSIQUES
        do i = 1, ndim
            cpt = cpt+1
            do j = 1, ndim
                grad(i, j) = grad(i, j)+dfdi(n, j)*zr(idepl-1+nn+cpt)
            end do
        end do
        if (axi) then
            r = r+ff(n)*zr(igeom-1+2*(n-1)+1)
            ur = ur+ff(n)*zr(idepl-1+nn+1)
        end if
!
! -- DDLS HEAVISIDE
        do ig = 1, nfh
            do i = 1, ndim
                cpt = cpt+1
                do j = 1, ndim
                    grad(i, j) = grad(i, j)+xcalc_heav(heavn(n, ig), hea_se, heavn(n, 5))*dfdi(n, j&
                                &)*zr(idepl-1+nn+cpt)
                end do
            end do
            if (axi) then
                ur = ur+ff(n)*zr(idepl-1+nn+ndim*ig+1)*xcalc_heav(heavn(n, ig), hea_se, heavn&
                     &(n, 5))
            end if
        end do
!
! -- DDL ENRICHIS EN FOND DE FISSURE
        do ig = 1, nfe
            do alp = 1, ndim
                cpt = cpt+1
                do i = 1, ndim
                    do j = 1, ndim
                        grad(i, j) = grad(i, j)+zr(idepl-1+nn+cpt)*dkdgl(n, alp, i, j)
                    end do
                end do
                if (axi) then
                    ur = ur+zr(idepl-1+nn+cpt)*fk(n, alp, 1)
                end if
            end do
        end do
!
    end do
!
    if (axi) then
        ASSERT(r .gt. zero)
    end if
!
! --- CALCUL DES DÉFORMATIONS : EPS
!
    do i = 1, ndim
        do j = 1, i
            tmp = grad(i, j)+grad(j, i)
            epstab(i, j) = 0.5d0*tmp
        end do
    end do
    eps(:) = zero
    eps(1) = epstab(1, 1)
    eps(2) = epstab(2, 2)
    eps(4) = epstab(2, 1)*rac2
    if (ndim .eq. 3) then
        eps(3) = epstab(3, 3)
        eps(5) = epstab(3, 1)*rac2
        eps(6) = epstab(3, 2)*rac2
    else if (axi) then
        eps(3) = ur/r
    end if
!
end subroutine
