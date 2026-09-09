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

subroutine nmbamb(ndim, nno, npg, geom, section, dxi_ff, w_ref, nddl, neps, b, w, ni2ldc)

    implicit none
#include "asterfort/dfdm1d_xd.h"

    integer(kind=8), intent(in) :: ndim, nno, npg
    real(kind=8), intent(in) :: geom(ndim, nno), dxi_ff(nno, npg), w_ref(npg), section
    integer(kind=8), intent(out) :: nddl, neps
    real(kind=8), intent(out), allocatable :: b(:, :, :)
    real(kind=8), intent(out), allocatable :: w(:, :), ni2ldc(:, :)
! -------------------------------------------------------------------------------------------------
!  calcul de la matrice B pour les éléments barre
! -------------------------------------------------------------------------------------------------
! in  ndim      dimension de l'espace
! in  nno       nombre de noeuds
! in  npg       nombre de points de gauss
! in  geom      coordonnees des noeuds
! in  section   aire de la section
! in  dxi_ff    dérivées des fonctions de forme dans l'élément de référence (nno,npg)
! in  w_ref     poids des pts de gauss de reference (npg)
! out nddl      nombre de ddl / element
! out neps      nbr de composante de deformation (generalisee)
! out b         matrice cinematique eps = b.u  (nddl, npg)
! out w         poids des points de gauss config initiale
! out ni2ldc    conversion contrainte stockee -> contrainte ldc
! -------------------------------------------------------------------------------------------------
    integer(kind=8) :: g, n, i, ddl
    real(kind=8) :: wg
    real(kind=8) :: ds_ff(nno), t(ndim)
! -------------------------------------------------------------------------------------------------

    nddl = nno*ndim
    neps = 1
    allocate (b(neps, npg, nddl), w(neps, npg), ni2ldc(neps, npg))

    do g = 1, npg

        ! Derivee des fonctions de forme no 2 (r et w non utilise)
        call dfdm1d_xd(dxi_ff(:, g), w_ref(g), geom, ds_ff, wg)

        ! Poids reel des points de Gauss (volume)
        w(:, g) = wg*section

        ! Vecteur tangent
        t = matmul(geom, ds_ff)

        ! Matrice B
        do concurrent(i=1:ndim)
            do concurrent(n=1:nno)
                ddl = (n-1)*ndim+i
                b(1, g, ddl) = t(i)*ds_ff(n)
            end do
        end do
    end do

    ! Fonction de transfert des contraintes utilisateur (N) -> ldc (sigma)
    ni2ldc(:, :) = 1.d0/section

end subroutine
