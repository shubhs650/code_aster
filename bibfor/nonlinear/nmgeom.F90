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
subroutine nmgeom(ndim, nno, axi, grand, geom, &
                  kpg, ipoids, ivf, idfde, depl, &
                  ldfdi, poids, dfdi, f, eps, &
                  r)
    use tenseur_dime_module, only: NDIM_TO_NDIMSI, matrix_to_voigt, identity
    implicit none

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"

    integer(kind=8) :: ndim
    integer(kind=8) :: nno
    aster_logical, intent(in) :: axi
    aster_logical, intent(in) :: grand
    real(kind=8) :: geom(ndim, nno)
    integer(kind=8) :: kpg
    integer(kind=8) :: ipoids
    integer(kind=8) :: ivf
    integer(kind=8) :: idfde
    real(kind=8) :: depl(ndim, nno)
    aster_logical :: ldfdi
    real(kind=8) :: poids
    real(kind=8) :: dfdi(nno, ndim)
    real(kind=8), intent(out) :: f(3, 3)
    real(kind=8), intent(out) :: eps(:)
    real(kind=8), intent(out) :: r
!
! --------------------------------------------------------------------------------------------------
!     but:  calcul des elements cinematiques (matrices f et e, rayon r)
!           en un point de gauss (eventuellement en grandes transform.)
! --------------------------------------------------------------------------------------------------
! in  ndim    : dimension de l'espace
! in  nno     : nombre de noeuds de l'element
! in  axi     : indicateur si axisymetrique
! in  grand   : indicateur si grandes transformations
! in  geom    : coordonees des noeuds
! in  kpg     : numero du point de gauss (pour l'acces aux fct. formes)
! in  ipoids  : poids du point de gauss de l'element de reference
! in  ivf     : valeur des fonctions de forme (en axisymetrique)
! in  idfde   : derivee des fonctions de forme de reference
! in  depl    : deplacement a partir de la conf de ref
! in  depl    : deplacement a partir de la conf de ref
! in  ldfdi   : veut-on calculer dfdi et poids
! out poids   : "poids" du point de gauss
! out dfdi    : derivee des fonctions de forme
! out f       : gradient de la transformation (identité si hpp)
! out eps     : deformations
! out r       : distance du point de gauss a l'axe (en axisymetrique)
! --------------------------------------------------------------------------------------------------
! remarque concernant l'argument ldfdi :
!  nmgeom est parfois appele 2 fois de suite avec u et delta_u (par
!  exemple dans nmpl3d). comme le calcul de dfdm3d est couteux et qu'il
!  est independant de u, on peut economiser le 2eme calcul en utilisant
!  l'argument ldfdi : 1er appel .true. ; 2eme appel .false.
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: ndimsi
    real(kind=8) :: grad(3, 3), id33(3, 3), eps33(3, 3), ur
! --------------------------------------------------------------------------------------------------

    !  Initialisations
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
    ndimsi = NDIM_TO_NDIMSI(ndim)
    ASSERT(size(eps) .ge. ndimsi)
    id33 = identity(3)

    ! Calcul de la distance a l'axe (axisymetrique) et du depl. radial
    if (axi) then
        r = dot_product(zr(ivf+(kpg-1)*nno:ivf+kpg*nno-1), geom(1, :))
        ur = dot_product(zr(ivf+(kpg-1)*nno:ivf+kpg*nno-1), depl(1, :))
    end if

    ! Calcul des derivees des fonctions de forme et jacobien
    if (ldfdi) then
        if (ndim .eq. 3) then
            call dfdm3d(nno, kpg, ipoids, idfde, geom, &
                        poids, dfdi(1, 1), dfdi(1, 2), dfdi(1, 3))
        else
            call dfdm2d(nno, kpg, ipoids, idfde, geom, &
                        poids, dfdi(1, 1), dfdi(1, 2))
            if (axi) poids = poids*r
        end if
    end if

    ! Calcul des déformations
    grad = 0
    grad(1:ndim, 1:ndim) = matmul(depl, dfdi)
    if (axi) grad(3, 3) = ur/r

    if (grand) then
        f = id33+grad
        eps33 = 0.5d0*(matmul(transpose(f), f)-id33)
    else
        f = id33
        eps33 = 0.5d0*(grad+transpose(grad))
    end if

    ! Storage
    eps = 0
    eps(1:ndimsi) = matrix_to_voigt(eps33, ndimsi)

end subroutine
