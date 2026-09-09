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
subroutine pipdef(typmod, &
                  ndim, nno, kpg, jv_poids, jv_vff, &
                  jv_dfde, geom, deplm, &
                  ddepl, depl0, depl1, &
                  epsm, deps_cst, deps_pil)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/nmgeom.h"
!
    character(len=8), intent(in) :: typmod(:)
    integer(kind=8), intent(in) :: ndim, nno, kpg
    integer(kind=8), intent(in) :: jv_poids, jv_vff, jv_dfde
    real(kind=8), intent(in) :: geom(:, :), deplm(:, :), ddepl(:, :), depl0(:, :), depl1(:, :)
    real(kind=8), intent(out) :: epsm(:), deps_cst(:), deps_pil(:)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE - PRED_ELAS/DEFORMATION)
!
! CALCUL DES DEFORMATIONS
!
! --------------------------------------------------------------------------------------------------
!
! in  typmod : type de modélisation
! in  ndim   : dimension de l'espace
! in  nno    : nombre de noeuds de l'element
! in  kpg    : numero du point de gauss
! in  jv_poids : poids des points de gauss
! in  jv_vff    : valeur des fonctions de forme
! in  jv_dfde  : derivee des fonctions de forme element de reference
! in  geom   : coordonnées des noeuds
! in  deplm  : deplacement en t-
! in  ddepl  : increment de deplacement a l'iteration newton courante
! in  depl0  : correction de deplacement pour forces fixes
! in  depl1  : correction de deplacement pour forces pilotees
! out epsm   : deformations au temps moins
! out deps_cst   : correction de deformations dues aux charges fixes
! out deps_pil   : correction de deformations dues aux charges pilotees
! --------------------------------------------------------------------------------------------------
    aster_logical :: axi
    real(kind=8) :: r, t9bid(3, 3), dfdi(nno, ndim)
    real(kind=8) :: poids
! --------------------------------------------------------------------------------------------------
    axi = typmod(1) .eq. 'AXIS'
!
! ----- eps(um)
    call nmgeom(ndim, nno, axi, ASTER_FALSE, geom, &
                kpg, jv_poids, jv_vff, jv_dfde, deplm, &
                ASTER_TRUE, poids, dfdi, t9bid, epsm, &
                r)

! ----- eps(du+du0)
    call nmgeom(ndim, nno, axi, ASTER_FALSE, geom, &
                kpg, jv_poids, jv_vff, jv_dfde, ddepl+depl0, &
                ASTER_FALSE, poids, dfdi, t9bid, deps_cst, &
                r)

! ----- eps(du1)
    call nmgeom(ndim, nno, axi, ASTER_FALSE, geom, &
                kpg, jv_poids, jv_vff, jv_dfde, depl1, &
                ASTER_FALSE, poids, dfdi, t9bid, deps_pil, &
                r)

!
end subroutine
