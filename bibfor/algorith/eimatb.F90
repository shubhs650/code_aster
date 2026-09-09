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
subroutine eimatb(nomte, ndim, axi, nno1, nno2, npg, &
                  wref, vff1, vff2, dffr2, geom, &
                  anglNautPg, b, wg, ni2ldc)

    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/eiinit.h"
#include "asterfort/matrot.h"
#include "asterfort/subaco.h"
#include "asterfort/sumetr.h"
!
    character(len=16) :: nomte
    aster_logical, intent(in) :: axi
    integer(kind=8), intent(in) :: ndim, nno1, nno2, npg
    real(kind=8), intent(in) :: vff1(nno1, npg), vff2(nno2, npg), geom(ndim, nno2)
    real(kind=8), intent(in) :: wref(npg)
    real(kind=8), intent(in) :: dffr2(ndim-1, nno2, npg)
    real(kind=8), intent(in) :: anglNautPg(merge(1, 3, ndim .eq. 2), nno2)
    real(kind=8), intent(out):: b(2*ndim, npg, 2*ndim*nno1+ndim*nno2)
    real(kind=8), intent(out):: wg(2*ndim, npg)
    real(kind=8), intent(out):: ni2ldc(2*ndim, npg)

!
! --------------------------------------------------------------------------------------------------
!
!   MATRICE B DDL -> DEFORMATIONS GENERALISEES POUR LES ELEMENTS D'INTERFACE
!
! --------------------------------------------------------------------------------------------------
! in  nomte   : nom de l'element
! in  ndim    : dimension de l'espace
! in  axi     : .true. si axisymetrie
! in  nno1    : nombre de noeuds (famille u)
! in  vff1    : valeur des fonctions de forme (famille u)
! in  nno2    : nombre de noeuds (famille x)
! in  vff2    : valeur des fonctions de forme (famille x)
! in  dffr2   : derivees des fonctions de forme de reference (famille l)
! in  npg     : nombre de points de gauss
! in  wref    : poids des points de gauss de reference
! in  geom    : coordonnees des noeuds
! out b       : matrice b
! out wg      : poids des points de gauss
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: n, i, j, nang, g
    integer(kind=8) :: iu(ndim, 2*nno1), il(ndim, nno2)
    real(kind=8) :: cova(3, 3), metr(2, 2), dfdx(9), cour, jac, cosa, sina
    real(kind=8) :: angloc(3), rot(3, 3), w
! --------------------------------------------------------------------------------------------------

    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
    nang = merge(1, 3, ndim .eq. 2)
    angloc = 0
    b = 0

    ! Decalage d'indice pour les elements d'interface
    call eiinit(nomte, iu, il)

    do g = 1, npg

        ! Jacobien
        if (ndim .eq. 3) then
            call subaco(nno2, dffr2(:, :, g), geom, cova)
            call sumetr(cova, metr, jac)
            w = wref(g)*jac
        else if (ndim .eq. 2) then
            call dfdm1d(nno2, wref(g), dffr2(:, :, g), geom, dfdx, &
                        cour, w, cosa, sina)
            if (axi) w = w*dot_product(geom(1, :), vff2(:, g))
        end if

        ! Contributions aux poids des points de Gauss
        wg(:, g) = w

        ! Angles nautiques au point d'integration
        angloc(1:nang) = matmul(anglNautPg, vff2(:, g))

        ! Matrice de rotation global -> local
        call matrot(angloc, rot)

        ! Contributions a la matrice B vis-a-vis de u
        forall (i=1:ndim, j=1:ndim, n=1:nno1)
            b(i, g, iu(j, n)) = -rot(i, j)*vff1(n, g)
            b(i, g, iu(j, n+nno1)) = rot(i, j)*vff1(n, g)
        end forall

        forall (i=1:ndim, n=1:nno2)
            b(i+ndim, g, il(i, n)) = vff2(n, g)
        end forall

    end do

    ! Profil des contraintes et des deformations : aucun correctif
    ni2ldc = 1

end subroutine
