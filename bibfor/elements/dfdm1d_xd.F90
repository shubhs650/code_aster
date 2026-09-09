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

subroutine dfdm1d_xd(dxi_ff, w_ref, geom, ds_ff, w)

    implicit none
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"

    real(kind=8), intent(in) :: dxi_ff(:), w_ref, geom(:, :)
    real(kind=8), intent(out):: ds_ff(:), w
! -------------------------------------------------------------------------------------------------
! Calcul des dérivées des fonctions de forme par rapport à l'abscisse
! curviligne et le poids des points de Gauss pour un élément 1D plongé
! dans un espace 2D ou 3D, pour le point de Gauss courant
! -------------------------------------------------------------------------------------------------
! dxi_ff    in  dérivées des fonctions de forme de l'élément de référence
! w_ref     in  poids des points de Gauss dans l'élément de référence
! geom      in  coordonnées des noeuds (ndim,nno)
! ds_ff     out dérivée des fonctions de forme par rapport à l'abscisse curviligne
! w         out poids des points de Gauss de l'élément réel
! -------------------------------------------------------------------------------------------------
    integer(kind=8) :: iadzi, iazk24, numail
    real(kind=8)    :: a(size(geom, 1)), jac
! -------------------------------------------------------------------------------------------------

    ! Vecteur tangent
    a = matmul(geom, dxi_ff)

    ! Jacobien de la transformation
    jac = norm2(a)

    if (abs(jac) .le. 1.d0/r8gaem()) then
        call tecael(iadzi, iazk24, 0)
        numail = zi(iadzi)
        call utmess('F', 'ALGORITH2_59', si=numail)
    end if

    ! Dérivées des fonctions de forme par rapport à l'abscisse curviligne
    ds_ff = dxi_ff/jac

    ! Poids du point de Gauss
    w = jac*w_ref

end subroutine
