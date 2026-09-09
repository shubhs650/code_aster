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
subroutine nmspci(nno_p, nno_s, vff_p, vff_s, b)
!
!
    implicit none
#include "blas/dcopy.h"
!
    integer(kind=8) :: nno_p, nno_s
    real(kind=8) :: vff_p(18), vff_s(nno_s)
    real(kind=8) :: b(3, nno_s*3+nno_p*6)
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!  MATRICE CINEMATIQUE POUR 3D_INTERF_POU (EN UN POINT DE GAUSS DONNE)
!-----------------------------------------------------------------------
! IN  NNO_P  NOMBRE DE NOEUDS DE L'ELEMENT DE POUTRE
! IN  NNO_S  NOMBRE DE NOEUDS DE L'ELEMENT DE SOL
! IN  VFF_P  VALEUR DES FONCTIONS DE FORME L'ELEMENT DE POUTRE
! IN  VFF_S  VALEUR DES FONCTIONS DE FORME L'ELEMENT DE SOL
! OUT B      MATRICE DE PASSAGE UNODAL -> SAUT DE U LOCAL
!-----------------------------------------------------------------------
    integer(kind=8) :: n, nddl_s
!-----------------------------------------------------------------------
!
    nddl_s = nno_s*3
    b(:, :) = 0.d0

!   Contribution du sol
    b_n = to_blas_int(nno_s)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(9)
    do n = 1, 3
        call dcopy(b_n, vff_s, b_incx, b(n, n), b_incy)
    end do

!   Contribution de la poutre
    b(1, nddl_s+(/1, 7/)) = -vff_p((/1, 4/))
    b(2, nddl_s+(/2, 6, 8, 12/)) = -vff_p((/2, 3, 5, 6/))
    b(3, nddl_s+(/3, 5, 9, 11/)) = -vff_p(11:14)

end subroutine
