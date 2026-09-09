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
subroutine dxtfor(plateOrie, global, xyzl, for, vecl)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/gtria3.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8) :: xyzl(3, *), for(6, *), vecl(*)
    aster_logical :: global
!
! --------------------------------------------------------------------------------------------------
!
!     CHARGEMENT FORCE_FACE DES ELEMENTS DE PLAQUE DKT ET DST
!
! --------------------------------------------------------------------------------------------------
!
!     IN  GLOBAL : VARIABLE LOGIQUE DE REPERE GLOBAL OU LOCAL
!     IN  XYZL   : COORDONNEES LOCALES DES TROIS NOEUDS
!     IN  FOR    : FORCE APPLIQUEE SUR LA FACE
!     OUT VECL   : CHARGEMENT NODAL RESULTANT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, nno
    real(kind=8) :: aire
    real(kind=8) :: fx, fy, carat3(21)
!
! --------------------------------------------------------------------------------------------------
!
    nno = 3

!     ----- CALCUL DES GRANDEURS GEOMETRIQUES SUR LE TRIANGLE ----------
    call gtria3(xyzl, carat3)

    if (.not. global) then
        do i = 1, nno
            fx = for(1, i)
            fy = for(2, i)
            for(1, i) = plateOrie%t2iu(1)*fx+plateOrie%t2iu(3)*fy
            for(2, i) = plateOrie%t2iu(2)*fx+plateOrie%t2iu(4)*fy
            fx = for(4, i)
            fy = for(5, i)
            for(4, i) = plateOrie%t2iu(1)*fx+plateOrie%t2iu(3)*fy
            for(5, i) = plateOrie%t2iu(2)*fx+plateOrie%t2iu(4)*fy
        end do
    end if
!
    aire = carat3(8)
!
    do i = 1, 6*nno
        vecl(i) = 0.d0
    end do
!
    do i = 1, 6
        vecl(i) = for(i, 1)*aire/3.d0
        vecl(i+6) = for(i, 2)*aire/3.d0
        vecl(i+12) = for(i, 3)*aire/3.d0
    end do
!
end subroutine
