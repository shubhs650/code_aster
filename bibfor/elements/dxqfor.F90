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
subroutine dxqfor(plateOrie, global, xyzl, for, vecl)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/gquad4.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    aster_logical :: global
    real(kind=8) :: xyzl(3, *)
    real(kind=8) :: for(6, *)
    real(kind=8) :: vecl(*)
!
! --------------------------------------------------------------------------------------------------
!
!     CHARGEMENT FORCE_FACE DES ELEMENTS DE PLAQUE DKQ ET DSQ
!
! --------------------------------------------------------------------------------------------------
!
!     IN  GLOBAL : VARIABLE LOGIQUE DE REPERE GLOBAL OU LOCAL
!     IN  XYZL   : COORDONNEES LOCALES DES QUATRE NOEUDS
!     IN  PGL    : MATRICE DE PASSAGE GLOBAL - LOCAL
!     IN  FOR    : FORCE APPLIQUE SUR LA FACE
!     OUT VECL   : CHARGEMENT NODAL RESULTANT
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: airetr(4), c1, c2, fno(6, 4, 4)
    real(kind=8) :: fx, fy
    real(kind=8) :: caraq4(25)
    integer(kind=8) :: i, ino, it, j, k, nno
!
! --------------------------------------------------------------------------------------------------
!
    nno = 4
!
!     ----- CALCUL DES GRANDEURS GEOMETRIQUES SUR LE QUADRANGLE --------
    call gquad4(xyzl, caraq4)
!
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
    do ino = 1, nno
        airetr(ino) = caraq4(21+ino)
    end do
!
    do i = 1, 6
        do j = 1, nno
            do k = 1, nno
                fno(i, j, k) = 0.d0
            end do
        end do
    end do
!
    do i = 1, 6*nno
        vecl(i) = 0.d0
    end do
!
    c1 = 1.d0/6.d0
    c2 = 1.d0/12.d0
!
    do i = 1, 6
        fno(i, 1, 1) = (c1*for(i, 1)+c2*for(i, 2)+c2*for(i, 4))*airetr(1)
        fno(i, 1, 2) = (c2*for(i, 1)+c1*for(i, 2)+c2*for(i, 4))*airetr(1)
        fno(i, 1, 4) = (c2*for(i, 1)+c2*for(i, 2)+c1*for(i, 4))*airetr(1)
        fno(i, 2, 2) = (c1*for(i, 2)+c2*for(i, 3)+c2*for(i, 1))*airetr(2)
        fno(i, 2, 3) = (c2*for(i, 2)+c1*for(i, 3)+c2*for(i, 1))*airetr(2)
        fno(i, 2, 1) = (c2*for(i, 2)+c2*for(i, 3)+c1*for(i, 1))*airetr(2)
        fno(i, 3, 3) = (c1*for(i, 3)+c2*for(i, 4)+c2*for(i, 2))*airetr(3)
        fno(i, 3, 4) = (c2*for(i, 3)+c1*for(i, 4)+c2*for(i, 2))*airetr(3)
        fno(i, 3, 2) = (c2*for(i, 3)+c2*for(i, 4)+c1*for(i, 2))*airetr(3)
        fno(i, 4, 4) = (c1*for(i, 4)+c2*for(i, 1)+c2*for(i, 3))*airetr(4)
        fno(i, 4, 1) = (c2*for(i, 4)+c1*for(i, 1)+c2*for(i, 3))*airetr(4)
        fno(i, 4, 3) = (c2*for(i, 4)+c2*for(i, 1)+c1*for(i, 3))*airetr(4)
        do ino = 1, nno
            do it = 1, nno
                vecl(i+6*(ino-1)) = vecl(i+6*(ino-1))+fno(i, it, ino)
            end do
            vecl(i+6*(ino-1)) = vecl(i+6*(ino-1))/2.d0
        end do
    end do
!
end subroutine
