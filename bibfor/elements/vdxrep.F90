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
subroutine vdxrep(plateOrie, nomte, epais, nodeCoor)
!
    use plate_type
    implicit none
!
#include "asterfort/jevete.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(inout) :: plateOrie
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: epais, nodeCoor(3, 9)
!
! --------------------------------------------------------------------------------------------------
!
! REMPLIR L'OBJET .DESR DANS LES ZONES 1090 ET 2000
! POUR POUVOIR CALCULER LES MATRICES DE PASSAGE AVEC VDREPE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ptType = 0
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: nb1, nb2, npgsr, i, j, k, intsr, lzi, lzr
    real(kind=8) ::  vectBaseKpg(3, 3)
    real(kind=8) :: vectNorm(9, 3), vectTang(9, 2, 3)
!
! --------------------------------------------------------------------------------------------------
!
    vectNorm = zero
    vectTang = zero

! - Access objects
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)

! - Compute local basis at nodes
    call vectan(nb1, nb2, &
                nodeCoor, zr(lzr), &
                vectNorm, vectTang)
    plateOrie%vectNorm = vectNorm
    plateOrie%vectTang = vectTang

! - Compute local basis at integration points
    k = 0
    do intsr = 1, npgsr
        call vectgt(plateOrie, ptType, nb1, &
                    nodeCoor, zero, intsr, &
                    epais, zr(lzr), &
                    vectBaseKpg)
        do j = 1, 3
            do i = 1, 3
                k = k+1
                zr(lzr+2000+k-1) = vectBaseKpg(i, j)
            end do
        end do
    end do
end subroutine
