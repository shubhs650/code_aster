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
subroutine vdefge(nomte, nb1, npgsr, xr, epais, &
                  sigmElno, efgeElno)
!
    implicit none
!
#include "asterfort/assert.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: nomte
    integer(kind=8), intent(in) :: nb1, npgsr
    real(kind=8), intent(in) :: xr(*), epais
    real(kind=8), intent(in) :: sigmElno(6, MT_NNOMAX3D)
    real(kind=8), intent(out) :: efgeElno(8, 9)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: l1 = 1452
    real(kind=8), parameter :: wnc1 = 0.33333333333333d0
    real(kind=8), parameter :: wnc2 = 1.33333333333333d0
    real(kind=8), parameter :: wnc3 = 0.33333333333333d0
    integer(kind=8) :: i, i1, k
    real(kind=8) :: effgtg(8, 8)
    real(kind=8) :: demiep, zic, zic1, zic2, zic3
!
! --------------------------------------------------------------------------------------------------
!
!     DANS LE CAS DE EFGE_ELNO, ON CALCULE AUX NB1 NOEUDS
!     ON NE CALCULE PAS AU NOEUD INTERNE
!
    efgeElno = 0.d0
    demiep = epais/2.d0
    zic = -epais/2.d0
!
    zic1 = zic
    zic2 = zic1+demiep
    zic3 = zic2+demiep
!
    if (nomte .eq. 'MEC3QU9H') then
        effgtg(1, 1) = demiep*(wnc1*sigmElno(1, 1)+wnc2*sigmElno(1, 5)+wnc3*sigmElno(1, 9))
        effgtg(1, 2) = demiep*(wnc1*sigmElno(1, 2)+wnc2*sigmElno(1, 6)+wnc3*sigmElno(1, 10))
        effgtg(1, 3) = demiep*(wnc1*sigmElno(1, 3)+wnc2*sigmElno(1, 7)+wnc3*sigmElno(1, 11))
        effgtg(1, 4) = demiep*(wnc1*sigmElno(1, 4)+wnc2*sigmElno(1, 8)+wnc3*sigmElno(1, 12))
!
        effgtg(2, 1) = demiep*(wnc1*sigmElno(2, 1)+wnc2*sigmElno(2, 5)+wnc3*sigmElno(2, 9))
        effgtg(2, 2) = demiep*(wnc1*sigmElno(2, 2)+wnc2*sigmElno(2, 6)+wnc3*sigmElno(2, 10))
        effgtg(2, 3) = demiep*(wnc1*sigmElno(2, 3)+wnc2*sigmElno(2, 7)+wnc3*sigmElno(2, 11))
        effgtg(2, 4) = demiep*(wnc1*sigmElno(2, 4)+wnc2*sigmElno(2, 8)+wnc3*sigmElno(2, 12))
!
        effgtg(3, 1) = demiep*(wnc1*sigmElno(4, 1)+wnc2*sigmElno(4, 5)+wnc3*sigmElno(4, 9))
        effgtg(3, 2) = demiep*(wnc1*sigmElno(4, 2)+wnc2*sigmElno(4, 6)+wnc3*sigmElno(4, 10))
        effgtg(3, 3) = demiep*(wnc1*sigmElno(4, 3)+wnc2*sigmElno(4, 7)+wnc3*sigmElno(4, 11))
        effgtg(3, 4) = demiep*(wnc1*sigmElno(4, 4)+wnc2*sigmElno(4, 8)+wnc3*sigmElno(4, 12))
!
        effgtg(4, 1) = demiep*(wnc1*zic1*sigmElno(1, 1)+wnc2*zic2*sigmElno(1, 5)+ &
                               wnc3*zic3*sigmElno(1, 9))
        effgtg(4, 2) = demiep*(wnc1*zic1*sigmElno(1, 2)+wnc2*zic2*sigmElno(1, 6)+ &
                               wnc3*zic3*sigmElno(1, 10))
        effgtg(4, 3) = demiep*(wnc1*zic1*sigmElno(1, 3)+wnc2*zic2*sigmElno(1, 7)+ &
                               wnc3*zic3*sigmElno(1, 11))
        effgtg(4, 4) = demiep*(wnc1*zic1*sigmElno(1, 4)+wnc2*zic2*sigmElno(1, 8)+ &
                               wnc3*zic3*sigmElno(1, 12))
!
        effgtg(5, 1) = demiep*(wnc1*zic1*sigmElno(2, 1)+wnc2*zic2*sigmElno(2, 5)+ &
                               wnc3*zic3*sigmElno(2, 9))
        effgtg(5, 2) = demiep*(wnc1*zic1*sigmElno(2, 2)+wnc2*zic2*sigmElno(2, 6)+ &
                               wnc3*zic3*sigmElno(2, 10))
        effgtg(5, 3) = demiep*(wnc1*zic1*sigmElno(2, 3)+wnc2*zic2*sigmElno(2, 7)+ &
                               wnc3*zic3*sigmElno(2, 11))
        effgtg(5, 4) = demiep*(wnc1*zic1*sigmElno(2, 4)+wnc2*zic2*sigmElno(2, 8)+ &
                               wnc3*zic3*sigmElno(2, 12))
!
        effgtg(6, 1) = demiep*(wnc1*zic1*sigmElno(4, 1)+wnc2*zic2*sigmElno(4, 5)+ &
                               wnc3*zic3*sigmElno(4, 9))
        effgtg(6, 2) = demiep*(wnc1*zic1*sigmElno(4, 2)+wnc2*zic2*sigmElno(4, 6)+ &
                               wnc3*zic3*sigmElno(4, 10))
        effgtg(6, 3) = demiep*(wnc1*zic1*sigmElno(4, 3)+wnc2*zic2*sigmElno(4, 7)+ &
                               wnc3*zic3*sigmElno(4, 11))
        effgtg(6, 4) = demiep*(wnc1*zic1*sigmElno(4, 4)+wnc2*zic2*sigmElno(4, 8)+ &
                               wnc3*zic3*sigmElno(4, 12))
!
        effgtg(7, 1) = demiep*(wnc1*sigmElno(5, 1)+wnc2*sigmElno(5, 5)+wnc3*sigmElno(5, 9))
        effgtg(7, 2) = demiep*(wnc1*sigmElno(5, 2)+wnc2*sigmElno(5, 6)+wnc3*sigmElno(5, 10))
        effgtg(7, 3) = demiep*(wnc1*sigmElno(5, 3)+wnc2*sigmElno(5, 7)+wnc3*sigmElno(5, 11))
        effgtg(7, 4) = demiep*(wnc1*sigmElno(5, 4)+wnc2*sigmElno(5, 8)+wnc3*sigmElno(5, 12))
!
        effgtg(8, 1) = demiep*(wnc1*sigmElno(6, 1)+wnc2*sigmElno(6, 5)+wnc3*sigmElno(6, 9))
        effgtg(8, 2) = demiep*(wnc1*sigmElno(6, 2)+wnc2*sigmElno(6, 6)+wnc3*sigmElno(6, 10))
        effgtg(8, 3) = demiep*(wnc1*sigmElno(6, 3)+wnc2*sigmElno(6, 7)+wnc3*sigmElno(6, 11))
        effgtg(8, 4) = demiep*(wnc1*sigmElno(6, 4)+wnc2*sigmElno(6, 8)+wnc3*sigmElno(6, 12))
!
        do i = 1, nb1
            i1 = l1+4*(i-1)
            do k = 1, npgsr
                efgeElno(1, i) = efgeElno(1, i)+effgtg(1, k)*xr(i1+k)
                efgeElno(2, i) = efgeElno(2, i)+effgtg(2, k)*xr(i1+k)
                efgeElno(3, i) = efgeElno(3, i)+effgtg(3, k)*xr(i1+k)
                efgeElno(4, i) = efgeElno(4, i)+effgtg(4, k)*xr(i1+k)
                efgeElno(5, i) = efgeElno(5, i)+effgtg(5, k)*xr(i1+k)
                efgeElno(6, i) = efgeElno(6, i)+effgtg(6, k)*xr(i1+k)
                efgeElno(7, i) = efgeElno(7, i)+effgtg(7, k)*xr(i1+k)
                efgeElno(8, i) = efgeElno(8, i)+effgtg(8, k)*xr(i1+k)
!
            end do
        end do
!
!     VALEURS AU NOEUD INTERNE OBTENUE PAR MOYENNE DES AUTRES
!
        efgeElno(1, 9) = (efgeElno(1, 5)+efgeElno(1, 6)+efgeElno(1, 7)+efgeElno(1, 8))/4.d0
        efgeElno(2, 9) = (efgeElno(2, 5)+efgeElno(2, 6)+efgeElno(2, 7)+efgeElno(2, 8))/4.d0
        efgeElno(3, 9) = (efgeElno(3, 5)+efgeElno(3, 6)+efgeElno(3, 7)+efgeElno(3, 8))/4.d0
        efgeElno(4, 9) = (efgeElno(4, 5)+efgeElno(4, 6)+efgeElno(4, 7)+efgeElno(4, 8))/4.d0
        efgeElno(5, 9) = (efgeElno(5, 5)+efgeElno(5, 6)+efgeElno(5, 7)+efgeElno(5, 8))/4.d0
        efgeElno(6, 9) = (efgeElno(6, 5)+efgeElno(6, 6)+efgeElno(6, 7)+efgeElno(6, 8))/4.d0
        efgeElno(7, 9) = (efgeElno(7, 5)+efgeElno(7, 6)+efgeElno(7, 7)+efgeElno(7, 8))/4.d0
        efgeElno(8, 9) = (efgeElno(8, 5)+efgeElno(8, 6)+efgeElno(8, 7)+efgeElno(8, 8))/4.d0
!
    else if (nomte .eq. 'MEC3TR7H') then
!
        effgtg(1, 1) = demiep*(wnc1*sigmElno(1, 1)+wnc2*sigmElno(1, 4)+wnc3*sigmElno(1, 7))
        effgtg(1, 2) = demiep*(wnc1*sigmElno(1, 2)+wnc2*sigmElno(1, 5)+wnc3*sigmElno(1, 8))
        effgtg(1, 3) = demiep*(wnc1*sigmElno(1, 3)+wnc2*sigmElno(1, 6)+wnc3*sigmElno(1, 9))
!
        effgtg(2, 1) = demiep*(wnc1*sigmElno(2, 1)+wnc2*sigmElno(2, 4)+wnc3*sigmElno(2, 7))
        effgtg(2, 2) = demiep*(wnc1*sigmElno(2, 2)+wnc2*sigmElno(2, 5)+wnc3*sigmElno(2, 8))
        effgtg(2, 3) = demiep*(wnc1*sigmElno(2, 3)+wnc2*sigmElno(2, 6)+wnc3*sigmElno(2, 9))
!
        effgtg(3, 1) = demiep*(wnc1*sigmElno(4, 1)+wnc2*sigmElno(4, 4)+wnc3*sigmElno(4, 7))
        effgtg(3, 2) = demiep*(wnc1*sigmElno(4, 2)+wnc2*sigmElno(4, 5)+wnc3*sigmElno(4, 8))
        effgtg(3, 3) = demiep*(wnc1*sigmElno(4, 3)+wnc2*sigmElno(4, 6)+wnc3*sigmElno(4, 9))
!
        effgtg(4, 1) = demiep*(wnc1*zic1*sigmElno(1, 1)+wnc2*zic2*sigmElno(1, 4)+ &
                               wnc3*zic3*sigmElno(1, 7))
        effgtg(4, 2) = demiep*(wnc1*zic1*sigmElno(1, 2)+wnc2*zic2*sigmElno(1, 5)+ &
                               wnc3*zic3*sigmElno(1, 8))
        effgtg(4, 3) = demiep*(wnc1*zic1*sigmElno(1, 3)+wnc2*zic2*sigmElno(1, 6)+ &
                               wnc3*zic3*sigmElno(1, 9))
!
        effgtg(5, 1) = demiep*(wnc1*zic1*sigmElno(2, 1)+wnc2*zic2*sigmElno(2, 4)+ &
                               wnc3*zic3*sigmElno(2, 7))
        effgtg(5, 2) = demiep*(wnc1*zic1*sigmElno(2, 2)+wnc2*zic2*sigmElno(2, 5)+ &
                               wnc3*zic3*sigmElno(2, 8))
        effgtg(5, 3) = demiep*(wnc1*zic1*sigmElno(2, 3)+wnc2*zic2*sigmElno(2, 6)+ &
                               wnc3*zic3*sigmElno(2, 9))
!
        effgtg(6, 1) = demiep*(wnc1*zic1*sigmElno(4, 1)+wnc2*zic2*sigmElno(4, 4)+ &
                               wnc3*zic3*sigmElno(4, 7))
        effgtg(6, 2) = demiep*(wnc1*zic1*sigmElno(4, 2)+wnc2*zic2*sigmElno(4, 5)+ &
                               wnc3*zic3*sigmElno(4, 8))
        effgtg(6, 3) = demiep*(wnc1*zic1*sigmElno(4, 3)+wnc2*zic2*sigmElno(4, 6)+ &
                               wnc3*zic3*sigmElno(4, 9))
!
        effgtg(7, 1) = demiep*(wnc1*sigmElno(5, 1)+wnc2*sigmElno(5, 4)+wnc3*sigmElno(5, 7))
        effgtg(7, 2) = demiep*(wnc1*sigmElno(5, 2)+wnc2*sigmElno(5, 5)+wnc3*sigmElno(5, 8))
        effgtg(7, 3) = demiep*(wnc1*sigmElno(5, 3)+wnc2*sigmElno(5, 6)+wnc3*sigmElno(5, 9))
!
        effgtg(8, 1) = demiep*(wnc1*sigmElno(6, 1)+wnc2*sigmElno(6, 4)+wnc3*sigmElno(6, 7))
        effgtg(8, 2) = demiep*(wnc1*sigmElno(6, 2)+wnc2*sigmElno(6, 5)+wnc3*sigmElno(6, 8))
        effgtg(8, 3) = demiep*(wnc1*sigmElno(6, 3)+wnc2*sigmElno(6, 6)+wnc3*sigmElno(6, 9))
!
        do i = 1, nb1
            i1 = l1+4*(i-1)
            do k = 1, npgsr
                efgeElno(1, i) = efgeElno(1, i)+effgtg(1, k)*xr(i1+k)
                efgeElno(2, i) = efgeElno(2, i)+effgtg(2, k)*xr(i1+k)
                efgeElno(3, i) = efgeElno(3, i)+effgtg(3, k)*xr(i1+k)
                efgeElno(4, i) = efgeElno(4, i)+effgtg(4, k)*xr(i1+k)
                efgeElno(5, i) = efgeElno(5, i)+effgtg(5, k)*xr(i1+k)
                efgeElno(6, i) = efgeElno(6, i)+effgtg(6, k)*xr(i1+k)
                efgeElno(7, i) = efgeElno(7, i)+effgtg(7, k)*xr(i1+k)
                efgeElno(8, i) = efgeElno(8, i)+effgtg(8, k)*xr(i1+k)
!
            end do
        end do
!
!     VALEURS AU NOEUD INTERNE OBTENUE PAR MOYENNE DES AUTRES
!
        efgeElno(1, 7) = (efgeElno(1, 1)+efgeElno(1, 2)+efgeElno(1, 3))/3.d0
        efgeElno(2, 7) = (efgeElno(2, 1)+efgeElno(2, 2)+efgeElno(2, 3))/3.d0
        efgeElno(3, 7) = (efgeElno(3, 1)+efgeElno(3, 2)+efgeElno(3, 3))/3.d0
        efgeElno(4, 7) = (efgeElno(4, 1)+efgeElno(4, 2)+efgeElno(4, 3))/3.d0
        efgeElno(5, 7) = (efgeElno(5, 1)+efgeElno(5, 2)+efgeElno(5, 3))/3.d0
        efgeElno(6, 7) = (efgeElno(6, 1)+efgeElno(6, 2)+efgeElno(6, 3))/3.d0
        efgeElno(7, 7) = (efgeElno(7, 1)+efgeElno(7, 2)+efgeElno(7, 3))/3.d0
        efgeElno(8, 7) = (efgeElno(8, 1)+efgeElno(8, 2)+efgeElno(8, 3))/3.d0
    else
        ASSERT(ASTER_FALSE)

    end if
!
end subroutine
