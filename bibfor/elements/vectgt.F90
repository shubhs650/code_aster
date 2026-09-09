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
subroutine vectgt(plateOrie, ptType, nb1, &
                  nodeCoor, ksi3s2, kpg, &
                  epais, desr, &
                  vectBaseKpg, &
                  vectTangKpg_)
!
    use plate_type
    implicit none
!
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: ptType, nb1
    real(kind=8), intent(in) :: nodeCoor(3, *), ksi3s2
    integer(kind=8), intent(in) :: kpg
    real(kind=8), intent(in) :: epais
    real(kind=8), intent(in) :: desr(*)
    real(kind=8), intent(out) :: vectBaseKpg(3, 3)
    real(kind=8), optional, intent(out) :: vectTangKpg_(2, 3)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i1, i2, j, k, l1, l2
    integer(kind=8) :: l3
    real(kind=8) :: rnorm, vectTangKpg(2, 3), vectNormKpg(3)
!
! --------------------------------------------------------------------------------------------------
!
    vectBaseKpg = 0.d0
    vectNormKpg = 0.d0
    vectTangKpg = 0.d0
    if (ptType .eq. 0) then
!     CALCULS AUX PTS D'INTEGRATION REDUITE
        l1 = 12
        l2 = 44
        l3 = 76
    else if (ptType .eq. 1) then
!     CALCULS AUX PTS D'INTEGRATION NORMALE
        l1 = 135
        l2 = 207
        l3 = 279
    end if

! - CONSTRUCTION DU VECTEUR N AUX X PTS DE GAUSS
    i1 = l1+8*(kpg-1)
    do k = 1, 3
        vectNormKpg(k) = 0
        do j = 1, nb1
            vectNormKpg(k) = vectNormKpg(k)+ &
                             desr(i1+j)*plateOrie%vectNorm(j, k)
        end do
    end do
    vectBaseKpg(3, :) = vectNormKpg(:)

! - CONSTRUCTION DES VECTEURS GA AUX X PTS DE GAUSS
    i1 = l2+8*(kpg-1)
    i2 = l3+8*(kpg-1)
    do k = 1, 3
        vectTangKpg(1, k) = 0.d0
        vectTangKpg(2, k) = 0.d0
        do j = 1, nb1
            vectTangKpg(1, k) = vectTangKpg(1, k)+ &
                                desr(i1+j)*(nodeCoor(k, j)+ksi3s2*epais*plateOrie%vectNorm(j, k))
            vectTangKpg(2, k) = vectTangKpg(2, k)+ &
                                desr(i2+j)*(nodeCoor(k, j)+ksi3s2*epais*plateOrie%vectNorm(j, k))
        end do
    end do
    rnorm = sqrt(vectTangKpg(1, 1)*vectTangKpg(1, 1)+ &
                 vectTangKpg(1, 2)*vectTangKpg(1, 2)+ &
                 vectTangKpg(1, 3)*vectTangKpg(1, 3))
!
    do k = 1, 3
        vectBaseKpg(1, k) = vectTangKpg(1, k)/rnorm
    end do
!
    vectBaseKpg(2, 1) = vectNormKpg(2)*vectBaseKpg(1, 3)-vectNormKpg(3)*vectBaseKpg(1, 2)
    vectBaseKpg(2, 2) = vectNormKpg(3)*vectBaseKpg(1, 1)-vectNormKpg(1)*vectBaseKpg(1, 3)
    vectBaseKpg(2, 3) = vectNormKpg(1)*vectBaseKpg(1, 2)-vectNormKpg(2)*vectBaseKpg(1, 1)

    if (present(vectTangKpg_)) then
        vectTangKpg_ = vectTangKpg
    end if
!
end subroutine
