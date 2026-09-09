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
subroutine vectan(nb1, nb2, &
                  nodeCoor, desr, &
                  vectNorm, vectTang)
!
    implicit none
!
    integer(kind=8), intent(in) :: nb1, nb2
    real(kind=8), intent(in) :: nodeCoor(3, *)
    real(kind=8), intent(inout) :: desr(*)
    real(kind=8), intent(out) :: vectNorm(9, 3), vectTang(9, 2, 3)
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Compute local basis at nodes
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: jvDFuncKsi = 828, jvDFuncEta = 900
    integer(kind=8) :: ib, i, j, k, l
    real(kind=8) :: rnorm, vecta(9, 2, 3)
!
! --------------------------------------------------------------------------------------------------
!
    do i = 1, nb2
! ----- Local tangents at nodes (manifold)
        do k = 1, 3
            vecta(i, 1, k) = 0.d0
            vecta(i, 2, k) = 0.d0
            do j = 1, nb1
                vecta(i, 1, k) = vecta(i, 1, k)+ &
                                 desr(jvDFuncKsi+8*(i-1)+j)*nodeCoor(k, j)
                vecta(i, 2, k) = vecta(i, 2, k)+ &
                                 desr(jvDFuncEta+8*(i-1)+j)*nodeCoor(k, j)
            end do
        end do

! ----- Local normal at nodes
        vectNorm(i, 1) = vecta(i, 1, 2)*vecta(i, 2, 3)-vecta(i, 1, 3)*vecta(i, 2, 2)
        vectNorm(i, 2) = vecta(i, 1, 3)*vecta(i, 2, 1)-vecta(i, 1, 1)*vecta(i, 2, 3)
        vectNorm(i, 3) = vecta(i, 1, 1)*vecta(i, 2, 2)-vecta(i, 1, 2)*vecta(i, 2, 1)
        rnorm = sqrt(vectNorm(i, 1)*vectNorm(i, 1)+ &
                     vectNorm(i, 2)*vectNorm(i, 2)+ &
                     vectNorm(i, 3)*vectNorm(i, 3))
        vectNorm(i, 1) = vectNorm(i, 1)/rnorm
        vectNorm(i, 2) = vectNorm(i, 2)/rnorm
        vectNorm(i, 3) = vectNorm(i, 3)/rnorm

! ----- Reconstruct orthornormals tangents
        rnorm = sqrt(vecta(i, 1, 1)*vecta(i, 1, 1)+ &
                     vecta(i, 1, 2)*vecta(i, 1, 2)+ &
                     vecta(i, 1, 3)*vecta(i, 1, 3))
        do k = 1, 3
            vectTang(i, 1, k) = vecta(i, 1, k)/rnorm
        end do
!
        vectTang(i, 2, 1) = vectNorm(i, 2)*vectTang(i, 1, 3)- &
                            vectNorm(i, 3)*vectTang(i, 1, 2)
        vectTang(i, 2, 2) = vectNorm(i, 3)*vectTang(i, 1, 1)- &
                            vectNorm(i, 1)*vectTang(i, 1, 3)
        vectTang(i, 2, 3) = vectNorm(i, 1)*vectTang(i, 1, 2)- &
                            vectNorm(i, 2)*vectTang(i, 1, 1)
    end do

!   STOCKAGE DES NB2 MATRICES DE PASSAGE LOCALES GLOBALE (3,3) DANS XR
    do ib = 1, nb2
        l = 9*(ib-1)
        do j = 1, 3
            do i = 1, 2
                k = l+(j-1)*3+i
                desr(1090+k) = vectTang(ib, i, j)
            end do
            k = l+(j-1)*3+3
            desr(1090+k) = vectNorm(ib, j)
        end do
    end do
!
end subroutine
