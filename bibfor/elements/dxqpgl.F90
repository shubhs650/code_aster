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
subroutine dxqpgl(xyzg, pgl)
!
    implicit none
!
    real(kind=8), intent(in) :: xyzg(3, *)
    real(kind=8), intent(out) :: pgl(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
! CONSTRUCTION DE LA MATRICE DE PASSAGE GLOBAL --> INTRINSEQUE
! POUR UNE MAILLE de plaque quadrangle
!
! --------------------------------------------------------------------------------------------------
!
!
!            I MILIEU DE 4 1                        3
!            J MILIEU DE 2 3                        *
!            K MILIEU DE 1 2                     L *  *
!            L MILIEU DE 3 4                      *     *
!                                                *        *
!        I : VECTEUR UNITAIRE PORTE PAR IJ    4 *           * J
!                                                *            *
!        K : PERPENDICULAIRE A IJ ET A KL        I*             *
!                                                  *              *
!        J : PRODUIT VECTORIEL K I                  *****************
!                                                  1        K        2
!
!
!     VERIFICATION QUE L'ELEMENT EST REELLEMENT PLAN
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: vx, vy, vz, xi, yi, zzi, xj
    real(kind=8) :: yj, zzj, xk, yk, zzk, xl, yl, zzl
    real(kind=8) :: norm
!
! --------------------------------------------------------------------------------------------------
!
    pgl = 0.d0
    xi = (xyzg(1, 1)+xyzg(1, 4))/2.d0
    yi = (xyzg(2, 1)+xyzg(2, 4))/2.d0
    zzi = (xyzg(3, 1)+xyzg(3, 4))/2.d0
    xj = (xyzg(1, 3)+xyzg(1, 2))/2.d0
    yj = (xyzg(2, 3)+xyzg(2, 2))/2.d0
    zzj = (xyzg(3, 3)+xyzg(3, 2))/2.d0
    xk = (xyzg(1, 2)+xyzg(1, 1))/2.d0
    yk = (xyzg(2, 2)+xyzg(2, 1))/2.d0
    zzk = (xyzg(3, 2)+xyzg(3, 1))/2.d0
    xl = (xyzg(1, 4)+xyzg(1, 3))/2.d0
    yl = (xyzg(2, 4)+xyzg(2, 3))/2.d0
    zzl = (xyzg(3, 4)+xyzg(3, 3))/2.d0
!
    norm = sqrt((xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zzj-zzi)*(zzj-zzi))
    pgl(1, 1) = (xj-xi)/norm
    pgl(1, 2) = (yj-yi)/norm
    pgl(1, 3) = (zzj-zzi)/norm
!
    vx = (yj-yi)*(zzl-zzk)-(zzj-zzi)*(yl-yk)
    vy = -(xj-xi)*(zzl-zzk)+(zzj-zzi)*(xl-xk)
    vz = (xj-xi)*(yl-yk)-(yj-yi)*(xl-xk)

    norm = sqrt(vx*vx+vy*vy+vz*vz)
    pgl(3, 1) = vx/norm
    pgl(3, 2) = vy/norm
    pgl(3, 3) = vz/norm
!
    pgl(2, 1) = pgl(3, 2)*pgl(1, 3)-pgl(3, 3)*pgl(1, 2)
    pgl(2, 2) = -pgl(3, 1)*pgl(1, 3)+pgl(3, 3)*pgl(1, 1)
    pgl(2, 3) = pgl(3, 1)*pgl(1, 2)-pgl(3, 2)*pgl(1, 1)
!
    norm = sqrt(pgl(2, 1)*pgl(2, 1)+pgl(2, 2)*pgl(2, 2)+pgl(2, 3)*pgl(2, 3))
    pgl(2, 1) = pgl(2, 1)/norm
    pgl(2, 2) = pgl(2, 2)/norm
    pgl(2, 3) = pgl(2, 3)/norm
!
end subroutine
