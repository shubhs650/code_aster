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
subroutine deflg2(gn, lamb, logl, pes, feta, &
                  xi, me)
    implicit none
!     CALCUL DES DES TERMES NECESSAIRES
!     AU POST TRAITEMENT DES CONTRAINTES
!     SUIVANT ARTICLE MIEHE APEL LAMBRECHT CMAME 2002
! ----------------------------------------------------------------------
!     IN   GN    directions propres du tenseur F
!     IN   LAMB  valeurs propres du tenseur F
!     IN   LOGL  log des valeurs propres du tenseur F
!     OUT  PES   tenseur P symetrise, pour le passage de T a PK2
!     OUT  FETA  utilitaires pour DEFLG3 : f_i=-2/lambda_i**2 puis eta
!     OUT  XI    utilitaires pour DEFLG3 : xi_ij
!     OUT  ME    utilitaires pour DEFLG3 : tenseur M d'ordre 4
! ----------------------------------------------------------------------
#include "asterc/r8miem.h"
#include "asterfort/symt46.h"
    real(kind=8) :: gn(3, 3), lamb(3), logl(3)
    real(kind=8) :: di(3), theta(3, 3)
    real(kind=8) :: pe(3, 3, 3, 3), me(3, 3, 3, 3), xi(3, 3), feta(4), pes(6, 6)
    real(kind=8) :: invla(3, 3), tt(3, 3)
    integer(kind=8) :: i, icas, j, k, l, a, b
! ----------------------------------------------------------------------
!
    feta = 0.d0
    do i = 1, 3
        di(i) = 1.d0/lamb(i)
        feta(i) = -2.d0/(lamb(i)*lamb(i))
    end do
!
    if (abs(lamb(1)-lamb(2)) .lt. r8miem()) then
        if (abs(lamb(1)-lamb(3)) .lt. r8miem()) then
            icas = 123
        else
            icas = 12
        end if
    else
        if (abs(lamb(1)-lamb(3)) .lt. r8miem()) then
            icas = 13
        else if (abs(lamb(2)-lamb(3)) .lt. r8miem()) then
            icas = 23
        else
            icas = 1
        end if
    end if
!
    theta = 0.d0
    xi = 0.d0
!
    if (icas .eq. 1) then
        invla(1, 2) = 1.d0/(lamb(1)-lamb(2))
        invla(1, 3) = 1.d0/(lamb(1)-lamb(3))
        invla(2, 3) = 1.d0/(lamb(2)-lamb(3))
        invla(2, 1) = -invla(1, 2)
        invla(3, 1) = -invla(1, 3)
        invla(3, 2) = -invla(2, 3)
!
        theta(1, 2) = (logl(1)-logl(2))*invla(1, 2)
        theta(2, 1) = theta(1, 2)
        theta(3, 2) = (logl(3)-logl(2))*invla(3, 2)
        theta(2, 3) = theta(3, 2)
        theta(1, 3) = (logl(1)-logl(3))*invla(1, 3)
        theta(3, 1) = theta(1, 3)
        xi(1, 2) = (theta(1, 2)-0.5d0*di(2))*invla(1, 2)
        xi(2, 1) = (theta(2, 1)-0.5d0*di(1))*invla(2, 1)
        xi(3, 2) = (theta(3, 2)-0.5d0*di(2))*invla(3, 2)
        xi(2, 3) = (theta(2, 3)-0.5d0*di(3))*invla(2, 3)
        xi(1, 3) = (theta(1, 3)-0.5d0*di(3))*invla(1, 3)
        xi(3, 1) = (theta(3, 1)-0.5d0*di(1))*invla(3, 1)

        feta(4) &
            = logl(1)*invla(1, 2)*invla(1, 3) &
              +logl(2)*invla(2, 1)*invla(2, 3) &
              +logl(3)*invla(3, 1)*invla(3, 2)

    else if (icas .eq. 12) then
        invla(1, 3) = 1.d0/(lamb(1)-lamb(3))
        invla(2, 3) = 1.d0/(lamb(2)-lamb(3))
        invla(3, 1) = -invla(1, 3)
        invla(3, 2) = -invla(2, 3)
!
        theta(1, 2) = di(1)*0.5d0
        theta(2, 1) = theta(1, 2)
        theta(1, 3) = (logl(1)-logl(3))*invla(1, 3)
        theta(3, 1) = theta(1, 3)
        theta(3, 2) = (logl(3)-logl(2))*invla(3, 2)
        theta(2, 3) = theta(3, 2)
        xi(1, 2) = feta(1)/8.d0
        xi(2, 1) = xi(1, 2)
        xi(3, 2) = (theta(3, 2)-0.5d0*di(2))*invla(3, 2)
        xi(2, 3) = (theta(2, 3)-0.5d0*di(3))*invla(2, 3)
        xi(1, 3) = (theta(1, 3)-0.5d0*di(3))*invla(1, 3)
        xi(3, 1) = (theta(3, 1)-0.5d0*di(1))*invla(3, 1)
        feta(4) = xi(3, 1)
    else if (icas .eq. 13) then
        invla(1, 2) = 1.d0/(lamb(1)-lamb(2))
        invla(2, 3) = 1.d0/(lamb(2)-lamb(3))
        invla(2, 1) = -invla(1, 2)
        invla(3, 2) = -invla(2, 3)
!
        theta(1, 2) = (logl(1)-logl(2))*invla(1, 2)
        theta(2, 1) = theta(1, 2)
        theta(1, 3) = di(1)*0.5d0
        theta(3, 1) = theta(1, 3)
        theta(3, 2) = (logl(3)-logl(2))*invla(3, 2)
        theta(2, 3) = theta(3, 2)
        xi(1, 2) = (theta(1, 2)-0.5d0*di(2))*invla(1, 2)
        xi(2, 1) = (theta(2, 1)-0.5d0*di(1))*invla(2, 1)
        xi(3, 2) = (theta(3, 2)-0.5d0*di(2))*invla(3, 2)
        xi(2, 3) = (theta(2, 3)-0.5d0*di(3))*invla(2, 3)
        xi(1, 3) = feta(1)/8.d0
        xi(3, 1) = xi(1, 3)
        feta(4) = xi(1, 2)
    else if (icas .eq. 23) then
        invla(1, 2) = 1.d0/(lamb(1)-lamb(2))
        invla(1, 3) = 1.d0/(lamb(1)-lamb(3))
        invla(2, 1) = -invla(1, 2)
        invla(3, 1) = -invla(1, 3)
!
        theta(1, 2) = (logl(1)-logl(2))*invla(1, 2)
        theta(2, 1) = theta(1, 2)
        theta(1, 3) = (logl(1)-logl(3))*invla(1, 3)
        theta(3, 1) = theta(1, 3)
        theta(3, 2) = di(1)*0.5d0
        theta(2, 3) = theta(3, 2)
        xi(1, 2) = (theta(1, 2)-0.5d0*di(2))*invla(1, 2)
        xi(2, 1) = (theta(2, 1)-0.5d0*di(1))*invla(2, 1)
        xi(3, 2) = feta(2)/8.d0
        xi(2, 3) = feta(2)/8.d0
        xi(1, 3) = (theta(1, 3)-0.5d0*di(3))*invla(1, 3)
        xi(3, 1) = (theta(3, 1)-0.5d0*di(1))*invla(3, 1)
        feta(4) = xi(2, 3)
    else if (icas .eq. 123) then
        theta(1, 2) = di(1)*0.5d0
        theta(2, 1) = theta(1, 2)
        theta(1, 3) = theta(1, 2)
        theta(3, 1) = theta(1, 2)
        theta(3, 2) = theta(1, 2)
        theta(2, 3) = theta(1, 2)
        xi(1, 2) = feta(2)/8.d0
        xi(2, 1) = xi(1, 2)
        xi(3, 2) = xi(1, 2)
        xi(2, 3) = xi(1, 2)
        xi(1, 3) = xi(1, 2)
        xi(3, 1) = xi(1, 2)
        feta(4) = xi(1, 2)
    end if
!
!
!     calcul de M_E (Lagrangien) pour chaque direction propre
    do i = 1, 3
        do j = 1, 3
            do a = 1, 3
                do b = 1, 3
                    me(a, b, i, j) = gn(a, i)*gn(b, j)+gn(a, j)*gn(b, i)
                end do
            end do
        end do
    end do
!
    do a = 1, 3
        do b = 1, 3
            tt(1, 2) = theta(1, 2)*me(a, b, 1, 2)
            tt(2, 1) = theta(2, 1)*me(a, b, 2, 1)
            tt(1, 3) = theta(1, 3)*me(a, b, 1, 3)
            tt(3, 1) = theta(3, 1)*me(a, b, 3, 1)
            tt(2, 3) = theta(2, 3)*me(a, b, 2, 3)
            tt(3, 2) = theta(3, 2)*me(a, b, 3, 2)
            tt(1, 1) = di(1)*me(a, b, 1, 1)*0.5d0
            tt(2, 2) = di(2)*me(a, b, 2, 2)*0.5d0
            tt(3, 3) = di(3)*me(a, b, 3, 3)*0.5d0
            do k = 1, 3
                do l = 1, 3
                    pe(k, l, a, b) &
                        = gn(k, 1)*(tt(1, 2)*gn(l, 2)+tt(1, 3)*gn(l, 3)+tt(1, 1)*gn(l, 1)) &
                          +gn(k, 2)*(tt(2, 1)*gn(l, 1)+tt(2, 3)*gn(l, 3)+tt(2, 2)*gn(l, 2)) &
                          +gn(k, 3)*(tt(3, 1)*gn(l, 1)+tt(3, 2)*gn(l, 2)+tt(3, 3)*gn(l, 3))
                end do
            end do
        end do
    end do
!
    call symt46(pe, pes)
!
end subroutine
