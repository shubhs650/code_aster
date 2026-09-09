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
subroutine calsta(gamma, dh, def, nno, &
                  kpg, sig, tmp, kk, kkd, &
                  matuu, dsidep, jac)
!
!     CALCUL DES TERMES DE STABILISATION POUR LE QUAD4 SOUS INTEGRE
!     STABILISE PAR LA METHODE ASSUMED STRAIN => QUAS4
!-----------------------------------------------------------------------
!
    implicit none
!
    integer(kind=8) :: kpg, kk, kkd, n, i, m, j, kl, nno, j1, ifiltr
    real(kind=8) :: dsidep(6, 6), f(3, 3)
    real(kind=8) :: tmp, sig(6)
    real(kind=8) :: gamma(4), dh(8)
    real(kind=8) :: kron(3, 3), kron2(2, 3), kron3(2, 3)
    real(kind=8) :: filtr1(2, 3), filtr2(2, 3), jac
    real(kind=8) :: matuu(*), defst1(4, 4, 2)
    real(kind=8) :: defst2(4, 4, 2), def(4, 4, 2)
    real(kind=8) :: rac2
!
    data kron/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
    data kron2/1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 1.d0/
    data kron3/1.d0, 0.d0, 0.d0, 1.d0, 1.d0, 0.d0/
!
!
!    PROJ : INDICATEUR DE LA PROJECTION
!           0 AUCUNE
!           1 OPTIMAL BENDING
!           2 INCOMPRESSIBLE
!
    rac2 = sqrt(2.d0)
!
    do i = 1, 3
        do j = 1, 3
            f(i, j) = kron(i, j)
        end do
    end do
!
!
!
    do i = 1, 2
        do j = 1, 3
            filtr1(i, j) = kron2(i, j)
            filtr2(i, j) = kron3(i, j)
        end do
    end do
!
    do ifiltr = 1, 3
        do n = 1, nno
            do i = 1, 2

!
                defst1(1, n, i) = f(i, 1)*gamma(n)*dh(2*kpg-1)*0.5d0*filtr1(1, ifiltr)+ &
                                  def(1, n, i)*filtr1(2, ifiltr)+ &
                                  f(i, 2)*(-0.5d0)*dh(2*kpg)*gamma(n)
!
                defst1(2, n, i) = f(i, 2)*gamma(n)*dh(2*kpg)*0.5d0*filtr1(1, ifiltr)+ &
                                  def(2, n, i)*filtr1(2, ifiltr)+ &
                                  f(i, 1)*(-0.5d0)*dh(2*kpg-1)*gamma(n)
!
                defst1(3, n, i) = 0.d0
!
                defst1(4, n, i) = def(4, n, i)*filtr1(2, ifiltr)
!
!
!
                defst2(1, n, i) = f(i, 1)*gamma(n)*dh(2*kpg-1)*0.5d0*filtr2(1, ifiltr)+ &
                                  def(1, n, i)*filtr2(2, ifiltr)+ &
                                  f(i, 2)*gamma(n)*dh(2*kpg)*(-0.5d0)
!
                defst2(2, n, i) = f(i, 2)*gamma(n)*dh(2*kpg)*0.5d0*filtr2(1, ifiltr)+ &
                                  def(2, n, i)*filtr2(2, ifiltr)+ &
                                  f(i, 1)*gamma(n)*dh(2*kpg-1)*(-0.5d0)
!
                defst2(3, n, i) = 0.d0
!
                defst2(4, n, i) = def(4, n, i)*filtr2(2, ifiltr)

            end do
        end do
!
!
! - CALCUL DE LA MATRICE DE RIGIDITE
!
!
        do n = 1, nno
            do i = 1, 2
                do kl = 1, 4
                    sig(kl) = 0.d0
                    sig(kl) = sig(kl)+defst1(1, n, i)*dsidep(1, kl)
                    sig(kl) = sig(kl)+defst1(2, n, i)*dsidep(2, kl)
                    sig(kl) = sig(kl)+defst1(3, n, i)*dsidep(3, kl)
                    sig(kl) = sig(kl)+defst1(4, n, i)*dsidep(4, kl)
                end do
!
                do j = 1, 2
                    do m = 1, n
                        if (m .eq. n) then
                            j1 = i
                        else
                            j1 = 2
                        end if
!
!                 RIGIDITE ELASTIQUE
                        tmp = 0.d0
                        tmp = tmp+sig(1)*defst2(1, m, j)
                        tmp = tmp+sig(2)*defst2(2, m, j)
                        tmp = tmp+sig(3)*defst2(3, m, j)
                        tmp = tmp+sig(4)*defst2(4, m, j)
!
!                 STOCKAGE EN TENANT COMPTE DE LA SYMETRIE
                        if (j .le. j1) then
                            kkd = (2*(n-1)+i-1)*(2*(n-1)+i)/2
                            kk = kkd+2*(m-1)+j
                            matuu(kk) = matuu(kk)+tmp*jac
                        end if
!
                    end do
                end do
            end do
        end do
    end do
!
end subroutine
