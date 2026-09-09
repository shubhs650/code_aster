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
subroutine deflg3(gn, feta, xi, me, t, &
                  tl)
    implicit none
!     CALCUL DES DEFORMATIONS LOGARITHMIQUES ET DES TERMES NECESSAIRES
!     AU POST TRAITEMENT DES CONTRAINTES ET A LA RIGIDITE TANGENTE
!     SUIVANT ARTICLE MIEHE APEL LAMBRECHT CMAME 2002
! ----------------------------------------------------------------------
!     IN GN    directions propres du tenseur F
!     IN FETA  utilitaires issus de DEFLG2  f_i=-2/lambda_i**2 puis eta
!     IN XI    utilitaires issus de DEFLG2  xi_ij
!     IN ME    utilitaires issus de DEFLG2  tenseur M d'ordre 4
!     IN T     tenseur des contraintesissu de NMCOMP (avec sqrt(2))
!     OUT TL   tenseur d'ordre 4 T:L
! ----------------------------------------------------------------------
#include "asterfort/tnsvec.h"
    real(kind=8) :: gn(3, 3), t(6), tl(3, 3, 3, 3)
    real(kind=8) :: dzeta(3, 3), t33(3, 3), me(3, 3, 3, 3), xi(3, 3), feta(4)
    real(kind=8) :: dzf(3, 3)
    integer(kind=8) :: a, b, c, d
! ----------------------------------------------------------------------
!
!     CALCUL DU TERME T.L
!
    dzeta = 0.d0
    call tnsvec(6, 3, t33, t, 1.d0/sqrt(2.d0))
!
!     A,B sont les composantes, J,I sont les modes propres

    do a = 1, 3
        do b = 1, 3
            dzeta(1, 1) = dzeta(1, 1)+t33(a, b)*gn(a, 1)*gn(b, 1)
            dzeta(1, 2) = dzeta(1, 2)+t33(a, b)*gn(a, 1)*gn(b, 2)
            dzeta(1, 3) = dzeta(1, 3)+t33(a, b)*gn(a, 1)*gn(b, 3)
            dzeta(2, 1) = dzeta(2, 1)+t33(a, b)*gn(a, 2)*gn(b, 1)
            dzeta(2, 2) = dzeta(2, 2)+t33(a, b)*gn(a, 2)*gn(b, 2)
            dzeta(2, 3) = dzeta(2, 3)+t33(a, b)*gn(a, 2)*gn(b, 3)
            dzeta(3, 1) = dzeta(3, 1)+t33(a, b)*gn(a, 3)*gn(b, 1)
            dzeta(3, 2) = dzeta(3, 2)+t33(a, b)*gn(a, 3)*gn(b, 2)
            dzeta(3, 3) = dzeta(3, 3)+t33(a, b)*gn(a, 3)*gn(b, 3)
        end do
    end do
!
    do a = 1, 3
        do b = 1, 3
            dzf(1, 1) = 0.25d0*feta(1)*dzeta(1, 1)*me(a, b, 1, 1)+ &
                        2.d0*xi(2, 1)*dzeta(2, 1)*me(a, b, 2, 1) &
                        +2.d0*xi(3, 1)*dzeta(3, 1)*me(a, b, 3, 1)
            dzf(2, 2) = 0.25d0*feta(2)*dzeta(2, 2)*me(a, b, 2, 2)+ &
                        2.d0*xi(1, 2)*dzeta(1, 2)*me(a, b, 1, 2) &
                        +2.d0*xi(3, 2)*dzeta(3, 2)*me(a, b, 3, 2)
            dzf(3, 3) = 0.25d0*feta(3)*dzeta(3, 3)*me(a, b, 3, 3)+ &
                        2.d0*xi(1, 3)*dzeta(1, 3)*me(a, b, 1, 3) &
                        +2.d0*xi(2, 3)*dzeta(2, 3)*me(a, b, 2, 3)
            dzf(2, 3) = 2.d0*feta(4)*dzeta(1, 2)*me(a, b, 1, 3)+ &
                        2.d0*xi(2, 3)*dzeta(2, 3)*me(a, b, 3, 3) &
                        +2.d0*xi(2, 3)*dzeta(3, 3)*me(a, b, 2, 3)
            dzf(3, 2) = 2.d0*feta(4)*dzeta(1, 3)*me(a, b, 1, 2)+ &
                        2.d0*xi(3, 2)*dzeta(3, 2)*me(a, b, 2, 2) &
                        +2.d0*xi(3, 2)*dzeta(2, 2)*me(a, b, 3, 2)
            dzf(1, 2) = 2.d0*feta(4)*dzeta(3, 1)*me(a, b, 3, 2)+ &
                        2.d0*xi(1, 2)*dzeta(1, 2)*me(a, b, 2, 2) &
                        +2.d0*xi(1, 2)*dzeta(2, 2)*me(a, b, 1, 2)
            dzf(2, 1) = 2.d0*feta(4)*dzeta(3, 2)*me(a, b, 3, 1)+ &
                        2.d0*xi(2, 1)*dzeta(2, 1)*me(a, b, 1, 1) &
                        +2.d0*xi(2, 1)*dzeta(1, 1)*me(a, b, 2, 1)
            dzf(1, 3) = 2.d0*feta(4)*dzeta(2, 1)*me(a, b, 2, 3)+ &
                        2.d0*xi(1, 3)*dzeta(1, 3)*me(a, b, 3, 3) &
                        +2.d0*xi(1, 3)*dzeta(3, 3)*me(a, b, 1, 3)
            dzf(3, 1) = 2.d0*feta(4)*dzeta(2, 3)*me(a, b, 2, 1)+ &
                        2.d0*xi(3, 1)*dzeta(3, 1)*me(a, b, 1, 1) &
                        +2.d0*xi(3, 1)*dzeta(1, 1)*me(a, b, 3, 1)
!
            do c = 1, 3
                do d = 1, 3
                    tl(a, b, c, d) &
                        = dzf(1, 1)*me(c, d, 1, 1) &
                          +dzf(2, 2)*me(c, d, 2, 2) &
                          +dzf(3, 3)*me(c, d, 3, 3) &
                          !
                          +dzf(2, 3)*me(c, d, 2, 3) &
                          +dzf(3, 2)*me(c, d, 3, 2) &
                          +dzf(1, 3)*me(c, d, 1, 3) &
                          +dzf(3, 1)*me(c, d, 3, 1) &
                          +dzf(1, 2)*me(c, d, 1, 2) &
                          +dzf(2, 1)*me(c, d, 2, 1)
                end do
            end do
        end do
    end do

!
end subroutine
