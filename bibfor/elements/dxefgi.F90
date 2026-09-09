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
subroutine dxefgi(plateCara, plateOrie, &
                  npg, epsini, efge)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dxmate.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: npg
    real(kind=8), intent(in) :: epsini(6)
    real(kind=8), intent(out) :: efge(*)
!
! --------------------------------------------------------------------------------------------------
!
! --- EFFORTS GENERALISES D'ORIGINE THERMIQUE AUX POINTS D'INTEGRATION
! --- POUR LES ELEMENTS COQUES A FACETTES PLANES :
! --- DST, DKT, DSQ, DKQ, Q4G
! --- CALCULES A PARTIR D'UN CHAMP DE DEFORMATIONS INITIALES QUI EST
! --- POUR L'INSTANT CONSTANT PAR ELEMENT ET QUI NE PREND PAS EN
! --- COMPTE LES DEFORMATIONS INITIALES DE CISAILLEMENT TRANSVERSE.
!
! --------------------------------------------------------------------------------------------------
!
!     IN  NOMTE        : NOM DU TYPE D'ELEMENT
!     IN  PGL(3,3)     : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                        LOCAL
!     IN  EPSINI(6)    : DEFORMATIONS INITIALES CONSTANTES SUR L'ELEMENT
!                        DANS L'ORDRE : EPXX, EPYY, EPXY, KXX, KYY, KXY
!     OUT SIGT(1)      : EFFORTS  GENERALISES D'ORIGINE THERMIQUE
!                        AUX POINTS D'INTEGRATION
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8), parameter :: ncomp = 6
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: kxx, kyy, kxy
    aster_logical :: coupmf
    integer(kind=8) :: kpg, multic
    real(kind=8) :: epxx, epxy, epyy
!
! --------------------------------------------------------------------------------------------------
!
    efge(1:32) = zero

! - Get elementary matrix of rigidity
    call dxmate(plateCara, plateOrie, &
                'RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)

    do kpg = 1, npg
        epxx = epsini(ncomp*(kpg-1)+1)
        epyy = epsini(ncomp*(kpg-1)+2)
        epxy = 2.d0*epsini(ncomp*(kpg-1)+3)
        kxx = epsini(ncomp*(kpg-1)+4)
        kyy = epsini(ncomp*(kpg-1)+5)
        kxy = 2.d0*epsini(ncomp*(kpg-1)+6)
        efge(1+8*(kpg-1)) = dm(1, 1)*epxx+dm(1, 2)*epyy+dm(1, 3)*epxy
        efge(2+8*(kpg-1)) = dm(2, 1)*epxx+dm(2, 2)*epyy+dm(2, 3)*epxy
        efge(3+8*(kpg-1)) = dm(3, 1)*epxx+dm(3, 2)*epyy+dm(3, 3)*epxy
        efge(4+8*(kpg-1)) = df(1, 1)*kxx+df(1, 2)*kyy+df(1, 3)*kxy
        efge(5+8*(kpg-1)) = df(2, 1)*kxx+df(2, 2)*kyy+df(2, 3)*kxy
        efge(6+8*(kpg-1)) = df(3, 1)*kxx+df(3, 2)*kyy+df(3, 3)*kxy
        efge(7+8*(kpg-1)) = zero
        efge(8+8*(kpg-1)) = zero
    end do
!
end subroutine
