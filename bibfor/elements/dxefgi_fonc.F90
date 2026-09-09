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
subroutine dxefgi_fonc(plateCara, plateOrie, &
                       nno, npg, &
                       epsinif, xyz, ni, efge)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dxmate.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: nno, npg
    character(len=8), intent(in) :: epsinif(6)
    real(kind=8), intent(in) :: xyz(*), ni(*)
    real(kind=8), intent(out) :: efge(*)
!
! --------------------------------------------------------------------------------------------------
!
! --- EFFORTS GENERALISES DE DEFORMATION INITIALE AUX POINTS D'INTEGRATION
! --- POUR LES ELEMENTS COQUES A FACETTES PLANES :
! --- DST, DKT, DSQ, DKQ, Q4G
! --- CALCULES A PARTIR D'UN CHAMP DE DEFORMATIONS INITIALES SOUS FORME
! --- DE FONCTION ET QUI NE PREND PAS EN
! --- COMPTE LES DEFORMATIONS INITIALES DE CISAILLEMENT TRANSVERSE.
!
! --------------------------------------------------------------------------------------------------
!
!     IN  NOMTE        : NOM DU TYPE D'ELEMENT
!     IN  PGL(3,3)     : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                        LOCAL
!     IN  EPSINIF(6)   : FONCTIONS DE DEFORMATIONS INITIALES
!                        DANS L'ORDRE : EPXX, EPYY, EPXY, KXX, KYY, KXY
!     XYZ              : COORDONNEES DES CONNECTIVITES
!     NI               : FONCTIONS DE FORME
!     OUT SIGT(1)      : EFFORTS  GENERALISES D'ORIGINE THERMIQUE
!                        AUX POINTS D'INTEGRATION
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8), parameter :: nbPara = 4
    character(len=8), parameter :: paraName(nbPara) = (/"X   ", "Y   ", "Z   ", "INST"/)
    real(kind=8) :: paraVale(nbPara)
    integer(kind=8) :: multic, jvInstr, paraCode
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    aster_logical :: coupmf
    real(kind=8) :: xgau, ygau, zgau
    integer(kind=8) :: i, kpg
    real(kind=8) :: epxx, epxy, epyy, kxx, kyy, kxy
!
! --------------------------------------------------------------------------------------------------
!
    efge(1:32) = zero

! - Get elementary matrix of rigidity
    call dxmate(plateCara, plateOrie, &
                'RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)

! - Get current time
    call jevech('PINSTR', 'L', jvInstr)
    paraVale(4) = zr(jvInstr)

    do kpg = 1, npg
        xgau = zero
        ygau = zero
        zgau = zero
        do i = 1, nno
            xgau = xgau+ni(i+nno*(kpg-1))*xyz(1+3*(i-1))
            ygau = ygau+ni(i+nno*(kpg-1))*xyz(2+3*(i-1))
            zgau = zgau+ni(i+nno*(kpg-1))*xyz(3+3*(i-1))
        end do
        paraVale(1) = xgau
        paraVale(2) = ygau
        paraVale(3) = zgau

!  --   INTERPOLATION
        call fointe('FM', epsinif(1), nbPara, paraName, paraVale, epxx, paraCode)
        call fointe('FM', epsinif(2), nbPara, paraName, paraVale, epyy, paraCode)
        call fointe('FM', epsinif(3), nbPara, paraName, paraVale, epxy, paraCode)
        call fointe('FM', epsinif(4), nbPara, paraName, paraVale, kxx, paraCode)
        call fointe('FM', epsinif(5), nbPara, paraName, paraVale, kyy, paraCode)
        call fointe('FM', epsinif(6), nbPara, paraName, paraVale, kxy, paraCode)

        efge(1+8*(kpg-1)) = dm(1, 1)*epxx+dm(1, 2)*epyy+dm(1, 3)*epxy
        efge(2+8*(kpg-1)) = dm(2, 1)*epxx+dm(2, 2)*epyy+dm(2, 3)*epxy
        efge(3+8*(kpg-1)) = dm(3, 1)*epxx+dm(3, 2)*epyy+dm(3, 3)*2.d0*epxy
        efge(4+8*(kpg-1)) = df(1, 1)*kxx+df(1, 2)*kyy+df(1, 3)*kxy
        efge(5+8*(kpg-1)) = df(2, 1)*kxx+df(2, 2)*kyy+df(2, 3)*kxy
        efge(6+8*(kpg-1)) = df(3, 1)*kxx+df(3, 2)*kyy+df(3, 3)*2.d0*kxy
!
        efge(7+8*(kpg-1)) = zero
        efge(8+8*(kpg-1)) = zero
    end do
!
end subroutine
