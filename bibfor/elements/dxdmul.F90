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
subroutine dxdmul(plateCara, plateOrie, &
                  lcalct, iLayer, iniv, &
                  h, d1i, d2i, x3i, hLayer)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/codent.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvala.h"
#include "asterfort/utbtab.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    aster_logical, intent(in) :: lcalct
    integer(kind=8), intent(in) :: iLayer, iniv
    real(kind=8), intent(out) :: h(3, 3), d1i(2, 2), d2i(2, 4)
    real(kind=8), intent(out) :: x3i, hLayer
!
! --------------------------------------------------------------------------------------------------
!
!     MATRICES D1I ET D2I POUR LE CALCUL DES CONTRAINTES EN MULTICOUCHE
!     (ELLES SONT FOURNIES DANS LE REPERE DE L'ELEMENT) CF BATOZ-DHATT
!
! --------------------------------------------------------------------------------------------------
!
!     IN  ICOU   : NUMERO DE LA COUCHE
!     IN  INIV   : NIVEAU DANS LA COUCHE (-1:INF , 0:MOY , 1:SUP)
!     OUT H      : MATRICE D'ELASTICITE DE LA COUCHE, REPERE INTRINSEQUE
!     OUT D1I    : MATRICE D1I REPERE INTRINSEQUE
!     OUT D2I    : MATRICE D2I REPERE INTRINSEQUE
!     OUT X3I    : Z DE CALCUL DE LA CONTRAINTE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, k, l, jLayer, jvMaterc, propCode(27)
    character(len=2) :: iValeStr
    character(len=3) :: iLayerStr
    character(len=16) :: propName(27)
    real(kind=8) :: propVale(27), r8bid
    real(kind=8) :: dx1i(2, 2), dx2i(2, 4)
    real(kind=8) :: da1i(2, 2), da2i(2, 4)
    real(kind=8) :: ordi, ai(3, 3)
    real(kind=8) :: xab1(3, 3), excen
    real(kind=8) :: t1ve(9)
!
! --------------------------------------------------------------------------------------------------
!
    excen = plateCara%offset
    t1ve = plateOrie%t1ve
!
    call jevech('PMATERC', 'L', jvMaterc)

    call codent(iLayer, 'G', iLayerStr)
    do i = 1, 27
        call codent(i, 'G', iValeStr)
        propName(i) = 'C'//iLayerStr//'_V'//iValeStr
    end do
    r8bid = 0.d0
    call rcvala(zi(jvMaterc), ' ', 'ELAS_COQMU', 0, ' ', &
                [r8bid], 9, propName, propVale, propCode, &
                1)
    hLayer = propVale(1)
    ordi = propVale(3)
    h(1, 1) = propVale(4)
    h(1, 2) = propVale(5)
    h(1, 3) = propVale(6)
    h(2, 2) = propVale(7)
    h(2, 3) = propVale(8)
    h(3, 3) = propVale(9)
    h(2, 1) = h(1, 2)
    h(3, 1) = h(1, 3)
    h(3, 2) = h(2, 3)

!     ----- CALCUL DE Z ------------------------------------------------
    x3i = ordi+excen
    if (iniv .lt. 0) then
        x3i = x3i-hLayer/2.d0
    else if (iniv .gt. 0) then
        x3i = x3i+hLayer/2.d0
    end if
!     ----- CALCUL DE D1I ET D2I ---------------------------------------
!
    dx1i = 0.d0
    dx2i = 0.d0
!
    if (lcalct) then
        do jLayer = 1, iLayer
            call codent(jLayer, 'G', iLayerStr)
            do i = 1, 27
                call codent(i, 'G', iValeStr)
                propName(i) = 'C'//iLayerStr//'_V'//iValeStr
            end do
            call rcvala(zi(jvMaterc), ' ', 'ELAS_COQMU', 0, ' ', &
                        [r8bid], 1, propName, propVale, propCode, &
                        1)
            hLayer = propVale(1)
            call rcvala(zi(jvMaterc), ' ', 'ELAS_COQMU', 0, ' ', &
                        [r8bid], 1, propName(3), propVale(3), propCode(3), &
                        1)
            ordi = propVale(3)
            call rcvala(zi(jvMaterc), ' ', 'ELAS_COQMU', 0, ' ', &
                        [r8bid], 12, propName(16), propVale(16), propCode(16), &
                        1)
!
!      RECUP MATRICE AI = H(Z).HF-1
!
            ai(1, 1) = propVale(16)
            ai(2, 1) = propVale(17)
            ai(3, 1) = propVale(18)
            ai(1, 2) = propVale(19)
            ai(2, 2) = propVale(20)
            ai(3, 2) = propVale(21)
            ai(1, 3) = propVale(22)
            ai(2, 3) = propVale(23)
            ai(3, 3) = propVale(24)
!
!      PASSAGE DANS LE REPERE INTRINSEQUE A L'ELEMENT
!
            call utbtab('ZERO', 3, 3, ai, t1ve, xab1, ai)
!
!         TERMES DE LA MATRICE INTERVENANT DANS D1(Z)
!      CF. DHAT-BATOZ VOL 2 PAGE 243
!      DAI1 MATRICE (2,2) CONSTANTE PAR COUCHE
!      TERME 1,1 : A11+A33 TERME 1,2 : A13+A32
!      TERME 2,1 : A31+A23 TERME 2,2 : A22+A33
            da1i(1, 1) = ai(1, 1)+ai(3, 3)
            da1i(1, 2) = ai(1, 3)+ai(3, 2)
            da1i(2, 1) = ai(3, 1)+ai(2, 3)
            da1i(2, 2) = ai(2, 2)+ai(3, 3)
!         TERMES DE LA MATRICE INTERVENANT DANS D2(Z)
!      CF. DHAT-BATOZ VOL 2 PAGE 243
!      DAI2 MATRICE (2,4) CONSTANTE PAR COUCHE
            da2i(1, 1) = ai(1, 1)-ai(3, 3)
            da2i(1, 2) = ai(1, 3)-ai(3, 2)
            da2i(1, 3) = ai(1, 2)*2.d0
            da2i(1, 4) = ai(3, 1)*2.d0
            da2i(2, 1) = ai(3, 1)-ai(2, 3)
            da2i(2, 2) = ai(3, 3)-ai(2, 2)
            da2i(2, 3) = ai(3, 2)*2.d0
            da2i(2, 4) = ai(2, 1)*2.d0
!
!
!      D1(Z)=SOMME(-T,Z)(-Z/2*DA1I DZ)
!      TOUS CALCULS FAITS : K INDICE MAX TEL QUE ZK < X3I
!      D1=SOMME(I=1,K-1)(-1/2*DA1I*EPI*ORDI)+DAIK(ZK**2-X3I**2)
!
            do k = 1, 2
                do l = 1, 2
                    d1i(k, l) = dx1i(k, l)+((ordi-hLayer/2.d0)**2-x3i*x3i)*da1i(k, l)/4.d0
                    dx1i(k, l) = dx1i(k, l)-hLayer*ordi*da1i(k, l)/2.d0
                end do
            end do
!
!      D2(Z)=SOMME(-T,Z)(-Z/2*DA2I DZ)
!      TOUS CALCULS FAITS : K INDICE MAX TEL QUE ZK < X3I
!      D2=SOMME(I=1,K-1)(-1/2*DA2I*EPI*ORDI)+DA2K(ZK**2-X3I**2)
!
            do k = 1, 2
                do l = 1, 4
                    d2i(k, l) = dx2i(k, l)+((ordi-hLayer/2.d0)**2-x3i*x3i)*da2i(k, l)/4.d0
                    dx2i(k, l) = dx2i(k, l)-hLayer*ordi*da2i(k, l)/2.d0
                end do
            end do
        end do
!
    end if
!
!     MATRICE H DANS LE REPERE INTRINSEQUE DE L'ELEMENT
!
    call utbtab('ZERO', 3, 3, h, t1ve, xab1, h)
!
end subroutine
