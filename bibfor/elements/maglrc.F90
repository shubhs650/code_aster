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
! aslint: disable=W0413
subroutine maglrc(plateCara, plateOrie, jvMaterc, &
                  matr, matrElas, ecr)
!
    use plate_type
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: jvMaterc
    real(kind=8), intent(out) :: matr(50), matrElas(6, 6)
    real(kind=8), intent(inout) :: ecr(*)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropMaxi = 15
    character(len=16) :: propName(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: nbProp, i
    real(kind=8) :: vglob(3)
    real(kind=8) :: alpha, beta, vel
    character(len=16), parameter :: nonLinKeyword = ('GLRC_DAMAGE')
    character(len=16), parameter :: elasKeyword = ('ELAS_GLRC')
    real(kind=8) :: epais
!
! --------------------------------------------------------------------------------------------------
!
    matr = 0.D0
    matrElas = 0.d0

! - Check consistency of thickness
    epais = plateCara%thick
    propName(1) = 'EPAIS'
    call rcvala(jvMaterc, ' ', nonLinKeyword, &
                0, ' ', [0.d0], &
                1, propName, propVale, &
                propCode, 1)
    if (propVale(1) .ne. epais) then
        propVale(2) = epais
        call utmess('F', 'ELEMENTS5_42', nr=2, valr=propVale)
    end if

! - Elasticity (bending)
    propName(1) = 'E_F'
    propName(2) = 'NU_F'
    call rcvala(jvMaterc, ' ', elasKeyword, &
                0, ' ', [0.d0], &
                2, propName, propVale, &
                propCode, 1)
    matr(6) = propVale(1)
    matr(7) = propVale(2)

! - MATRICE ELASTIQUE MEMBRANE/CISAILLEMENT
    propName(1) = 'BN11'
    propName(2) = 'BN12'
    propName(3) = 'BN22'
    propName(4) = 'BN33'
    propName(5) = 'BT1'
    propName(6) = 'BT2'
    propName(7) = 'BM11'
    propName(8) = 'BM12'
    propName(9) = 'BM22'
    propName(10) = 'BM33'
    nbProp = 10
    call rcvala(jvMaterc, ' ', nonLinKeyword, &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    matr(1) = 1.0d0
    matr(2) = propVale(1)
    matr(3) = propVale(2)
    matr(4) = propVale(3)
    matr(5) = propVale(4)
    matrElas(4, 4) = propVale(7)
    matrElas(4, 5) = propVale(8)
    matrElas(5, 4) = matrElas(4, 5)
    matrElas(5, 5) = propVale(9)
    matrElas(6, 6) = propVale(10)
    matr(14) = propVale(5)
    matr(15) = propVale(6)

! - SEUILS ET PENTES
    propName(1) = 'MF1'
    propName(2) = 'MF2'
    propName(3) = 'QP1'
    propName(4) = 'QP2'
    propName(5) = 'GAMMA'
    nbProp = 5
    call rcvala(jvMaterc, ' ', nonLinKeyword, &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    matr(8) = propVale(1)
    matr(9) = propVale(2)
    matr(10) = propVale(3)
    matr(11) = propVale(4)
    matr(12) = propVale(5)

! - PARAMETRES TENSEUR DE PRAGER/MEMBRANE
    propName(1) = 'C1N1'
    propName(2) = 'C1N2'
    propName(3) = 'C1N3'
    propName(4) = 'C2N1'
    propName(5) = 'C2N2'
    propName(6) = 'C2N3'
    nbProp = 6
    call rcvala(jvMaterc, ' ', nonLinKeyword, &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    matr(16) = propVale(1)
    matr(17) = propVale(2)
    matr(18) = propVale(3)
    matr(22) = propVale(4)
    matr(23) = propVale(5)
    matr(24) = propVale(6)

! - PARAMETRES TENSEUR DE PRAGER/FLEXION
    propName(1) = 'C1M1'
    propName(2) = 'C1M2'
    propName(3) = 'C1M3'
    propName(4) = 'C2M1'
    propName(5) = 'C2M2'
    propName(6) = 'C2M3'
    nbProp = 6
    call rcvala(jvMaterc, ' ', nonLinKeyword, &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    matr(19) = propVale(1)
    matr(20) = propVale(2)
    matr(21) = propVale(3)
    matr(25) = propVale(4)
    matr(26) = propVale(5)
    matr(27) = propVale(6)

! - Elastic matrix
    matrElas(1, 1) = matr(2)
    matrElas(1, 2) = matr(3)
    matrElas(2, 1) = matrElas(1, 2)
    matrElas(2, 2) = matr(4)
    matrElas(3, 3) = matr(5)
!
    ASSERT(plateOrie%lRead)
    alpha = plateOrie%alpha
    beta = plateOrie%beta
!
    vglob(1) = cos(beta)*cos(alpha)
    vglob(2) = cos(beta)*sin(alpha)
    vglob(3) = -sin(beta)
    vel = vglob(1)*vglob(1)+vglob(2)*vglob(2)
    vel = vel+vglob(3)*vglob(3)
    vel = sqrt(vel)
    do i = 1, 3
        ecr(10+i) = vglob(i)/vel
    end do
!
end subroutine
