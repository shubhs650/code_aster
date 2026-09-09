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
! aslint: disable=W1501
!
subroutine elraga(elrefz, fapz, ndim, nbpg, coopg, poipg)
!
    implicit none
!
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/elraca.h"
#include "MeshTypes_type.h"
!
    character(len=*), intent(in) :: elrefz, fapz
    integer(kind=8), intent(out) :: nbpg, ndim
    real(kind=8), intent(out) :: coopg(*), poipg(*)
!
! --------------------------------------------------------------------------------------------------
!
! Finite elements management
!
! Get parameters of integration scheme
!
! --------------------------------------------------------------------------------------------------
!
! In  elrefe           : name of geometric support for finite element
! In  fapg             : name of Gauss integration scheme
! Out ndim             : topological dimension (0/1/2/3)
! Out nbpg             : number of points of integration schemes
! Out coopg            : coordinataes of Gauss points
! Out poipg            : weight of Gauss points
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elrefa, fapg, nofpg(MT_NBFAMX)
    integer(kind=8) :: i, npar, npi, ix, iy, iz, npx, npyz
    integer(kind=8) :: nno, nnos, nbfpg, nbpg1(MT_NBFAMX), iNode, ifam
    real(kind=8) :: xpg(MT_NBPGMX), ypg(MT_NBPGMX), zpg(MT_NBPGMX), hpg(MT_NBPGMX)
    real(kind=8) :: h(8), a(8)
    real(kind=8) :: aty(40), ht(40), atz(40)
    real(kind=8) :: lobWeight(7), lobCoor(7), quadWeight(3), quadCoor(3)
    real(kind=8) :: a1, a2, b1, b2, c1, c2, d1, e1
    real(kind=8) :: h1, h2, h3, h5
    real(kind=8) :: p1, p2, p3, p4, p5
    real(kind=8) :: xno(3*MT_NNOMAX), vol
    real(kind=8), parameter :: zero = 0.d0, undemi = 0.5d0
    real(kind=8), parameter :: un = 1.d0, deux = 2.d0

    real(kind=8), parameter :: rac5 = sqrt(5.d0), rac15 = sqrt(15.d0), rac30 = sqrt(30.d0)
    real(kind=8), parameter :: rac_1div3 = sqrt(1.d0/3.d0), rac_3div5 = sqrt(3.d0/5.d0)
    real(kind=8), parameter :: rac_3div7 = sqrt(3.d0/7.d0)

    real(kind=8), parameter :: gauss4p12 = sqrt(3.d0/7.d0-2.d0/7.d0*sqrt(6.d0/5.d0))
    real(kind=8), parameter :: gauss4p34 = sqrt(3.d0/7.d0+2.d0/7.d0*sqrt(6.d0/5.d0))

    real(kind=8), parameter :: lobatto7p35 = sqrt(5.d0/11.d0-2.d0/11.d0*sqrt(5.d0/3.d0))
    real(kind=8), parameter :: lobatto7p26 = sqrt(5.d0/11.d0+2.d0/11.d0*sqrt(5.d0/3.d0))

!
#define t(u) 2.0d0*(u) - 1.0d0
!
! --------------------------------------------------------------------------------------------------
!
    elrefa = elrefz
    fapg = fapz

! - Get list of integration schemes of geometric support
    call elraca(elrefa, &
                nbfpg, nofpg, nbpg1, &
                ndim, nno, nnos, &
                xno, vol)
    ASSERT((ndim .ge. 0) .and. (ndim .le. 3))

! - Get index for integration scheme
    ifam = indik8(nofpg, fapg, 1, nbfpg)
    ASSERT(ifam .gt. 0)

! - Get number of Gauss points
    nbpg = nbpg1(ifam)

! - For 'NOEU' scheme
    if (fapg .eq. 'NOEU') then
        ASSERT(nbpg .eq. nno)
        do iNode = 1, nno
            hpg(iNode) = vol/nno
            if (ndim .ge. 1) xpg(iNode) = xno(ndim*(iNode-1)+1)
            if (ndim .ge. 2) ypg(iNode) = xno(ndim*(iNode-1)+2)
            if (ndim .eq. 3) zpg(iNode) = xno(ndim*(iNode-1)+3)
        end do
        goto 170
    end if

! - For 'NOEU_S' scheme
    if (fapg .eq. 'NOEU_S') then
        ASSERT(nbpg .eq. nnos)
        do iNode = 1, nnos
            hpg(iNode) = vol/nnos
            if (ndim .ge. 1) xpg(iNode) = xno(ndim*(iNode-1)+1)
            if (ndim .ge. 2) ypg(iNode) = xno(ndim*(iNode-1)+2)
            if (ndim .eq. 3) zpg(iNode) = xno(ndim*(iNode-1)+3)
        end do
        goto 170
    end if

! - For 'FPG1' scheme
    if (fapg .eq. 'FPG1') then
        ASSERT(nbpg .eq. 1)
        xpg(1) = zero
        if (ndim .ge. 1) xpg(1) = zero
        if (ndim .ge. 2) ypg(1) = zero
        if (ndim .eq. 3) zpg(1) = zero
        do iNode = 1, nno
            if (ndim .ge. 1) xpg(1) = xpg(1)+xno(ndim*(iNode-1)+1)
            if (ndim .ge. 2) ypg(1) = ypg(1)+xno(ndim*(iNode-1)+2)
            if (ndim .eq. 3) zpg(1) = zpg(1)+xno(ndim*(iNode-1)+3)
        end do
        if (ndim .ge. 1) xpg(1) = xpg(1)/nno
        if (ndim .ge. 2) ypg(1) = ypg(1)/nno
        if (ndim .eq. 3) zpg(1) = zpg(1)/nno
        hpg(1) = vol
        goto 170
    end if

! - For other schemes
    if (elrefa .eq. 'HE8' .or. elrefa .eq. 'H20' .or. elrefa .eq. 'H27') then
        npar = 0
        if (fapg .eq. 'FPG1') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 1 POINTS ( ORDRE 1 )
            xpg(1) = zero
            ypg(1) = zero
            zpg(1) = zero
            hpg(1) = 8.d0

        else if (fapg .eq. 'FPG8') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 2 POINTS DANS CHAQUE DIRECTION ( ORDRE 3 )
            npar = 2
            a(1) = -rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un

        else if (fapg .eq. 'FPG27') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 3 POINTS DANS CHAQUE DIRECTION ( ORDRE 5 )
            npar = 3
            a(1) = -rac_3div5
            a(2) = zero
            a(3) = -a(1)
            h(1) = 5.d0/9.d0
            h(2) = 8.d0/9.d0
            h(3) = h(1)

        else if (fapg .eq. 'FPG64') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 4 POINTS DANS CHAQUE DIRECTION ( ORDRE 7 )
            npar = 4
            a(1) = -gauss4p12
            a(2) = -a(1)
            a(3) = -gauss4p34
            a(4) = -a(3)
            h(1) = (18.d0+rac30)/36.d0
            h(2) = h(1)
            h(3) = (18.d0-rac30)/36.d0
            h(4) = h(3)

        else if (fapg .eq. 'FPG125') then
            ! order 9 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            a(1) = 0.d0
            a(2) = 0.53846931010568309103d0
            a(3) = -a(2)
            a(4) = 0.90617984593866399279d0
            a(5) = -a(4)
            h(1) = 0.56888888888888888888d0
            h(2) = 0.47862867049936646804d0
            h(3) = h(2)
            h(4) = 0.23692688505618908751d0
            h(5) = h(4)
            npar = 5

        else if (fapg .eq. 'FPG216') then
            ! order 11 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            a(1) = 0.661209386466264513661d0
            a(2) = -a(1)
            a(3) = 0.238619186083196908630d0
            a(4) = -a(3)
            a(5) = 0.932469514203152027812d0
            a(6) = -a(5)
            h(1) = 0.360761573048138607569d0
            h(2) = h(1)
            h(3) = 0.467913934572691047389d0
            h(4) = h(3)
            h(5) = 0.171324492379170345040d0
            h(6) = h(5)
            npar = 6

        else if (fapg .eq. 'FPG343') then
            ! order 13 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            a(1) = 0.d0
            a(2) = 0.405845151377397166906d0
            a(3) = -a(2)
            a(4) = 0.741531185599394439863d0
            a(5) = -a(4)
            a(6) = 0.949107912342758524526d0
            a(7) = -a(6)
            h(1) = 0.417959183673469387755d0
            h(2) = 0.381830050505118944950d0
            h(3) = h(2)
            h(4) = 0.279705391489276667901d0
            h(5) = h(4)
            h(6) = 0.129484966168869693270d0
            h(7) = h(6)
            npar = 7

        else if (fapg .eq. 'FPG512') then
            ! order 15 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            a(1) = 0.183434642495649804939d0
            a(2) = -a(1)
            a(3) = 0.525532409916328985817d0
            a(4) = -a(3)
            a(5) = 0.796666477413626739591d0
            a(6) = -a(5)
            a(7) = 0.960289856497536231683d0
            a(8) = -a(7)
            h(1) = 0.362683783378361982965d0
            h(2) = h(1)
            h(3) = 0.313706645877887287337d0
            h(4) = h(3)
            h(5) = 0.222381034453374470544d0
            h(6) = h(5)
            h(7) = 0.101228536290376259152d0
            h(8) = h(7)
            npar = 8

        else if (fapg .eq. 'FPG8NOS') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 2 POINTS DANS CHAQUE DIRECTION ( ORDRE 3 )
            npar = 2
            a(1) = -rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+8) = vol/nnos
                xpg(iNode+8) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+8) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+8) = xno(ndim*(iNode-1)+3)
            end do

        else
            ASSERT(ASTER_FALSE)

        end if
        npi = 0
        do ix = 1, npar
            do iy = 1, npar
                do iz = 1, npar
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = a(iy)
                    zpg(npi) = a(iz)
                    hpg(npi) = h(ix)*h(iy)*h(iz)
                end do
            end do
        end do

    else if (elrefa .eq. 'HE9') then
        if (fapg .eq. 'LOB5') then
! --------- FORMULE DE QUADRATURE DE GAUSS-LOBATTO A 5 POINTS DANS
!           L EPAISSEUR AU CENTRE DE L'ELEMENT
            lobCoor(1) = -un
            lobCoor(2) = -rac_3div7
            lobCoor(3) = zero
            lobCoor(4) = -lobCoor(2)
            lobCoor(5) = -lobCoor(1)
            lobWeight(1) = 0.1d0
            lobWeight(2) = 49.d0/90.d0
            lobWeight(3) = 32.0/45.d0
            lobWeight(4) = lobWeight(2)
            lobWeight(5) = lobWeight(1)
            do iz = 1, 5
                xpg(iz) = zero
                ypg(iz) = zero
                zpg(iz) = lobCoor(iz)
                hpg(iz) = lobWeight(iz)*4.d0
            end do

        else if (fapg .eq. 'LOB7') then
! --------- FORMULE DE QUADRATURE DE GAUSS-LOBATTO A 7 POINTS DANS
!           L EPAISSEUR AU CENTRE DE L'ELEMENT
            lobCoor(1) = -un
            lobCoor(2) = -lobatto7p26
            lobCoor(3) = -lobatto7p35
            lobCoor(4) = zero
            lobCoor(5) = -lobCoor(3)
            lobCoor(6) = -lobCoor(2)
            lobCoor(7) = -lobCoor(1)
            lobWeight(1) = un/21.d0
            lobWeight(2) = (124.d0-7.d0*rac15)/350.d0
            lobWeight(3) = (124.d0+7.d0*rac15)/350.d0
            lobWeight(4) = 256.d0/525.d0
            lobWeight(5) = lobWeight(3)
            lobWeight(6) = lobWeight(2)
            lobWeight(7) = lobWeight(1)
            do iz = 1, 7
                xpg(iz) = zero
                ypg(iz) = zero
                zpg(iz) = lobCoor(iz)
                hpg(iz) = lobWeight(iz)*4.d0
            end do

        else if (fapg .eq. 'FPG8') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 2 POINTS DANS CHAQUE DIRECTION ( ORDRE 3 )
            npar = 2
            a(1) = -rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un
            npi = 0
            do ix = 1, npar
                do iy = 1, npar
                    do iz = 1, npar
                        npi = npi+1
                        xpg(npi) = a(ix)
                        ypg(npi) = a(iy)
                        zpg(npi) = a(iz)
                        hpg(npi) = h(ix)*h(iy)*h(iz)
                    end do
                end do
            end do

        else
            ASSERT(ASTER_FALSE)
        end if
!
    else if (elrefa .eq. 'H10') then
        if (fapg .eq. 'FIS2') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 2 POINTS DANS
!           LA DIRECTION DU PREMIER AXE ET AU CENTRE DE L'ELEMENT
            ! order 3
            xpg(1) = -rac_1div3
            xpg(2) = -xpg(1)
            hpg(1) = un
            hpg(2) = hpg(1)
            do iz = 1, 2
                ypg(iz) = zero
                zpg(iz) = zero
            end do

        else if (fapg .eq. 'FIS3') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 3 POINTS DANS
!           LA DIRECTION DU PREMIER AXE ET AU CENTRE DE L'ELEMENT
            ! order 5
            quadCoor(1) = -rac_3div5
            quadCoor(2) = zero
            quadCoor(3) = -quadCoor(1)
            quadWeight(1) = 5.d0/9.d0
            quadWeight(2) = 8.d0/9.d0
            quadWeight(3) = quadWeight(1)
            do iz = 1, 3
                xpg(iz) = quadCoor(iz)
                ypg(iz) = zero
                zpg(iz) = zero
                hpg(iz) = quadWeight(iz)
            end do

        else if (fapg .eq. 'FIS4') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 4 POINTS DANS
!           LA DIRECTION DU PREMIER AXE ET AU CENTRE DE L'ELEMENT
            ! order 7
            xpg(1) = gauss4p12
            xpg(2) = -xpg(1)
            xpg(3) = gauss4p34
            xpg(4) = -xpg(3)
            hpg(1) = (18.d0+rac30)/36.d0
            hpg(2) = hpg(1)
            hpg(3) = (18.d0-rac30)/36.d0
            hpg(4) = hpg(3)
            do iz = 1, 4
                ypg(iz) = zero
                zpg(iz) = zero
            end do

        else
            ASSERT(ASTER_FALSE)
        end if
!
    else if (elrefa .eq. 'PE6' .or. elrefa .eq. 'P15' .or. elrefa .eq. 'P18' .or. &
             elrefa .eq. 'P21') then
        if (fapg .eq. 'FPG6') then
! --------- FORMULE A 4 * 2 POINTS (CF TOUZOT PAGE 297) -> ORDRE 2
! --------- FORMULE DE GAUSS - 2 POINTS DE GAUSS  EN X (ORDRE 3)
            npx = 2
            a(1) = -rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un

! --------- FORMULE DE HAMMER - 3 POINTS DE HAMMER EN Y Z (ORDRE 2 EN Y Z)
            npyz = 3
            aty(1) = undemi
            aty(2) = zero
            aty(3) = undemi
            atz(1) = undemi
            atz(2) = undemi
            atz(3) = zero

            ht(1) = un/6.d0
            ht(2) = ht(1)
            ht(3) = ht(1)

        else if (fapg .eq. 'FPG6B') then
! --------- FORMULE A 4 * 2 POINTS (CF TOUZOT PAGE 297) -> ORDRE 2
! --------- FORMULE DE GAUSS - 2 POINTS DE GAUSS  EN X (ORDRE 3)
            npx = 2
            a(1) = -rac_1div3
            a(2) = -a(1)

            h(1) = un
            h(2) = un

! --------- FORMULE DE HAMMER - 3 POINTS DE HAMMER EN Y Z (ORDRE 2 EN Y Z)
            npyz = 3
            aty(1) = un/6.d0
            aty(2) = deux/3.d0
            aty(3) = un/6.d0
            atz(1) = un/6.d0
            atz(2) = un/6.d0
            atz(3) = deux/3.d0

            ht(1) = un/6.d0
            ht(2) = ht(1)
            ht(3) = ht(1)

        else if (fapg .eq. 'FPG8') then
! --------- FORMULE A 4 * 2 POINTS (CF TOUZOT PAGE 297) -> ORDRE 3
! --------- FORMULE DE GAUSS - 2 POINTS DE GAUSS  EN X (ORDRE 3)
            npx = 2
            a(1) = -rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un

! --------- FORMULE DE HAMMER - 4 POINTS DE HAMMER EN Y Z (ORDRE 3 EN Y Z)
            npyz = 4
            aty(1) = un/3.d0
            aty(2) = 0.6d0
            aty(3) = 0.2d0
            aty(4) = 0.2d0
            atz(1) = un/3.d0
            atz(2) = 0.2d0
            atz(3) = 0.6d0
            atz(4) = 0.2d0
            ht(1) = -27.d0/96.d0
            ht(2) = 25.d0/96.d0
            ht(3) = ht(2)
            ht(4) = ht(2)

        else if (fapg .eq. 'FPG21') then
! --------- FORMULE A 7 * 3 POINTS :   (CF TOUZOT PAGE 298)
! --------- FORMULE DE GAUSS - 3 POINTS DE GAUSS EN X (ORDRE 5)
            npx = 3
            a(1) = -rac_3div5
            a(2) = zero
            a(3) = -a(1)
            h(1) = 5.d0/9.d0
            h(2) = 8.d0/9.d0
            h(3) = h(1)

! --------- FORMULE DE HAMMER - 7 POINTS DE HAMMER EN Y Z (ORDRE 5 EN Y Z)
            npyz = 7
            aty(1) = un/3.d0
            atz(1) = un/3.d0
            aty(2) = (6.d0+rac15)/21.d0
            atz(2) = aty(2)
            aty(3) = un-deux*aty(2)
            atz(3) = aty(2)
            aty(4) = aty(2)
            atz(4) = un-deux*aty(2)
            aty(5) = (6.d0-rac15)/21.d0
            atz(5) = aty(5)
            aty(6) = un-deux*aty(5)
            atz(6) = aty(5)
            aty(7) = aty(5)
            atz(7) = un-deux*aty(5)
            ht(1) = 9.d0/80.d0
            ht(2) = (155.d0+rac15)/2400.d0
            ht(3) = ht(2)
            ht(4) = ht(2)
            ht(5) = (155.d0-rac15)/2400.d0
            ht(6) = ht(5)
            ht(7) = ht(5)

        else if (fapg .eq. 'FPG29') then
            ! ORDRE 6
            ! x = zeta
            ! y = (xi + 1) / 2
            ! z = (eta + 1) / 2
            ! h = weight / 4

            hpg(1) = 0.027191062410231d0
            hpg(2) = 0.027191062410231d0
            hpg(3) = 0.027191062410231d0
            hpg(4) = 0.040636041641220d0
            hpg(5) = 0.040636041641220d0
            hpg(6) = 0.040636041641220d0
            hpg(7) = 0.040636041641220d0
            hpg(8) = 0.040636041641220d0
            hpg(9) = 0.040636041641220d0
            hpg(10) = 0.050275140937507d0
            hpg(11) = 0.011774414962347d0
            hpg(12) = 0.011774414962347d0
            hpg(13) = 0.011774414962347d0
            hpg(14) = 0.041951149272741d0
            hpg(15) = 0.041951149272741d0
            hpg(16) = 0.041951149272741d0
            hpg(17) = 0.041951149272741d0
            hpg(18) = 0.041951149272741d0
            hpg(19) = 0.041951149272741d0
            hpg(20) = 0.050275140937507d0
            hpg(21) = 0.011774414962347d0
            hpg(22) = 0.011774414962347d0
            hpg(23) = 0.011774414962347d0
            hpg(24) = 0.041951149272741d0
            hpg(25) = 0.041951149272741d0
            hpg(26) = 0.041951149272741d0
            hpg(27) = 0.041951149272741d0
            hpg(28) = 0.041951149272741d0
            hpg(29) = 0.041951149272741d0

            xpg(1) = 0.000000000000000d0
            xpg(2) = 0.000000000000000d0
            xpg(3) = 0.000000000000000d0
            xpg(4) = 0.000000000000000d0
            xpg(5) = 0.000000000000000d0
            xpg(6) = 0.000000000000000d0
            xpg(7) = 0.000000000000000d0
            xpg(8) = 0.000000000000000d0
            xpg(9) = 0.000000000000000d0
            xpg(10) = 0.936241512371697d0
            xpg(11) = 0.948681147283254d0
            xpg(12) = 0.948681147283254d0
            xpg(13) = 0.948681147283254d0
            xpg(14) = 0.600638052820557d0
            xpg(15) = 0.600638052820557d0
            xpg(16) = 0.600638052820557d0
            xpg(17) = 0.600638052820557d0
            xpg(18) = 0.600638052820557d0
            xpg(19) = 0.600638052820557d0
            xpg(20) = -0.936241512371697d0
            xpg(21) = -0.948681147283254d0
            xpg(22) = -0.948681147283254d0
            xpg(23) = -0.948681147283254d0
            xpg(24) = -0.600638052820557d0
            xpg(25) = -0.600638052820557d0
            xpg(26) = -0.600638052820557d0
            xpg(27) = -0.600638052820557d0
            xpg(28) = -0.600638052820557d0
            xpg(29) = -0.600638052820557d0

            ypg(1) = 0.895512822481133d0
            ypg(2) = 0.052243588759434d0
            ypg(3) = 0.052243588759434d0
            ypg(4) = 0.198304865473555d0
            ypg(5) = 0.198304865473555d0
            ypg(6) = 0.270635256143164d0
            ypg(7) = 0.531059878383280d0
            ypg(8) = 0.531059878383280d0
            ypg(9) = 0.270635256143164d0
            ypg(10) = 0.333333333333333d0
            ypg(11) = 0.841699897299232d0
            ypg(12) = 0.079150051350384d0
            ypg(13) = 0.079150051350384d0
            ypg(14) = 0.054831294873304d0
            ypg(15) = 0.054831294873304d0
            ypg(16) = 0.308513201856883d0
            ypg(17) = 0.636655503269814d0
            ypg(18) = 0.636655503269814d0
            ypg(19) = 0.308513201856883d0
            ypg(20) = 0.333333333333333d0
            ypg(21) = 0.841699897299232d0
            ypg(22) = 0.079150051350384d0
            ypg(23) = 0.079150051350384d0
            ypg(24) = 0.054831294873304d0
            ypg(25) = 0.054831294873304d0
            ypg(26) = 0.308513201856883d0
            ypg(27) = 0.636655503269814d0
            ypg(28) = 0.636655503269814d0
            ypg(29) = 0.308513201856883d0

            zpg(1) = 0.052243588759434d0
            zpg(2) = 0.895512822481133d0
            zpg(3) = 0.052243588759434d0
            zpg(4) = 0.270635256143164d0
            zpg(5) = 0.531059878383280d0
            zpg(6) = 0.531059878383280d0
            zpg(7) = 0.270635256143164d0
            zpg(8) = 0.198304865473555d0
            zpg(9) = 0.198304865473555d0
            zpg(10) = 0.333333333333333d0
            zpg(11) = 0.079150051350384d0
            zpg(12) = 0.841699897299232d0
            zpg(13) = 0.079150051350384d0
            zpg(14) = 0.308513201856883d0
            zpg(15) = 0.636655503269814d0
            zpg(16) = 0.636655503269814d0
            zpg(17) = 0.308513201856883d0
            zpg(18) = 0.054831294873304d0
            zpg(19) = 0.054831294873304d0
            zpg(20) = 0.333333333333333d0
            zpg(21) = 0.079150051350384d0
            zpg(22) = 0.841699897299232d0
            zpg(23) = 0.079150051350384d0
            zpg(24) = 0.308513201856883d0
            zpg(25) = 0.636655503269814d0
            zpg(26) = 0.636655503269814d0
            zpg(27) = 0.308513201856883d0
            zpg(28) = 0.054831294873304d0
            zpg(29) = 0.054831294873304d0
!
        else if (fapg .eq. 'FPG52') then
            ! order 7
            npx = 4
            a(1) = gauss4p12
            a(2) = -a(1)
            a(3) = gauss4p34
            a(4) = -a(3)
            h(1) = (18.d0+rac30)/36.d0
            h(2) = h(1)
            h(3) = (18.d0-rac30)/36.d0
            h(4) = h(3)

!         FORMULE A 13 POINTS : ORDRE 7  (CF BATHE, PAGE 280)
            npyz = 13
            aty(1) = 0.0651301029022d0
            atz(1) = 0.0651301029022d0
            aty(2) = 0.8697397941956d0
            atz(2) = 0.0651301029022d0
            aty(3) = 0.0651301029022d0
            atz(3) = 0.8697397941956d0
            aty(4) = 0.3128654960049d0
            atz(4) = 0.0486903154253d0
            aty(5) = 0.6384441885698d0
            atz(5) = 0.3128654960049d0
            aty(6) = 0.0486903154253d0
            atz(6) = 0.6384441885698d0
            aty(7) = 0.6384441885698d0
            atz(7) = 0.0486903154253d0
            aty(8) = 0.3128654960049d0
            atz(8) = 0.6384441885698d0
            aty(9) = 0.0486903154253d0
            atz(9) = 0.3128654960049d0
            aty(10) = 0.2603459660790d0
            atz(10) = 0.2603459660790d0
            aty(11) = 0.4793080678419d0
            atz(11) = 0.2603459660790d0
            aty(12) = 0.2603459660790d0
            atz(12) = 0.4793080678419d0
            aty(13) = 0.3333333333333d0
            atz(13) = 0.3333333333333d0
            p1 = 0.0533472356088d0/deux
            p2 = 0.0771137608903d0/deux
            p3 = 0.1756152574332d0/deux
            p4 = -0.1495700444677d0/deux
            ht(1) = p1
            ht(2) = p1
            ht(3) = p1
            ht(4) = p2
            ht(5) = p2
            ht(6) = p2
            ht(7) = p2
            ht(8) = p2
            ht(9) = p2
            ht(10) = p3
            ht(11) = p3
            ht(12) = p3
            ht(13) = p4
!
        else if (fapg .eq. 'FPG95') then
            ! order 9 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            npx = 5
            a(1) = 0.d0
            a(2) = 0.53846931010568309103d0
            a(3) = -a(2)
            a(4) = 0.90617984593866399279d0
            a(5) = -a(4)
            h(1) = 0.56888888888888888888d0
            h(2) = 0.47862867049936646804d0
            h(3) = h(2)
            h(4) = 0.23692688505618908751d0
            h(5) = h(4)

! --------- Order 9 : Very high-order symmetric positive-interior
!                      quadrature rules on triangles and tetrahedra
!                      Zelalem Arega Worku  and Jason E. Hicken
!           https://github.com/OptimalDesignLab/SummationByParts.jl
            npyz = 19
            aty(1) = 0.0447295133944527098d0
            atz(1) = 0.9105409732110945803d0
            ht(1) = 0.0127888378293490156d0
            aty(2) = 0.0447295133944527098d0
            atz(2) = 0.0447295133944527098d0
            ht(2) = 0.0127888378293490156d0
            aty(3) = 0.9105409732110945803d0
            atz(3) = 0.0447295133944527098d0
            ht(3) = 0.0127888378293490156d0
            aty(4) = 0.1882035356190327302d0
            atz(4) = 0.6235929287619345395d0
            ht(4) = 0.0398238694636051265d0
            aty(5) = 0.1882035356190327302d0
            atz(5) = 0.1882035356190327302d0
            ht(5) = 0.0398238694636051265d0
            aty(6) = 0.6235929287619345395d0
            atz(6) = 0.1882035356190327302d0
            ht(6) = 0.0398238694636051265d0
            aty(7) = 0.4896825191987376277d0
            atz(7) = 0.0206349616025247444d0
            ht(7) = 0.0156673501135695352d0
            aty(8) = 0.4896825191987376277d0
            atz(8) = 0.4896825191987376277d0
            ht(8) = 0.0156673501135695352d0
            aty(9) = 0.0206349616025247444d0
            atz(9) = 0.4896825191987376277d0
            ht(9) = 0.0156673501135695352d0
            aty(10) = 0.4370895914929366372d0
            atz(10) = 0.1258208170141267254d0
            ht(10) = 0.0389137705023871396d0
            aty(11) = 0.4370895914929366372d0
            atz(11) = 0.4370895914929366372d0
            ht(11) = 0.0389137705023871396d0
            aty(12) = 0.1258208170141267254d0
            atz(12) = 0.4370895914929366372d0
            ht(12) = 0.0389137705023871396d0
            aty(13) = 0.0368384120547362836d0
            atz(13) = 0.2219629891607656956d0
            ht(13) = 0.0216417696886446886d0
            aty(14) = 0.7411985987844980207d0
            atz(14) = 0.2219629891607656956d0
            ht(14) = 0.0216417696886446886d0
            aty(15) = 0.7411985987844980207d0
            atz(15) = 0.0368384120547362836d0
            ht(15) = 0.0216417696886446886d0
            aty(16) = 0.0368384120547362836d0
            atz(16) = 0.7411985987844980207d0
            ht(16) = 0.0216417696886446886d0
            aty(17) = 0.2219629891607656956d0
            atz(17) = 0.7411985987844980207d0
            ht(17) = 0.0216417696886446886d0
            aty(18) = 0.2219629891607656956d0
            atz(18) = 0.0368384120547362836d0
            ht(18) = 0.0216417696886446886d0
            aty(19) = 0.3333333333333333333d0
            atz(19) = 0.3333333333333333333d0
            ht(19) = 0.0485678981413994169d0
!
        else if (fapg .eq. 'FPG168') then
            ! order 11 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            npx = 6
            a(1) = 0.661209386466264513661d0
            a(2) = -a(1)
            a(3) = 0.238619186083196908630d0
            a(4) = -a(3)
            a(5) = 0.932469514203152027812d0
            a(6) = -a(5)
            h(1) = 0.360761573048138607569d0
            h(2) = h(1)
            h(3) = 0.467913934572691047389d0
            h(4) = h(3)
            h(5) = 0.171324492379170345040d0
            h(6) = h(5)

! --------- Order 11 : Very high-order symmetric positive-interior
!                      quadrature rules on triangles and tetrahedra
!                      Zelalem Arega Worku  and Jason E. Hicken
!           https://github.com/OptimalDesignLab/SummationByParts.jl
            npyz = 28
            aty(1) = 0.0306255245353121603d0
            atz(1) = 0.9387489509293756793d0
            ht(1) = 0.0060348721385088178d0
            aty(2) = 0.0306255245353121603d0
            atz(2) = 0.0306255245353121603d0
            ht(2) = 0.0060348721385088178d0
            aty(3) = 0.9387489509293756793d0
            atz(3) = 0.0306255245353121603d0
            ht(3) = 0.0060348721385088178d0
            aty(4) = 0.1120794384974225991d0
            atz(4) = 0.7758411230051548016d0
            ht(4) = 0.0199610731868902660d0
            aty(5) = 0.1120794384974225991d0
            atz(5) = 0.1120794384974225991d0
            ht(5) = 0.0199610731868902660d0
            aty(6) = 0.7758411230051548016d0
            atz(6) = 0.1120794384974225991d0
            ht(6) = 0.0199610731868902660d0
            aty(7) = 0.2139994793245382346d0
            atz(7) = 0.5720010413509235306d0
            ht(7) = 0.0340877415313096921d0
            aty(8) = 0.2139994793245382346d0
            atz(8) = 0.2139994793245382346d0
            ht(8) = 0.0340877415313096921d0
            aty(9) = 0.5720010413509235306d0
            atz(9) = 0.2139994793245382346d0
            ht(9) = 0.0340877415313096921d0
            aty(10) = 0.4983342186162269636d0
            atz(10) = 0.0033315627675460727d0
            ht(10) = 0.0065100046990440848d0
            aty(11) = 0.4983342186162269636d0
            atz(11) = 0.4983342186162269636d0
            ht(11) = 0.0065100046990440848d0
            aty(12) = 0.0033315627675460727d0
            atz(12) = 0.4983342186162269636d0
            ht(12) = 0.0065100046990440848d0
            aty(13) = 0.4369864234941642685d0
            atz(13) = 0.1260271530116714628d0
            ht(13) = 0.0318464179605561488d0
            aty(14) = 0.4369864234941642685d0
            atz(14) = 0.4369864234941642685d0
            ht(14) = 0.0318464179605561488d0
            aty(15) = 0.1260271530116714628d0
            atz(15) = 0.4369864234941642685d0
            ht(15) = 0.0318464179605561488d0
            aty(16) = 0.0137073703144710521d0
            atz(16) = 0.1582595826844803475d0
            ht(16) = 0.0070380771122731376d0
            aty(17) = 0.8280330470010486002d0
            atz(17) = 0.1582595826844803475d0
            ht(17) = 0.0070380771122731376d0
            aty(18) = 0.8280330470010486002d0
            atz(18) = 0.0137073703144710521d0
            ht(18) = 0.0070380771122731376d0
            aty(19) = 0.0137073703144710521d0
            atz(19) = 0.8280330470010486002d0
            ht(19) = 0.0070380771122731376d0
            aty(20) = 0.1582595826844803475d0
            atz(20) = 0.8280330470010486002d0
            ht(20) = 0.0070380771122731376d0
            aty(21) = 0.1582595826844803475d0
            atz(21) = 0.0137073703144710521d0
            ht(21) = 0.0070380771122731376d0
            aty(22) = 0.0474766004458184237d0
            atz(22) = 0.3079833895041687134d0
            ht(22) = 0.0202382836943673032d0
            aty(23) = 0.6445400100500128628d0
            atz(23) = 0.3079833895041687134d0
            ht(23) = 0.0202382836943673032d0
            aty(24) = 0.6445400100500128628d0
            atz(24) = 0.0474766004458184237d0
            ht(24) = 0.0202382836943673032d0
            aty(25) = 0.0474766004458184237d0
            atz(25) = 0.6445400100500128628d0
            ht(25) = 0.0202382836943673032d0
            aty(26) = 0.3079833895041687134d0
            atz(26) = 0.6445400100500128628d0
            ht(26) = 0.0202382836943673032d0
            aty(27) = 0.3079833895041687134d0
            atz(27) = 0.0474766004458184237d0
            ht(27) = 0.0202382836943673032d0
            aty(28) = 0.3333333333333333333d0
            atz(28) = 0.3333333333333333333d0
            ht(28) = 0.0410215066112303255d0

        else if (fapg .eq. 'FPG6NOS') then
! --------- POUR LES POINTS DE GAUSS
            npx = 2
            npyz = 3
            a(1) = -rac_1div3
            a(2) = rac_1div3
            aty(1) = undemi
            aty(2) = zero
            aty(3) = undemi
            atz(1) = undemi
            atz(2) = undemi
            atz(3) = zero
            h(1) = un
            h(2) = un
            ht(1) = un/6.d0
            ht(2) = ht(1)
            ht(3) = ht(1)
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+6) = vol/nnos
                xpg(iNode+6) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+6) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+6) = xno(ndim*(iNode-1)+3)
            end do
        else
            ASSERT(ASTER_FALSE)

        end if

        if (fapg .ne. 'FPG29') then
            npi = 0
            do ix = 1, npx
                do iy = 1, npyz
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = aty(iy)
                    zpg(npi) = atz(iy)
                    hpg(npi) = h(ix)*ht(iy)
                end do
            end do
        end if

    else if (elrefa .eq. 'PE7') then
        if (fapg .eq. 'LOB5') then
! --------- FORMULE DE QUADRATURE DE GAUSS-LOBATTO A 5 POINTS DANS
!           L EPAISSEUR AU CENTRE DE L'ELEMENT
            lobCoor(1) = -un
            lobCoor(2) = -rac_3div7
            lobCoor(3) = zero
            lobCoor(4) = -lobCoor(2)
            lobCoor(5) = un
            lobWeight(1) = 0.1d0
            lobWeight(2) = 49.d0/90.d0
            lobWeight(3) = 32.0/45.d0
            lobWeight(4) = lobWeight(2)
            lobWeight(5) = lobWeight(1)
            do iz = 1, 5
                xpg(iz) = un/3.d0
                ypg(iz) = un/3.d0
                zpg(iz) = lobCoor(iz)
                hpg(iz) = lobWeight(iz)*undemi
            end do

        else if (fapg .eq. 'LOB7') then
! --------- FORMULE DE QUADRATURE DE GAUSS-LOBATTO A 7 POINTS DANS
!           L EPAISSEUR AU CENTRE DE L'ELEMENT
            lobCoor(1) = -un
            lobCoor(2) = -lobatto7p26
            lobCoor(3) = -lobatto7p35
            lobCoor(4) = zero
            lobCoor(5) = -lobCoor(3)
            lobCoor(6) = -lobCoor(2)
            lobCoor(7) = -lobCoor(1)
            lobWeight(1) = un/21.d0
            lobWeight(2) = (124.d0-7.d0*rac15)/350.d0
            lobWeight(3) = (124.d0+7.d0*rac15)/350.d0
            lobWeight(4) = 256.d0/525.d0
            lobWeight(5) = lobWeight(3)
            lobWeight(6) = lobWeight(2)
            lobWeight(7) = lobWeight(1)
            do iz = 1, 7
                xpg(iz) = un/3.d0
                ypg(iz) = un/3.d0
                zpg(iz) = lobCoor(iz)
                hpg(iz) = lobWeight(iz)
            end do

        else
            ASSERT(ASTER_FALSE)
        end if

    else if (elrefa .eq. 'TE4' .or. elrefa .eq. 'T10' &
             .or. elrefa .eq. 'T15' .or. elrefa .eq. 'T20') then
        if (fapg .eq. 'FPG4') then
! --------- FORMULE A 4 POINTS :  (CF TOUZOT PAGE 300) - ORDRE 2 EN X Y Z
            a1 = (5.d0-rac5)/20.d0
            b1 = (5.d0+3.d0*rac5)/20.d0
            h5 = un/24.d0
            npi = 0
            do i = 1, 4
                npi = npi+1
                xpg(npi) = a1
                ypg(npi) = a1
                zpg(npi) = a1
                hpg(npi) = h5
            end do
            zpg(2) = b1
            ypg(3) = b1
            xpg(4) = b1

        else if (fapg .eq. 'FPG5') then
! --------- FORMULE A 5 POINTS :  (CF TOUZOT PAGE 300) - ORDRE 3 EN X Y Z
            a1 = 0.25d0
            b1 = un/6.d0
            c1 = undemi
            h1 = -deux/15.d0
            h2 = 3.d0/40.d0
            xpg(1) = a1
            ypg(1) = a1
            zpg(1) = a1
            hpg(1) = h1
            do i = 2, 5
                xpg(i) = b1
                ypg(i) = b1
                zpg(i) = b1
                hpg(i) = h2
            end do
            zpg(3) = c1
            ypg(4) = c1
            xpg(5) = c1

        else if (fapg .eq. 'FPG11') then
! --------- FORMULE A 11 POINTS :  (CF DUNAVANT) - ORDRE 4 EN X Y Z
            xpg(1) = 0.097204644587583d0
            ypg(1) = 0.106604172561993d0
            zpg(1) = 0.684390415453040d0
            hpg(1) = 0.106468034155490d0/6.d0
            xpg(2) = 0.029569495206479d0
            ypg(2) = 0.329232959742646d0
            zpg(2) = 0.317903560213394d0
            hpg(2) = 0.110234232428497d0/6.d0
            xpg(3) = 0.432710239047768d0
            ypg(3) = 0.103844116410993d0
            zpg(3) = 0.353823239209297d0
            hpg(3) = 0.154976116016246d0/6.d0
            xpg(4) = 0.240276664928072d0
            ypg(4) = 0.304448402434496d0
            zpg(4) = 0.126801725915392d0
            hpg(4) = 0.193410812049634d0/6.d0
            xpg(5) = 0.129411373788910d0
            ypg(5) = 0.538007203916185d0
            zpg(5) = 0.330190414837464d0
            hpg(5) = 0.076162715245558d0/6.d0
            xpg(6) = 0.121541991333927d0
            ypg(6) = 0.008991260093335d0
            zpg(6) = 0.306493988429690d0
            hpg(6) = 0.079426680068025d0/6.d0
            xpg(7) = 0.450765876091276d0
            ypg(7) = 0.432953490481355d0
            zpg(7) = 0.059456616299433d0
            hpg(7) = 0.069469965937635d0/6.d0
            xpg(8) = 0.419266313879513d0
            ypg(8) = 0.053341239535745d0
            zpg(8) = 0.047781435559086d0
            hpg(8) = 0.059933185146559d0/6.d0
            xpg(9) = 0.067223294893383d0
            ypg(9) = 0.741228882093622d0
            zpg(9) = 0.035183929773598d0
            hpg(9) = 0.055393798871576d0/6.d0
            xpg(10) = 0.752508507009654d0
            ypg(10) = 0.081404918402859d0
            zpg(10) = 0.068099370938206d0
            hpg(10) = 0.055273369155936d0/6.d0
            xpg(11) = 0.040490506727590d0
            ypg(11) = 0.174694058697230d0
            zpg(11) = 0.013560701879802d0
            hpg(11) = 0.039251090924839d0/6.d0

        else if (fapg .eq. 'FPG15') then
! --------- FORMULE A 15 POINTS :  (CF TOUZOT PAGE 300) - ORDRE 5 EN X Y Z
            a1 = 0.25d0
            b1 = (7.0d0+rac15)/34.0d0
            b2 = (7.0d0-rac15)/34.0d0
            c1 = (13.0d0-3.0d0*rac15)/34.0d0
            c2 = (13.0d0+3.0d0*rac15)/34.0d0
            d1 = (5.0d0-rac15)/20.0d0
            e1 = (5.0d0+rac15)/20.0d0
            h5 = 8.0d0/405.0d0
            h1 = (2665.0d0-14.0d0*rac15)/226800.0d0
            h2 = (2665.0d0+14.0d0*rac15)/226800.0d0
            h3 = 5.0d0/567.0d0

            xpg(1) = a1
            ypg(1) = a1
            zpg(1) = a1
            hpg(1) = h5
            xpg(2) = b1
            ypg(2) = b1
            zpg(2) = b1
            hpg(2) = h1
            xpg(3) = b1
            ypg(3) = b1
            zpg(3) = c1
            hpg(3) = h1
            xpg(4) = b1
            ypg(4) = c1
            zpg(4) = b1
            hpg(4) = h1
            xpg(5) = c1
            ypg(5) = b1
            zpg(5) = b1
            hpg(5) = h1
            xpg(6) = b2
            ypg(6) = b2
            zpg(6) = b2
            hpg(6) = h2
            xpg(7) = b2
            ypg(7) = b2
            zpg(7) = c2
            hpg(7) = h2
            xpg(8) = b2
            ypg(8) = c2
            zpg(8) = b2
            hpg(8) = h2
            xpg(9) = c2
            ypg(9) = b2
            zpg(9) = b2
            hpg(9) = h2
            xpg(10) = d1
            ypg(10) = d1
            zpg(10) = e1
            hpg(10) = h3
            xpg(11) = d1
            ypg(11) = e1
            zpg(11) = d1
            hpg(11) = h3
            xpg(12) = e1
            ypg(12) = d1
            zpg(12) = d1
            hpg(12) = h3
            xpg(13) = d1
            ypg(13) = e1
            zpg(13) = e1
            hpg(13) = h3
            xpg(14) = e1
            ypg(14) = d1
            zpg(14) = e1
            hpg(14) = h3
            xpg(15) = e1
            ypg(15) = e1
            zpg(15) = d1
            hpg(15) = h3

        else if (fapg .eq. 'FPG24') then
! --------- Order 6 : High-order symmetric cubature rules for tetrahedra and pyramids
!                     10.1002/nme.6528
            xpg(1) = 0.040673958534611353d0
            ypg(1) = 0.040673958534611353d0
            zpg(1) = 0.040673958534611353d0
            hpg(1) = 0.001679535175886774d0
            xpg(2) = 0.877978124396165940d0
            ypg(2) = 0.040673958534611353d0
            zpg(2) = 0.040673958534611353d0
            hpg(2) = 0.001679535175886774d0
            xpg(3) = 0.040673958534611353d0
            ypg(3) = 0.877978124396165940d0
            zpg(3) = 0.040673958534611353d0
            hpg(3) = 0.001679535175886774d0
            xpg(4) = 0.040673958534611353d0
            ypg(4) = 0.040673958534611353d0
            zpg(4) = 0.877978124396165940d0
            hpg(4) = 0.001679535175886774d0
            xpg(5) = 0.322337890142275510d0
            ypg(5) = 0.322337890142275510d0
            zpg(5) = 0.322337890142275510d0
            hpg(5) = 0.009226196923942455d0
            xpg(6) = 0.032986329573173468d0
            ypg(6) = 0.322337890142275510d0
            zpg(6) = 0.322337890142275510d0
            hpg(6) = 0.009226196923942455d0
            xpg(7) = 0.322337890142275510d0
            ypg(7) = 0.032986329573173468d0
            zpg(7) = 0.322337890142275510d0
            hpg(7) = 0.009226196923942455d0
            xpg(8) = 0.322337890142275510d0
            ypg(8) = 0.322337890142275510d0
            zpg(8) = 0.032986329573173468d0
            hpg(8) = 0.009226196923942455d0
            xpg(9) = 0.214602871259152029d0
            ypg(9) = 0.214602871259152029d0
            zpg(9) = 0.214602871259152029d0
            hpg(9) = 0.006653791709694583d0
            xpg(10) = 0.356191386222543912d0
            ypg(10) = 0.214602871259152029d0
            zpg(10) = 0.214602871259152029d0
            hpg(10) = 0.006653791709694583d0
            xpg(11) = 0.214602871259152029d0
            ypg(11) = 0.356191386222543912d0
            zpg(11) = 0.214602871259152029d0
            hpg(11) = 0.006653791709694583d0
            xpg(12) = 0.214602871259152029d0
            ypg(12) = 0.214602871259152029d0
            zpg(12) = 0.356191386222543912d0
            hpg(12) = 0.006653791709694583d0
            xpg(13) = 0.603005664791649141d0
            ypg(13) = 0.063661001875017525d0
            zpg(13) = 0.063661001875017525d0
            hpg(13) = 0.008035714285714287d0
            xpg(14) = 0.063661001875017525d0
            ypg(14) = 0.603005664791649141d0
            zpg(14) = 0.063661001875017525d0
            hpg(14) = 0.008035714285714287d0
            xpg(15) = 0.063661001875017525d0
            ypg(15) = 0.063661001875017525d0
            zpg(15) = 0.603005664791649141d0
            hpg(15) = 0.008035714285714287d0
            xpg(16) = 0.269672331458315808d0
            ypg(16) = 0.063661001875017525d0
            zpg(16) = 0.063661001875017525d0
            hpg(16) = 0.008035714285714287d0
            xpg(17) = 0.269672331458315808d0
            ypg(17) = 0.603005664791649141d0
            zpg(17) = 0.063661001875017525d0
            hpg(17) = 0.008035714285714287d0
            xpg(18) = 0.269672331458315808d0
            ypg(18) = 0.063661001875017525d0
            zpg(18) = 0.603005664791649141d0
            hpg(18) = 0.008035714285714287d0
            xpg(19) = 0.063661001875017525d0
            ypg(19) = 0.269672331458315808d0
            zpg(19) = 0.063661001875017525d0
            hpg(19) = 0.008035714285714287d0
            xpg(20) = 0.603005664791649141d0
            ypg(20) = 0.269672331458315808d0
            zpg(20) = 0.063661001875017525d0
            hpg(20) = 0.008035714285714287d0
            xpg(21) = 0.063661001875017525d0
            ypg(21) = 0.269672331458315808d0
            zpg(21) = 0.603005664791649141d0
            hpg(21) = 0.008035714285714287d0
            xpg(22) = 0.063661001875017525d0
            ypg(22) = 0.063661001875017525d0
            zpg(22) = 0.269672331458315808d0
            hpg(22) = 0.008035714285714287d0
            xpg(23) = 0.603005664791649141d0
            ypg(23) = 0.063661001875017525d0
            zpg(23) = 0.269672331458315808d0
            hpg(23) = 0.008035714285714287d0
            xpg(24) = 0.063661001875017525d0
            ypg(24) = 0.603005664791649141d0
            zpg(24) = 0.269672331458315808d0
            hpg(24) = 0.008035714285714287d0

        else if (fapg .eq. 'FPG35') then
! --------- Order 7 : High-order symmetric cubature rules for tetrahedra and pyramids
!                     10.1002/nme.6528
            xpg(1) = 0.25d0
            ypg(1) = 0.25d0
            zpg(1) = 0.25d0
            hpg(1) = 0.015914214910688475d0
            xpg(2) = 0.315701149778202799d0
            ypg(2) = 0.315701149778202799d0
            zpg(2) = 0.315701149778202799d0
            hpg(2) = 0.007054930201661171d0
            xpg(3) = 0.052896550665391601d0
            ypg(3) = 0.315701149778202799d0
            zpg(3) = 0.315701149778202799d0
            hpg(3) = 0.007054930201661171d0
            xpg(4) = 0.315701149778202799d0
            ypg(4) = 0.052896550665391601d0
            zpg(4) = 0.315701149778202799d0
            hpg(4) = 0.007054930201661171d0
            xpg(5) = 0.315701149778202799d0
            ypg(5) = 0.315701149778202799d0
            zpg(5) = 0.052896550665391601d0
            hpg(5) = 0.007054930201661171d0
            xpg(6) = 0.050489822598396368d0
            ypg(6) = 0.449510177401603631d0
            zpg(6) = 0.449510177401603631d0
            hpg(6) = 0.005316154638809596d0
            xpg(7) = 0.449510177401603631d0
            ypg(7) = 0.050489822598396368d0
            zpg(7) = 0.449510177401603631d0
            hpg(7) = 0.005316154638809596d0
            xpg(8) = 0.050489822598396368d0
            ypg(8) = 0.050489822598396368d0
            zpg(8) = 0.449510177401603631d0
            hpg(8) = 0.005316154638809596d0
            xpg(9) = 0.449510177401603631d0
            ypg(9) = 0.449510177401603631d0
            zpg(9) = 0.050489822598396368d0
            hpg(9) = 0.005316154638809596d0
            xpg(10) = 0.050489822598396368d0
            ypg(10) = 0.449510177401603631d0
            zpg(10) = 0.050489822598396368d0
            hpg(10) = 0.005316154638809596d0
            xpg(11) = 0.449510177401603631d0
            ypg(11) = 0.050489822598396368d0
            zpg(11) = 0.050489822598396368d0
            hpg(11) = 0.005316154638809596d0
            xpg(12) = 0.810830241098548561d0
            ypg(12) = 0.021265472541483245d0
            zpg(12) = 0.021265472541483245d0
            hpg(12) = 0.001351795138317223d0
            xpg(13) = 0.021265472541483245d0
            ypg(13) = 0.810830241098548561d0
            zpg(13) = 0.021265472541483245d0
            hpg(13) = 0.001351795138317223d0
            xpg(14) = 0.021265472541483245d0
            ypg(14) = 0.021265472541483245d0
            zpg(14) = 0.810830241098548561d0
            hpg(14) = 0.001351795138317223d0
            xpg(15) = 0.146638813818484946d0
            ypg(15) = 0.021265472541483245d0
            zpg(15) = 0.021265472541483245d0
            hpg(15) = 0.001351795138317223d0
            xpg(16) = 0.146638813818484946d0
            ypg(16) = 0.810830241098548561d0
            zpg(16) = 0.021265472541483245d0
            hpg(16) = 0.001351795138317223d0
            xpg(17) = 0.146638813818484946d0
            ypg(17) = 0.021265472541483245d0
            zpg(17) = 0.810830241098548561d0
            hpg(17) = 0.001351795138317223d0
            xpg(18) = 0.021265472541483245d0
            ypg(18) = 0.146638813818484946d0
            zpg(18) = 0.021265472541483245d0
            hpg(18) = 0.001351795138317223d0
            xpg(19) = 0.810830241098548561d0
            ypg(19) = 0.146638813818484946d0
            zpg(19) = 0.021265472541483245d0
            hpg(19) = 0.001351795138317223d0
            xpg(20) = 0.021265472541483245d0
            ypg(20) = 0.146638813818484946d0
            zpg(20) = 0.810830241098548561d0
            hpg(20) = 0.001351795138317223d0
            xpg(21) = 0.021265472541483245d0
            ypg(21) = 0.021265472541483245d0
            zpg(21) = 0.146638813818484946d0
            hpg(21) = 0.001351795138317223d0
            xpg(22) = 0.810830241098548561d0
            ypg(22) = 0.021265472541483245d0
            zpg(22) = 0.146638813818484946d0
            hpg(22) = 0.001351795138317223d0
            xpg(23) = 0.021265472541483245d0
            ypg(23) = 0.810830241098548561d0
            zpg(23) = 0.146638813818484946d0
            hpg(23) = 0.001351795138317223d0
            xpg(24) = 0.575171637587000023d0
            ypg(24) = 0.188833831026001047d0
            zpg(24) = 0.188833831026001047d0
            hpg(24) = 0.006201188454722437d0
            xpg(25) = 0.188833831026001047d0
            ypg(25) = 0.575171637587000023d0
            zpg(25) = 0.188833831026001047d0
            hpg(25) = 0.006201188454722437d0
            xpg(26) = 0.188833831026001047d0
            ypg(26) = 0.188833831026001047d0
            zpg(26) = 0.575171637587000023d0
            hpg(26) = 0.006201188454722437d0
            xpg(27) = 0.047160700360997881d0
            ypg(27) = 0.188833831026001047d0
            zpg(27) = 0.188833831026001047d0
            hpg(27) = 0.006201188454722437d0
            xpg(28) = 0.047160700360997881d0
            ypg(28) = 0.575171637587000023d0
            zpg(28) = 0.188833831026001047d0
            hpg(28) = 0.006201188454722437d0
            xpg(29) = 0.047160700360997881d0
            ypg(29) = 0.188833831026001047d0
            zpg(29) = 0.575171637587000023d0
            hpg(29) = 0.006201188454722437d0
            xpg(30) = 0.188833831026001047d0
            ypg(30) = 0.047160700360997881d0
            zpg(30) = 0.188833831026001047d0
            hpg(30) = 0.006201188454722437d0
            xpg(31) = 0.575171637587000023d0
            ypg(31) = 0.047160700360997881d0
            zpg(31) = 0.188833831026001047d0
            hpg(31) = 0.006201188454722437d0
            xpg(32) = 0.188833831026001047d0
            ypg(32) = 0.047160700360997881d0
            zpg(32) = 0.575171637587000023d0
            hpg(32) = 0.006201188454722437d0
            xpg(33) = 0.188833831026001047d0
            ypg(33) = 0.188833831026001047d0
            zpg(33) = 0.047160700360997881d0
            hpg(33) = 0.006201188454722437d0
            xpg(34) = 0.575171637587000023d0
            ypg(34) = 0.188833831026001047d0
            zpg(34) = 0.047160700360997881d0
            hpg(34) = 0.006201188454722437d0
            xpg(35) = 0.188833831026001047d0
            ypg(35) = 0.575171637587000023d0
            zpg(35) = 0.047160700360997881d0
            hpg(35) = 0.006201188454722437d0

        else if (fapg .eq. 'FPG46') then
! --------- Order 8 : High-order symmetric cubature rules for tetrahedra and pyramids
!                     10.1002/nme.6528
            xpg(1) = 0.186812758270720971d0
            ypg(1) = 0.186812758270720971d0
            zpg(1) = 0.186812758270720971d0
            hpg(1) = 0.00796883133905111d0
            xpg(2) = 0.439561725187837086d0
            ypg(2) = 0.186812758270720971d0
            zpg(2) = 0.186812758270720971d0
            hpg(2) = 0.00796883133905111d0
            xpg(3) = 0.186812758270720971d0
            ypg(3) = 0.439561725187837086d0
            zpg(3) = 0.186812758270720971d0
            hpg(3) = 0.00796883133905111d0
            xpg(4) = 0.186812758270720971d0
            ypg(4) = 0.186812758270720971d0
            zpg(4) = 0.439561725187837086d0
            hpg(4) = 0.00796883133905111d0
            xpg(5) = 0.114462406761230491d0
            ypg(5) = 0.114462406761230491d0
            zpg(5) = 0.114462406761230491d0
            hpg(5) = 0.004954642347799037d0
            xpg(6) = 0.656612779716308524d0
            ypg(6) = 0.114462406761230491d0
            zpg(6) = 0.114462406761230491d0
            hpg(6) = 0.004954642347799037d0
            xpg(7) = 0.114462406761230491d0
            ypg(7) = 0.656612779716308524d0
            zpg(7) = 0.114462406761230491d0
            hpg(7) = 0.004954642347799037d0
            xpg(8) = 0.114462406761230491d0
            ypg(8) = 0.114462406761230491d0
            zpg(8) = 0.656612779716308524d0
            hpg(8) = 0.004954642347799037d0
            xpg(9) = 0.313806592272401620d0
            ypg(9) = 0.313806592272401620d0
            zpg(9) = 0.313806592272401620d0
            hpg(9) = 0.007254554036495878d0
            xpg(10) = 0.058580223182795139d0
            ypg(10) = 0.313806592272401620d0
            zpg(10) = 0.313806592272401620d0
            hpg(10) = 0.007254554036495878d0
            xpg(11) = 0.313806592272401620d0
            ypg(11) = 0.058580223182795139d0
            zpg(11) = 0.313806592272401620d0
            hpg(11) = 0.007254554036495878d0
            xpg(12) = 0.313806592272401620d0
            ypg(12) = 0.313806592272401620d0
            zpg(12) = 0.058580223182795139d0
            hpg(12) = 0.007254554036495878d0
            xpg(13) = 0.044577379448846248d0
            ypg(13) = 0.044577379448846248d0
            zpg(13) = 0.044577379448846248d0
            hpg(13) = 0.001438789336802495d0
            xpg(14) = 0.866267861653461255d0
            ypg(14) = 0.044577379448846248d0
            zpg(14) = 0.044577379448846248d0
            hpg(14) = 0.001438789336802495d0
            xpg(15) = 0.044577379448846248d0
            ypg(15) = 0.866267861653461255d0
            zpg(15) = 0.044577379448846248d0
            hpg(15) = 0.001438789336802495d0
            xpg(16) = 0.044577379448846248d0
            ypg(16) = 0.044577379448846248d0
            zpg(16) = 0.866267861653461255d0
            hpg(16) = 0.001438789336802495d0
            xpg(17) = 0.434516509456157605d0
            ypg(17) = 0.065483490543842394d0
            zpg(17) = 0.065483490543842394d0
            hpg(17) = 0.006148423544977727d0
            xpg(18) = 0.065483490543842394d0
            ypg(18) = 0.434516509456157605d0
            zpg(18) = 0.065483490543842394d0
            hpg(18) = 0.006148423544977727d0
            xpg(19) = 0.434516509456157605d0
            ypg(19) = 0.434516509456157605d0
            zpg(19) = 0.065483490543842394d0
            hpg(19) = 0.006148423544977727d0
            xpg(20) = 0.065483490543842394d0
            ypg(20) = 0.065483490543842394d0
            zpg(20) = 0.434516509456157605d0
            hpg(20) = 0.006148423544977727d0
            xpg(21) = 0.434516509456157605d0
            ypg(21) = 0.065483490543842394d0
            zpg(21) = 0.434516509456157605d0
            hpg(21) = 0.006148423544977727d0
            xpg(22) = 0.065483490543842394d0
            ypg(22) = 0.434516509456157605d0
            zpg(22) = 0.434516509456157605d0
            hpg(22) = 0.006148423544977727d0
            xpg(23) = 0.005073320750421590d0
            ypg(23) = 0.203893174662110274d0
            zpg(23) = 0.203893174662110274d0
            hpg(23) = 0.002416171118476584d0
            xpg(24) = 0.203893174662110274d0
            ypg(24) = 0.005073320750421590d0
            zpg(24) = 0.203893174662110274d0
            hpg(24) = 0.002416171118476584d0
            xpg(25) = 0.203893174662110274d0
            ypg(25) = 0.203893174662110274d0
            zpg(25) = 0.005073320750421590d0
            hpg(25) = 0.002416171118476584d0
            xpg(26) = 0.587140329925357860d0
            ypg(26) = 0.203893174662110274d0
            zpg(26) = 0.203893174662110274d0
            hpg(26) = 0.002416171118476584d0
            xpg(27) = 0.587140329925357860d0
            ypg(27) = 0.005073320750421590d0
            zpg(27) = 0.203893174662110274d0
            hpg(27) = 0.002416171118476584d0
            xpg(28) = 0.587140329925357860d0
            ypg(28) = 0.203893174662110274d0
            zpg(28) = 0.005073320750421590d0
            hpg(28) = 0.002416171118476584d0
            xpg(29) = 0.203893174662110274d0
            ypg(29) = 0.587140329925357860d0
            zpg(29) = 0.203893174662110274d0
            hpg(29) = 0.002416171118476584d0
            xpg(30) = 0.005073320750421590d0
            ypg(30) = 0.587140329925357860d0
            zpg(30) = 0.203893174662110274d0
            hpg(30) = 0.002416171118476584d0
            xpg(31) = 0.203893174662110274d0
            ypg(31) = 0.587140329925357860d0
            zpg(31) = 0.005073320750421590d0
            hpg(31) = 0.002416171118476584d0
            xpg(32) = 0.203893174662110274d0
            ypg(32) = 0.203893174662110274d0
            zpg(32) = 0.587140329925357860d0
            hpg(32) = 0.002416171118476584d0
            xpg(33) = 0.005073320750421590d0
            ypg(33) = 0.203893174662110274d0
            zpg(33) = 0.587140329925357860d0
            hpg(33) = 0.002416171118476584d0
            xpg(34) = 0.203893174662110274d0
            ypg(34) = 0.005073320750421590d0
            zpg(34) = 0.587140329925357860d0
            hpg(34) = 0.002416171118476584d0
            xpg(35) = 0.714676926393304887d0
            ypg(35) = 0.021247779669571673d0
            zpg(35) = 0.021247779669571673d0
            hpg(35) = 0.001192900311207267d0
            xpg(36) = 0.021247779669571673d0
            ypg(36) = 0.714676926393304887d0
            zpg(36) = 0.021247779669571673d0
            hpg(36) = 0.001192900311207267d0
            xpg(37) = 0.021247779669571673d0
            ypg(37) = 0.021247779669571673d0
            zpg(37) = 0.714676926393304887d0
            hpg(37) = 0.001192900311207267d0
            xpg(38) = 0.242827514267551765d0
            ypg(38) = 0.021247779669571673d0
            zpg(38) = 0.021247779669571673d0
            hpg(38) = 0.001192900311207267d0
            xpg(39) = 0.242827514267551765d0
            ypg(39) = 0.714676926393304887d0
            zpg(39) = 0.021247779669571673d0
            hpg(39) = 0.001192900311207267d0
            xpg(40) = 0.242827514267551765d0
            ypg(40) = 0.021247779669571673d0
            zpg(40) = 0.714676926393304887d0
            hpg(40) = 0.001192900311207267d0
            xpg(41) = 0.021247779669571673d0
            ypg(41) = 0.242827514267551765d0
            zpg(41) = 0.021247779669571673d0
            hpg(41) = 0.001192900311207267d0
            xpg(42) = 0.714676926393304887d0
            ypg(42) = 0.242827514267551765d0
            zpg(42) = 0.021247779669571673d0
            hpg(42) = 0.001192900311207267d0
            xpg(43) = 0.021247779669571673d0
            ypg(43) = 0.242827514267551765d0
            zpg(43) = 0.714676926393304887d0
            hpg(43) = 0.001192900311207267d0
            xpg(44) = 0.021247779669571673d0
            ypg(44) = 0.021247779669571673d0
            zpg(44) = 0.242827514267551765d0
            hpg(44) = 0.001192900311207267d0
            xpg(45) = 0.714676926393304887d0
            ypg(45) = 0.021247779669571673d0
            zpg(45) = 0.242827514267551765d0
            hpg(45) = 0.001192900311207267d0
            xpg(46) = 0.021247779669571673d0
            ypg(46) = 0.714676926393304887d0
            zpg(46) = 0.242827514267551765d0
            hpg(46) = 0.001192900311207267d0

        else if (fapg .eq. 'FPG59') then
! --------- Order 9 : High-order symmetric cubature rules for tetrahedra and pyramids
!           10.1002/nme.6528
            xpg(1) = 0.25d0
            ypg(1) = 0.25d0
            zpg(1) = 0.25d0
            hpg(1) = 0.009477770708925288d0
            xpg(2) = 0.167906605236742789d0
            ypg(2) = 0.167906605236742789d0
            zpg(2) = 0.167906605236742789d0
            hpg(2) = 0.004282379448253863d0
            xpg(3) = 0.496280184289771631d0
            ypg(3) = 0.167906605236742789d0
            zpg(3) = 0.167906605236742789d0
            hpg(3) = 0.004282379448253863d0
            xpg(4) = 0.167906605236742789d0
            ypg(4) = 0.496280184289771631d0
            zpg(4) = 0.167906605236742789d0
            hpg(4) = 0.004282379448253863d0
            xpg(5) = 0.167906605236742789d0
            ypg(5) = 0.167906605236742789d0
            zpg(5) = 0.496280184289771631d0
            hpg(5) = 0.004282379448253863d0
            xpg(6) = 0.091436270514076258d0
            ypg(6) = 0.091436270514076258d0
            zpg(6) = 0.091436270514076258d0
            hpg(6) = 0.000383153892171395d0
            xpg(7) = 0.725691188457771224d0
            ypg(7) = 0.091436270514076258d0
            zpg(7) = 0.091436270514076258d0
            hpg(7) = 0.000383153892171395d0
            xpg(8) = 0.091436270514076258d0
            ypg(8) = 0.725691188457771224d0
            zpg(8) = 0.091436270514076258d0
            hpg(8) = 0.000383153892171395d0
            xpg(9) = 0.091436270514076258d0
            ypg(9) = 0.091436270514076258d0
            zpg(9) = 0.725691188457771224d0
            hpg(9) = 0.000383153892171395d0
            xpg(10) = 0.321855664817653266d0
            ypg(10) = 0.321855664817653266d0
            zpg(10) = 0.321855664817653266d0
            hpg(10) = 0.005076006444598501d0
            xpg(11) = 0.034433005547040201d0
            ypg(11) = 0.321855664817653266d0
            zpg(11) = 0.321855664817653266d0
            hpg(11) = 0.005076006444598501d0
            xpg(12) = 0.321855664817653266d0
            ypg(12) = 0.034433005547040201d0
            zpg(12) = 0.321855664817653266d0
            hpg(12) = 0.005076006444598501d0
            xpg(13) = 0.321855664817653266d0
            ypg(13) = 0.321855664817653266d0
            zpg(13) = 0.034433005547040201d0
            hpg(13) = 0.005076006444598501d0
            xpg(14) = 0.041837695900365602d0
            ypg(14) = 0.041837695900365602d0
            zpg(14) = 0.041837695900365602d0
            hpg(14) = 0.001187287056706480d0
            xpg(15) = 0.874486912298903193d0
            ypg(15) = 0.041837695900365602d0
            zpg(15) = 0.041837695900365602d0
            hpg(15) = 0.001187287056706480d0
            xpg(16) = 0.041837695900365602d0
            ypg(16) = 0.874486912298903193d0
            zpg(16) = 0.041837695900365602d0
            hpg(16) = 0.001187287056706480d0
            xpg(17) = 0.041837695900365602d0
            ypg(17) = 0.041837695900365602d0
            zpg(17) = 0.874486912298903193d0
            hpg(17) = 0.001187287056706480d0
            xpg(18) = 0.393270561425153639d0
            ypg(18) = 0.106729438574846360d0
            zpg(18) = 0.106729438574846360d0
            hpg(18) = 0.006140609246823980d0
            xpg(19) = 0.106729438574846360d0
            ypg(19) = 0.393270561425153639d0
            zpg(19) = 0.106729438574846360d0
            hpg(19) = 0.006140609246823980d0
            xpg(20) = 0.393270561425153639d0
            ypg(20) = 0.393270561425153639d0
            zpg(20) = 0.106729438574846360d0
            hpg(20) = 0.006140609246823980d0
            xpg(21) = 0.106729438574846360d0
            ypg(21) = 0.106729438574846360d0
            zpg(21) = 0.393270561425153639d0
            hpg(21) = 0.006140609246823980d0
            xpg(22) = 0.393270561425153639d0
            ypg(22) = 0.106729438574846360d0
            zpg(22) = 0.393270561425153639d0
            hpg(22) = 0.006140609246823980d0
            xpg(23) = 0.106729438574846360d0
            ypg(23) = 0.393270561425153639d0
            zpg(23) = 0.393270561425153639d0
            hpg(23) = 0.006140609246823980d0
            xpg(24) = 0.718393093981424399d0
            ypg(24) = 0.033957168183088885d0
            zpg(24) = 0.033957168183088885d0
            hpg(24) = 0.001720342797034974d0
            xpg(25) = 0.033957168183088885d0
            ypg(25) = 0.718393093981424399d0
            zpg(25) = 0.033957168183088885d0
            hpg(25) = 0.001720342797034974d0
            xpg(26) = 0.033957168183088885d0
            ypg(26) = 0.033957168183088885d0
            zpg(26) = 0.718393093981424399d0
            hpg(26) = 0.001720342797034974d0
            xpg(27) = 0.213692569652397829d0
            ypg(27) = 0.033957168183088885d0
            zpg(27) = 0.033957168183088885d0
            hpg(27) = 0.001720342797034974d0
            xpg(28) = 0.213692569652397829d0
            ypg(28) = 0.718393093981424399d0
            zpg(28) = 0.033957168183088885d0
            hpg(28) = 0.001720342797034974d0
            xpg(29) = 0.213692569652397829d0
            ypg(29) = 0.033957168183088885d0
            zpg(29) = 0.718393093981424399d0
            hpg(29) = 0.001720342797034974d0
            xpg(30) = 0.033957168183088885d0
            ypg(30) = 0.213692569652397829d0
            zpg(30) = 0.033957168183088885d0
            hpg(30) = 0.001720342797034974d0
            xpg(31) = 0.718393093981424399d0
            ypg(31) = 0.213692569652397829d0
            zpg(31) = 0.033957168183088885d0
            hpg(31) = 0.001720342797034974d0
            xpg(32) = 0.033957168183088885d0
            ypg(32) = 0.213692569652397829d0
            zpg(32) = 0.718393093981424399d0
            hpg(32) = 0.001720342797034974d0
            xpg(33) = 0.033957168183088885d0
            ypg(33) = 0.033957168183088885d0
            zpg(33) = 0.213692569652397829d0
            hpg(33) = 0.001720342797034974d0
            xpg(34) = 0.718393093981424399d0
            ypg(34) = 0.033957168183088885d0
            zpg(34) = 0.213692569652397829d0
            hpg(34) = 0.001720342797034974d0
            xpg(35) = 0.033957168183088885d0
            ypg(35) = 0.718393093981424399d0
            zpg(35) = 0.213692569652397829d0
            hpg(35) = 0.001720342797034974d0
            xpg(36) = 0.078683634476687930d0
            ypg(36) = 0.460658181054777632d0
            zpg(36) = 0.460658181054777632d0
            hpg(36) = 0.001272481192330142d0
            xpg(37) = 0.460658181054777632d0
            ypg(37) = 0.078683634476687930d0
            zpg(37) = 0.460658181054777632d0
            hpg(37) = 0.001272481192330142d0
            xpg(38) = 0.460658181054777632d0
            ypg(38) = 0.460658181054777632d0
            zpg(38) = 0.078683634476687930d0
            hpg(38) = 0.001272481192330142d0
            xpg(39) = 0.000000003413756803d0
            ypg(39) = 0.460658181054777632d0
            zpg(39) = 0.460658181054777632d0
            hpg(39) = 0.001272481192330142d0
            xpg(40) = 0.000000003413756803d0
            ypg(40) = 0.078683634476687930d0
            zpg(40) = 0.460658181054777632d0
            hpg(40) = 0.001272481192330142d0
            xpg(41) = 0.000000003413756803d0
            ypg(41) = 0.460658181054777632d0
            zpg(41) = 0.078683634476687930d0
            hpg(41) = 0.001272481192330142d0
            xpg(42) = 0.460658181054777632d0
            ypg(42) = 0.000000003413756803d0
            zpg(42) = 0.460658181054777632d0
            hpg(42) = 0.001272481192330142d0
            xpg(43) = 0.078683634476687930d0
            ypg(43) = 0.000000003413756803d0
            zpg(43) = 0.460658181054777632d0
            hpg(43) = 0.001272481192330142d0
            xpg(44) = 0.460658181054777632d0
            ypg(44) = 0.000000003413756803d0
            zpg(44) = 0.078683634476687930d0
            hpg(44) = 0.001272481192330142d0
            xpg(45) = 0.460658181054777632d0
            ypg(45) = 0.460658181054777632d0
            zpg(45) = 0.000000003413756803d0
            hpg(45) = 0.001272481192330142d0
            xpg(46) = 0.078683634476687930d0
            ypg(46) = 0.460658181054777632d0
            zpg(46) = 0.000000003413756803d0
            hpg(46) = 0.001272481192330142d0
            xpg(47) = 0.460658181054777632d0
            ypg(47) = 0.078683634476687930d0
            zpg(47) = 0.000000003413756803d0
            hpg(47) = 0.001272481192330142d0
            xpg(48) = 0.597210722761849931d0
            ypg(48) = 0.184387943500015203d0
            zpg(48) = 0.184387943500015203d0
            hpg(48) = 0.003393003769791262d0
            xpg(49) = 0.184387943500015203d0
            ypg(49) = 0.597210722761849931d0
            zpg(49) = 0.184387943500015203d0
            hpg(49) = 0.003393003769791262d0
            xpg(50) = 0.184387943500015203d0
            ypg(50) = 0.184387943500015203d0
            zpg(50) = 0.597210722761849931d0
            hpg(50) = 0.003393003769791262d0
            xpg(51) = 0.034013390238119661d0
            ypg(51) = 0.184387943500015203d0
            zpg(51) = 0.184387943500015203d0
            hpg(51) = 0.003393003769791262d0
            xpg(52) = 0.034013390238119661d0
            ypg(52) = 0.597210722761849931d0
            zpg(52) = 0.184387943500015203d0
            hpg(52) = 0.003393003769791262d0
            xpg(53) = 0.034013390238119661d0
            ypg(53) = 0.184387943500015203d0
            zpg(53) = 0.597210722761849931d0
            hpg(53) = 0.003393003769791262d0
            xpg(54) = 0.184387943500015203d0
            ypg(54) = 0.034013390238119661d0
            zpg(54) = 0.184387943500015203d0
            hpg(54) = 0.003393003769791262d0
            xpg(55) = 0.597210722761849931d0
            ypg(55) = 0.034013390238119661d0
            zpg(55) = 0.184387943500015203d0
            hpg(55) = 0.003393003769791262d0
            xpg(56) = 0.184387943500015203d0
            ypg(56) = 0.034013390238119661d0
            zpg(56) = 0.597210722761849931d0
            hpg(56) = 0.003393003769791262d0
            xpg(57) = 0.184387943500015203d0
            ypg(57) = 0.184387943500015203d0
            zpg(57) = 0.034013390238119661d0
            hpg(57) = 0.003393003769791262d0
            xpg(58) = 0.597210722761849931d0
            ypg(58) = 0.184387943500015203d0
            zpg(58) = 0.034013390238119661d0
            hpg(58) = 0.003393003769791262d0
            xpg(59) = 0.184387943500015203d0
            ypg(59) = 0.597210722761849931d0
            zpg(59) = 0.034013390238119661d0
            hpg(59) = 0.003393003769791262d0

        else if (fapg .eq. 'FPG74') then
! --------- Order 10 : High-order cubature rules for tetrahedra - 10.1002/nme.6313
            xpg(1) = 0.710821181163561004d0
            ypg(1) = 0.168093815148608258d0
            zpg(1) = 0.091043261506008295d0
            hpg(1) = 0.002182350019870431d0
            xpg(2) = 0.579007667001554496d0
            ypg(2) = 0.075413341909639651d0
            zpg(2) = 0.193173199165442126d0
            hpg(2) = 0.003872737807375740d0
            xpg(3) = 0.097894295280830007d0
            ypg(3) = 0.519326664905817731d0
            zpg(3) = 0.352023753338244266d0
            hpg(3) = 0.002647995019082234d0
            xpg(4) = 0.071174277974313488d0
            ypg(4) = 0.299085994269323218d0
            zpg(4) = 0.574864439546898916d0
            hpg(4) = 0.003283818998754858d0
            xpg(5) = 0.017880749089293777d0
            ypg(5) = 0.411156457218860277d0
            zpg(5) = 0.027100342875812305d0
            hpg(5) = 0.001014502585258677d0
            xpg(6) = 0.003714135933259649d0
            ypg(6) = 0.487027149727789748d0
            zpg(6) = 0.508667162441936434d0
            hpg(6) = 0.000280175184884064d0
            xpg(7) = 0.165301928762668738d0
            ypg(7) = 0.032214258648313131d0
            zpg(7) = 0.185606096904448124d0
            hpg(7) = 0.003095969596593190d0
            xpg(8) = 0.154834593876529642d0
            ypg(8) = 0.166588606435463824d0
            zpg(8) = 0.100852167920593681d0
            hpg(8) = 0.004310991417454283d0
            xpg(9) = 0.505365346710630347d0
            ypg(9) = 0.305274915800875908d0
            zpg(9) = 0.003529745262325717d0
            hpg(9) = 0.001238315304433764d0
            xpg(10) = 0.349721451595561549d0
            ypg(10) = 0.065040629197533067d0
            zpg(10) = 0.370487586728316691d0
            hpg(10) = 0.004663468429689951d0
            xpg(11) = 0.362129457118497282d0
            ypg(11) = 0.223047943728803069d0
            zpg(11) = 0.211122229756149666d0
            hpg(11) = 0.006629051767911569d0
            xpg(12) = 0.001519081155224517d0
            ypg(12) = 0.411267091538767502d0
            zpg(12) = 0.458570806845948091d0
            hpg(12) = 0.001263035790666400d0
            xpg(13) = 0.367671015888483266d0
            ypg(13) = 0.044230430969319516d0
            zpg(13) = 0.565336903160872039d0
            hpg(13) = 0.001557155171488173d0
            xpg(14) = 0.049709772054708907d0
            ypg(14) = 0.838858847721160245d0
            zpg(14) = 0.107604804604589760d0
            hpg(14) = 0.000517777452563824d0
            xpg(15) = 0.132539336531317979d0
            ypg(15) = 0.085196292116014844d0
            zpg(15) = 0.604412365542924088d0
            hpg(15) = 0.003716006543051116d0
            xpg(16) = 0.016793055088974239d0
            ypg(16) = 0.146737983006163192d0
            zpg(16) = 0.719894257786155718d0
            hpg(16) = 0.001554387915554572d0
            xpg(17) = 0.384644976510008888d0
            ypg(17) = 0.077830854428130452d0
            zpg(17) = 0.133820842318490240d0
            hpg(17) = 0.00566376575084548d0
            xpg(18) = 0.610719347541429062d0
            ypg(18) = 0.031493054229378961d0
            zpg(18) = 0.025645923986776844d0
            hpg(18) = 0.001290087563843645d0
            xpg(19) = 0.080009809180926627d0
            ypg(19) = 0.908233671118126695d0
            zpg(19) = 0.010606985116433330d0
            hpg(19) = 0.000170380133009282d0
            xpg(20) = 0.172211544573607231d0
            ypg(20) = 0.004171136945129776d0
            zpg(20) = 0.788140429262158342d0
            hpg(20) = 0.000526411877337597d0
            xpg(21) = 0.069735326902137974d0
            ypg(21) = 0.465324657892589599d0
            zpg(21) = 0.279743130325619603d0
            hpg(21) = 0.005086264077008329d0
            xpg(22) = 0.611370160489733183d0
            ypg(22) = 0.136777765750491082d0
            zpg(22) = 0.051140189028374378d0
            hpg(22) = 0.003473693066209297d0
            xpg(23) = 0.221116192418793230d0
            ypg(23) = 0.535674881793983028d0
            zpg(23) = 0.010811825572610978d0
            hpg(23) = 0.001798639039181593d0
            xpg(24) = 0.575443581149477263d0
            ypg(24) = 0.382823350763782858d0
            zpg(24) = 0.017764367573077392d0
            hpg(24) = 0.000886812529713127d0
            xpg(25) = 0.115839281892140834d0
            ypg(25) = 0.772788669239955026d0
            zpg(25) = 0.024219507940919968d0
            hpg(25) = 0.001449278087527263d0
            xpg(26) = 0.782798881531462869d0
            ypg(26) = 0.161928384574725930d0
            zpg(26) = 0.001617373723301612d0
            hpg(26) = 0.000593302967809375d0
            xpg(27) = 0.412163414906180645d0
            ypg(27) = 0.120290125453191509d0
            zpg(27) = 0.009527433995123112d0
            hpg(27) = 0.001443276892453364d0
            xpg(28) = 0.181426179587850997d0
            ypg(28) = 0.109629539169664374d0
            zpg(28) = 0.330173323416455686d0
            hpg(28) = 0.005548067340506629d0
            xpg(29) = 0.292868069344387546d0
            ypg(29) = 0.653701294496030294d0
            zpg(29) = 0.029147978424740672d0
            hpg(29) = 0.001093713271038079d0
            xpg(30) = 0.490712771415291857d0
            ypg(30) = 0.311437265416919243d0
            zpg(30) = 0.100667566467669548d0
            hpg(30) = 0.004084725426581246d0
            xpg(31) = 0.029321561613951272d0
            ypg(31) = 0.041402332917342573d0
            zpg(31) = 0.109514718958966341d0
            hpg(31) = 0.001236884663129069d0
            xpg(32) = 0.327603353093800045d0
            ypg(32) = 0.000074149488728507d0
            zpg(32) = 0.295487535102606277d0
            hpg(32) = 0.001397687787374997d0
            xpg(33) = 0.360017861338647169d0
            ypg(33) = 0.501297898701069466d0
            zpg(33) = 0.043250822680646380d0
            hpg(33) = 0.002346807846030675d0
            xpg(34) = 0.775531187938377343d0
            ypg(34) = 0.030148501866560255d0
            zpg(34) = 0.148187911146689746d0
            hpg(34) = 0.001347018140028771d0
            xpg(35) = 0.212622519012012208d0
            ypg(35) = 0.482535036274829435d0
            zpg(35) = 0.130354752098416560d0
            hpg(35) = 0.005386659641144906d0
            xpg(36) = 0.244349408253604301d0
            ypg(36) = 0.010318382922114762d0
            zpg(36) = 0.590349803755303313d0
            hpg(36) = 0.001512079937867526d0
            xpg(37) = 0.301650541545993227d0
            ypg(37) = 0.143001126409897532d0
            zpg(37) = 0.480719623973593007d0
            hpg(37) = 0.003820115194735675d0
            xpg(38) = 0.026205127676797405d0
            ypg(38) = 0.195686224271010268d0
            zpg(38) = 0.098161371735690487d0
            hpg(38) = 0.001711358660149380d0
            xpg(39) = 0.045499798757061205d0
            ypg(39) = 0.198825125666394756d0
            zpg(39) = 0.476864662926583267d0
            hpg(39) = 0.004286456909734536d0
            xpg(40) = 0.031438349393581904d0
            ypg(40) = 0.035932631783827859d0
            zpg(40) = 0.897709152837369758d0
            hpg(40) = 0.000654696883090330d0
            xpg(41) = 0.029122288726517093d0
            ypg(41) = 0.000575151627839405d0
            zpg(41) = 0.318843305523242203d0
            hpg(41) = 0.000564284111015018d0
            xpg(42) = 0.227222236394676046d0
            ypg(42) = 0.090205672555340825d0
            zpg(42) = 0.025161224041421096d0
            hpg(42) = 0.001647659837258332d0
            xpg(43) = 0.177369393271325314d0
            ypg(43) = 0.265942138332572336d0
            zpg(43) = 0.405367504954256102d0
            hpg(43) = 0.006021249539919202d0
            xpg(44) = 0.643148477427071416d0
            ypg(44) = 0.042827126064195619d0
            zpg(44) = 0.311367138822065673d0
            hpg(44) = 0.000630637466159878d0
            xpg(45) = 0.000265138775209478d0
            ypg(45) = 0.894349756040724680d0
            zpg(45) = 0.034217235095180989d0
            hpg(45) = 0.000371386056743092d0
            xpg(46) = 0.034361147868874548d0
            ypg(46) = 0.678005868883224272d0
            zpg(46) = 0.015142551369582735d0
            hpg(46) = 0.001005553900780592d0
            xpg(47) = 0.499987988264108349d0
            ypg(47) = 0.174203639875196185d0
            zpg(47) = 0.288600852957414206d0
            hpg(47) = 0.003869143307966543d0
            xpg(48) = 0.444450176577966701d0
            ypg(48) = 0.413505328487194281d0
            zpg(48) = 0.141928871152873136d0
            hpg(48) = 0.001203986873320328d0
            xpg(49) = 0.094338570241813714d0
            ypg(49) = 0.487240069581508865d0
            zpg(49) = 0.070813794150528460d0
            hpg(49) = 0.003612703301973658d0
            xpg(50) = 0.601187681268100956d0
            ypg(50) = 0.000168351201689568d0
            zpg(50) = 0.145441348608407450d0
            hpg(50) = 0.001086942693621432d0
            xpg(51) = 0.324968099558252989d0
            ypg(51) = 0.279627178020271888d0
            zpg(51) = 0.055728152910864672d0
            hpg(51) = 0.005386558789768295d0
            xpg(52) = 0.108923858645604586d0
            ypg(52) = 0.023167346067972153d0
            zpg(52) = 0.474948908362767983d0
            hpg(52) = 0.002405024565396667d0
            xpg(53) = 0.045372216771284115d0
            ypg(53) = 0.157701445237154310d0
            zpg(53) = 0.013406225503699345d0
            hpg(53) = 0.000903591390057696d0
            xpg(54) = 0.020627840314891339d0
            ypg(54) = 0.307195260583592497d0
            zpg(54) = 0.281509741910901848d0
            hpg(54) = 0.002257631630142710d0
            xpg(55) = 0.041661804502514335d0
            ypg(55) = 0.340292507411746880d0
            zpg(55) = 0.141944116896961888d0
            hpg(55) = 0.002189678767723374d0
            xpg(56) = 0.523224605694665751d0
            ypg(56) = 0.017193983506809406d0
            zpg(56) = 0.374786716728893110d0
            hpg(56) = 0.001714680538463014d0
            xpg(57) = 0.933382696541461643d0
            ypg(57) = 0.035566520008633332d0
            zpg(57) = 0.030717804976419337d0
            hpg(57) = 0.000221955508118992d0
            xpg(58) = 0.036879453168259424d0
            ypg(58) = 0.106873346342827945d0
            zpg(58) = 0.287827433681412943d0
            hpg(58) = 0.003345558259256166d0
            xpg(59) = 0.239069587380380989d0
            ypg(59) = 0.270775039711668146d0
            zpg(59) = 0.490031479403939609d0
            hpg(59) = 0.001415285833279145d0
            xpg(60) = 0.009877454016776290d0
            ypg(60) = 0.051594281605137175d0
            zpg(60) = 0.557651851334661132d0
            hpg(60) = 0.001171429165936078d0
            xpg(61) = 0.146872366907599138d0
            ypg(61) = 0.316319638644184553d0
            zpg(61) = 0.014096465583744309d0
            hpg(61) = 0.001971664228806746d0
            xpg(62) = 0.826009360652385629d0
            ypg(62) = 0.030616749579124643d0
            zpg(62) = 0.031463703602632070d0
            hpg(62) = 0.001085750735261891d0
            xpg(63) = 0.138011581734029365d0
            ypg(63) = 0.103133696653892385d0
            zpg(63) = 0.724684843262315566d0
            hpg(63) = 0.002199825850506972d0
            xpg(64) = 0.285235140396084042d0
            ypg(64) = 0.392704378842466758d0
            zpg(64) = 0.267778380570343726d0
            hpg(64) = 0.00467671077131169d0
            xpg(65) = 0.019480603322631985d0
            ypg(65) = 0.689164815745862530d0
            zpg(65) = 0.240662250891233125d0
            hpg(65) = 0.001327039440073762d0
            xpg(66) = 0.183247355053512386d0
            ypg(66) = 0.641098146347511995d0
            zpg(66) = 0.144071081394267842d0
            hpg(66) = 0.002476530267281989d0
            xpg(67) = 0.060549351184396778d0
            ypg(67) = 0.703599279613752695d0
            zpg(67) = 0.106403896783955549d0
            hpg(67) = 0.002498267628221633d0
            xpg(68) = 0.000591332247397936d0
            ypg(68) = 0.583419304082161433d0
            zpg(68) = 0.141904946321486920d0
            hpg(68) = 0.001140110042090766d0
            xpg(69) = 0.019738220220580314d0
            ypg(69) = 0.198339194163258624d0
            zpg(69) = 0.781488733462818138d0
            hpg(69) = 0.000361032571785597d0
            xpg(70) = 0.002021795121763506d0
            ypg(70) = 0.026560689677404987d0
            zpg(70) = 0.000256917470722789d0
            hpg(70) = 8.272224559376257d-05
            xpg(71) = 0.158367966418244876d0
            ypg(71) = 0.274068575612160929d0
            zpg(71) = 0.207880371717915191d0
            hpg(71) = 0.006074653578999147d0
            xpg(72) = 0.116460517332001445d0
            ypg(72) = 0.025537416473295411d0
            zpg(72) = 0.033545631729354729d0
            hpg(72) = 0.000931537843211918d0
            xpg(73) = 0.041433299417981430d0
            ypg(73) = 0.020555940628428632d0
            zpg(73) = 0.757109413373239706d0
            hpg(73) = 0.001156808227877411d0
            xpg(74) = 0.357232878847024485d0
            ypg(74) = 0.006146477733479475d0
            zpg(74) = 0.046577925570050455d0
            hpg(74) = 0.001023149007760126d0

        else if (fapg .eq. 'FPG94') then
! --------- Order 11 : High-order cubature rules for tetrahedra - 10.1002/nme.6313
            xpg(1) = 0.771184477823240208d0
            ypg(1) = 0.030149476635775709d0
            zpg(1) = 0.174079620317015995d0
            hpg(1) = 0.000925106164876252d0
            xpg(2) = 0.539439579554752608d0
            ypg(2) = 0.404452877552081177d0
            zpg(2) = 0.018874512805062669d0
            hpg(2) = 0.001024625535741073d0
            xpg(3) = 0.645757520465285451d0
            ypg(3) = 0.002468748755846532d0
            zpg(3) = 0.208388525680420460d0
            hpg(3) = 0.000859849262498925d0
            xpg(4) = 0.035206867250910605d0
            ypg(4) = 0.814674695357489335d0
            zpg(4) = 0.024980334075925634d0
            hpg(4) = 0.000660387079565831d0
            xpg(5) = 0.202143667520649522d0
            ypg(5) = 0.022980550715343346d0
            zpg(5) = 0.456996636702335749d0
            hpg(5) = 0.002040268151785345d0
            xpg(6) = 0.028709462926826182d0
            ypg(6) = 0.165595710768352266d0
            zpg(6) = 0.803273315480531127d0
            hpg(6) = 0.000347478072457070d0
            xpg(7) = 0.141931798748882889d0
            ypg(7) = 0.037474378133749342d0
            zpg(7) = 0.262394060027782822d0
            hpg(7) = 0.002501964533206585d0
            xpg(8) = 0.270122388743190010d0
            ypg(8) = 0.273271419323451026d0
            zpg(8) = 0.268925750693521426d0
            hpg(8) = 0.005192661065490372d0
            xpg(9) = 0.251018391068873245d0
            ypg(9) = 0.072094393510054992d0
            zpg(9) = 0.532678775541007763d0
            hpg(9) = 0.002030909349491641d0
            xpg(10) = 0.008015049556647654d0
            ypg(10) = 0.451735313172501196d0
            zpg(10) = 0.234948389704595852d0
            hpg(10) = 0.001347052252536158d0
            xpg(11) = 0.125611232214661904d0
            ypg(11) = 0.700825962749978613d0
            zpg(11) = 0.136562426608691700d0
            hpg(11) = 0.001754106014935642d0
            xpg(12) = 0.308299152461192647d0
            ypg(12) = 0.009474230997113106d0
            zpg(12) = 0.253901971967136447d0
            hpg(12) = 0.001165545974015295d0
            xpg(13) = 0.289450536294409354d0
            ypg(13) = 0.101287299579482284d0
            zpg(13) = 0.168094734924809082d0
            hpg(13) = 0.003569807153356659d0
            xpg(14) = 0.483149250392916685d0
            ypg(14) = 0.279296127839811641d0
            zpg(14) = 0.157830469624292334d0
            hpg(14) = 0.003308084812543372d0
            xpg(15) = 0.093926781513942231d0
            ypg(15) = 0.294975229465335345d0
            zpg(15) = 0.355331029405975530d0
            hpg(15) = 0.004415370609550707d0
            xpg(16) = 0.280881662159444558d0
            ypg(16) = 0.023518360349339504d0
            zpg(16) = 0.094788777825783315d0
            hpg(16) = 0.001537129235252944d0
            xpg(17) = 0.023439761203146644d0
            ypg(17) = 0.034535084403350958d0
            zpg(17) = 0.916134170617309608d0
            hpg(17) = 0.000332226847710768d0
            xpg(18) = 0.555645005196962878d0
            ypg(18) = 0.163650772131173611d0
            zpg(18) = 0.269604138503832510d0
            hpg(18) = 0.001361763968790416d0
            xpg(19) = 0.177322993507837634d0
            ypg(19) = 0.473784751619088624d0
            zpg(19) = 0.226422120325091236d0
            hpg(19) = 0.003650059064430291d0
            xpg(20) = 0.144175911492444186d0
            ypg(20) = 0.522067582444028574d0
            zpg(20) = 0.022572295722751675d0
            hpg(20) = 0.001995309775337051d0
            xpg(21) = 0.115211177549809887d0
            ypg(21) = 0.109012746648940875d0
            zpg(21) = 0.098725179093727562d0
            hpg(21) = 0.002438559137298117d0
            xpg(22) = 0.473317868651973181d0
            ypg(22) = 0.151714015847759890d0
            zpg(22) = 0.193813647753748833d0
            hpg(22) = 0.003618875285752321d0
            xpg(23) = 0.512222527188220078d0
            ypg(23) = 0.028719566427552803d0
            zpg(23) = 0.113591642215924660d0
            hpg(23) = 0.002165630219681029d0
            xpg(24) = 0.096204677945224707d0
            ypg(24) = 0.569886714540261817d0
            zpg(24) = 0.327391509621996107d0
            hpg(24) = 0.000920702991750362d0
            xpg(25) = 0.307657160391516008d0
            ypg(25) = 0.024816553006698880d0
            zpg(25) = 0.635581703244283174d0
            hpg(25) = 0.001075483256501068d0
            xpg(26) = 0.082480948888603683d0
            ypg(26) = 0.110375760758596518d0
            zpg(26) = 0.377603996814200783d0
            hpg(26) = 0.003129261448074559d0
            xpg(27) = 0.022312616178277880d0
            ypg(27) = 0.339637541774684139d0
            zpg(27) = 0.475148862347027000d0
            hpg(27) = 0.002099780833032553d0
            xpg(28) = 0.127194338733586801d0
            ypg(28) = 0.017097617470621552d0
            zpg(28) = 0.834306938660672796d0
            hpg(28) = 0.000422637062252629d0
            xpg(29) = 0.023580209723993884d0
            ypg(29) = 0.239954131595569925d0
            zpg(29) = 0.293686963185919185d0
            hpg(29) = 0.002191480262335345d0
            xpg(30) = 0.126132877485423499d0
            ypg(30) = 0.827729796998234878d0
            zpg(30) = 0.030366275652447153d0
            hpg(30) = 0.000515356515025525d0
            xpg(31) = 0.543386047040068635d0
            ypg(31) = 0.079018722164231409d0
            zpg(31) = 0.288965388291952754d0
            hpg(31) = 0.002630426986673036d0
            xpg(32) = 0.333098304330796073d0
            ypg(32) = 0.127022744562832199d0
            zpg(32) = 0.451595720885338791d0
            hpg(32) = 0.002961272789011408d0
            xpg(33) = 0.142868227350987772d0
            ypg(33) = 0.021174352303942007d0
            zpg(33) = 0.689763521380865357d0
            hpg(33) = 0.001515534819314825d0
            xpg(34) = 0.014293678102444836d0
            ypg(34) = 0.163745879306136074d0
            zpg(34) = 0.039500251493601511d0
            hpg(34) = 0.000719177690444211d0
            xpg(35) = 0.168257214981281685d0
            ypg(35) = 0.597540184230500027d0
            zpg(35) = 0.069705980169490571d0
            hpg(35) = 0.001499010803503053d0
            xpg(36) = 0.137196683065723231d0
            ypg(36) = 0.108501148082759486d0
            zpg(36) = 0.716471912756905051d0
            hpg(36) = 0.001965790601732406d0
            xpg(37) = 0.343489892636919368d0
            ypg(37) = 0.479058724871692859d0
            zpg(37) = 0.017255650079845479d0
            hpg(37) = 0.001695752122494603d0
            xpg(38) = 0.306691309067315600d0
            ypg(38) = 0.382289085879555073d0
            zpg(38) = 0.291666031120305485d0
            hpg(38) = 0.001717664871734031d0
            xpg(39) = 0.028807250464156035d0
            ypg(39) = 0.383312273233856527d0
            zpg(39) = 0.003557852602795268d0
            hpg(39) = 0.000528613492932897d0
            xpg(40) = 0.273446145397616084d0
            ypg(40) = 0.551327032694155171d0
            zpg(40) = 0.154324337815933376d0
            hpg(40) = 0.001508973654255826d0
            xpg(41) = 0.026001577959626183d0
            ypg(41) = 0.536365999276620396d0
            zpg(41) = 0.058101356509930102d0
            hpg(41) = 0.001574180379784949d0
            xpg(42) = 0.031310451133808527d0
            ypg(42) = 0.028501979442441433d0
            zpg(42) = 0.823648808085688970d0
            hpg(42) = 0.000671976908896233d0
            xpg(43) = 0.355061934545681610d0
            ypg(43) = 0.291707428657239658d0
            zpg(43) = 0.007813939611909787d0
            hpg(43) = 0.00143087116985064d0
            xpg(44) = 0.691605469266916012d0
            ypg(44) = 0.068197701691966879d0
            zpg(44) = 0.089906413442127108d0
            hpg(44) = 0.002557651000466822d0
            xpg(45) = 0.009405227228336286d0
            ypg(45) = 0.045746955017375170d0
            zpg(45) = 0.414099309081322681d0
            hpg(45) = 0.000812240352421740d0
            xpg(46) = 0.023255132453904748d0
            ypg(46) = 0.915019828059501526d0
            zpg(46) = 0.032806945745177203d0
            hpg(46) = 0.000356500974358299d0
            xpg(47) = 0.702440264961588700d0
            ypg(47) = 0.009157891497218531d0
            zpg(47) = 0.023111234039573174d0
            hpg(47) = 0.000494652625277812d0
            xpg(48) = 0.787983124705038746d0
            ypg(48) = 0.089763202417867738d0
            zpg(48) = 0.002915510329983943d0
            hpg(48) = 0.000586393480521284d0
            xpg(49) = 0.035053062960690433d0
            ypg(49) = 0.032893413844118061d0
            zpg(49) = 0.021959727940517558d0
            hpg(49) = 0.000425841683697435d0
            xpg(50) = 0.434530539486262845d0
            ypg(50) = 0.025292615209571488d0
            zpg(50) = 0.021065053356866531d0
            hpg(50) = 0.000875242950092945d0
            xpg(51) = 0.146876024522310893d0
            ypg(51) = 0.312993000146837229d0
            zpg(51) = 0.036325031737106768d0
            hpg(51) = 0.002793362898256485d0
            xpg(52) = 0.907654818128354503d0
            ypg(52) = 0.043836618078476730d0
            zpg(52) = 0.031432504376643333d0
            hpg(52) = 0.000387354910181323d0
            xpg(53) = 0.114922602677344669d0
            ypg(53) = 0.017372561595524409d0
            zpg(53) = 0.088193335056638151d0
            hpg(53) = 0.000855419574251798d0
            xpg(54) = 0.161942637306232740d0
            ypg(54) = 0.720116576093354988d0
            zpg(54) = 0.015926580976723308d0
            hpg(54) = 0.000984100484714994d0
            xpg(55) = 0.025233036380728764d0
            ypg(55) = 0.155149708646596402d0
            zpg(55) = 0.724273974425022178d0
            hpg(55) = 0.001606595643187324d0
            xpg(56) = 0.520052988544047412d0
            ypg(56) = 0.363665988935348739d0
            zpg(56) = 0.107961167137698022d0
            hpg(56) = 0.00105880131425778d0
            xpg(57) = 0.348598284107937657d0
            ypg(57) = 0.121467251205617964d0
            zpg(57) = 0.523443566060930415d0
            hpg(57) = 0.000971878458845587d0
            xpg(58) = 0.867729476024243266d0
            ypg(58) = 0.006252178054674931d0
            zpg(58) = 0.042409483427574491d0
            hpg(58) = 0.000357297834922962d0
            xpg(59) = 0.553760381384398693d0
            ypg(59) = 0.102776478420296773d0
            zpg(59) = 0.019537291657762027d0
            hpg(59) = 0.001619427095010102d0
            xpg(60) = 0.133582851361884819d0
            ypg(60) = 0.424904025908364372d0
            zpg(60) = 0.383608557846322173d0
            hpg(60) = 0.002538831217581424d0
            xpg(61) = 0.105431511160295359d0
            ypg(61) = 0.643534822154244835d0
            zpg(61) = 0.098706774763143021d0
            hpg(61) = 0.001500333540763712d0
            xpg(62) = 0.028117284133717343d0
            ypg(62) = 0.110890421021322858d0
            zpg(62) = 0.202061004868490930d0
            hpg(62) = 0.001811607813046539d0
            xpg(63) = 0.386947924866330985d0
            ypg(63) = 0.010909165218869565d0
            zpg(63) = 0.461826511210971549d0
            hpg(63) = 0.001272432187734095d0
            xpg(64) = 0.145259476356746304d0
            ypg(64) = 0.234911361430609445d0
            zpg(64) = 0.517425956285856767d0
            hpg(64) = 0.0037514820673859d0
            xpg(65) = 0.239521668999070826d0
            ypg(65) = 0.128354964775924320d0
            zpg(65) = 0.346340642065539261d0
            hpg(65) = 0.004695756987649625d0
            xpg(66) = 0.287790765267758187d0
            ypg(66) = 0.136771752113077852d0
            zpg(66) = 0.032575076367135517d0
            hpg(66) = 0.002589451082356837d0
            xpg(67) = 0.422238381839457519d0
            ypg(67) = 0.042917237666677621d0
            zpg(67) = 0.284973941983708864d0
            hpg(67) = 0.00283205792141347d0
            xpg(68) = 0.750772293893801715d0
            ypg(68) = 0.199552900709034594d0
            zpg(68) = 0.023912643073765792d0
            hpg(68) = 0.000751803261814584d0
            xpg(69) = 0.024298767767161547d0
            ypg(69) = 0.781292298193432658d0
            zpg(69) = 0.164364521007659802d0
            hpg(69) = 0.000852426235258589d0
            xpg(70) = 0.318909750901157812d0
            ypg(70) = 0.623796677577883756d0
            zpg(70) = 0.030969433913046913d0
            hpg(70) = 0.001142538097296711d0
            xpg(71) = 0.338391413024377089d0
            ypg(71) = 0.454050031755012838d0
            zpg(71) = 0.103285379536159045d0
            hpg(71) = 0.003137388928681668d0
            xpg(72) = 0.419393023646009756d0
            ypg(72) = 0.167853456304575201d0
            zpg(72) = 0.095826584901792632d0
            hpg(72) = 0.003360204988832838d0
            xpg(73) = 0.343125601943447928d0
            ypg(73) = 0.255238185521329699d0
            zpg(73) = 0.350746087406255896d0
            hpg(73) = 0.002600977037652207d0
            xpg(74) = 0.009312024647538243d0
            ypg(74) = 0.575720005798770144d0
            zpg(74) = 0.386868622874034281d0
            hpg(74) = 0.000561161471838098d0
            xpg(75) = 0.710551229203001273d0
            ypg(75) = 0.144349549150834578d0
            zpg(75) = 0.107742325195422503d0
            hpg(75) = 0.001698216049021417d0
            xpg(76) = 0.039007920817473344d0
            ypg(76) = 0.556871691287134417d0
            zpg(76) = 0.284649710486833009d0
            hpg(76) = 0.002551267017042547d0
            xpg(77) = 0.034413366386708231d0
            ypg(77) = 0.349065313744471406d0
            zpg(77) = 0.583784831929274290d0
            hpg(77) = 0.001492350871789632d0
            xpg(78) = 0.030482936103575377d0
            ypg(78) = 0.689354579857344228d0
            zpg(78) = 0.005360741363266648d0
            hpg(78) = 0.000402670519502796d0
            xpg(79) = 0.149359078771279288d0
            ypg(79) = 0.220344917465877018d0
            zpg(79) = 0.178549359594873500d0
            hpg(79) = 0.004884937106144391d0
            xpg(80) = 0.030716426468808937d0
            ypg(80) = 0.303979842649248436d0
            zpg(80) = 0.104100814728800395d0
            hpg(80) = 0.002325701632708601d0
            xpg(81) = 0.008511440408983542d0
            ypg(81) = 0.158711067239331793d0
            zpg(81) = 0.530954542707787487d0
            hpg(81) = 0.001207854443466630d0
            xpg(82) = 0.065399620983518019d0
            ypg(82) = 0.014497888533002320d0
            zpg(82) = 0.522226958738786577d0
            hpg(82) = 0.001039048579808548d0
            xpg(83) = 0.019822233385480147d0
            ypg(83) = 0.031932603045940630d0
            zpg(83) = 0.680658414281380335d0
            hpg(83) = 0.000798948890366117d0
            xpg(84) = 0.279865292460383972d0
            ypg(84) = 0.337317609456204603d0
            zpg(84) = 0.097495340016112871d0
            hpg(84) = 0.00483423181169034d0
            xpg(85) = 0.549764340110827705d0
            ypg(85) = 0.025929410713816946d0
            zpg(85) = 0.398878370569677286d0
            hpg(85) = 0.001056496344812016d0
            xpg(86) = 0.562025615609201907d0
            ypg(86) = 0.249436818062170008d0
            zpg(86) = 0.036028379709469166d0
            hpg(86) = 0.002942308804327285d0
            xpg(87) = 0.089216979965616418d0
            ypg(87) = 0.114638401437019911d0
            zpg(87) = 0.572182554059299974d0
            hpg(87) = 0.003184492496357301d0
            xpg(88) = 0.195549430066187139d0
            ypg(88) = 0.032545788012622124d0
            zpg(88) = 0.007923794324134563d0
            hpg(88) = 0.000550252409516352d0
            xpg(89) = 0.014164519920845215d0
            ypg(89) = 0.031370417332663679d0
            zpg(89) = 0.123305288904082602d0
            hpg(89) = 0.000468793528498282d0
            xpg(90) = 0.020019464138805672d0
            ypg(90) = 0.700121816352029338d0
            zpg(90) = 0.111653313014469755d0
            hpg(90) = 0.001309772262999681d0
            xpg(91) = 0.043832033736204037d0
            ypg(91) = 0.004485252771218634d0
            zpg(91) = 0.268197719995975189d0
            hpg(91) = 0.000530837438702655d0
            xpg(92) = 0.171529730438565469d0
            ypg(92) = 0.271807817871638560d0
            zpg(92) = 0.550761375264204434d0
            hpg(92) = 0.001091471802233006d0
            xpg(93) = 0.097898045170082096d0
            ypg(93) = 0.157052189386831146d0
            zpg(93) = 0.015266431471566375d0
            hpg(93) = 0.000964082274312674d0
            xpg(94) = 0.095203749770050551d0
            ypg(94) = 0.441026578450119737d0
            zpg(94) = 0.160565555733030674d0
            hpg(94) = 0.00417500003769135d0

        else if (fapg .eq. 'FPG117') then
! --------- Order 12 : High-order cubature rules for tetrahedra - 10.1002/nme.6313
            xpg(1) = 0.007573319150277204d0
            ypg(1) = 0.008102704337199682d0
            zpg(1) = 0.817847563349231734d0
            hpg(1) = 0.000160893581407839d0
            xpg(2) = 0.335639033905089418d0
            ypg(2) = 0.401866030805244767d0
            zpg(2) = 0.100805251809097971d0
            hpg(2) = 0.003502147895549203d0
            xpg(3) = 0.504687960523029788d0
            ypg(3) = 0.271791648582185284d0
            zpg(3) = 0.026615956617110240d0
            hpg(3) = 0.002177391770417587d0
            xpg(4) = 0.234348600377926096d0
            ypg(4) = 0.119172067273948856d0
            zpg(4) = 0.087484386570917726d0
            hpg(4) = 0.002506557962357166d0
            xpg(5) = 0.486780203021983850d0
            ypg(5) = 0.141002320608257690d0
            zpg(5) = 0.000418766475255139d0
            hpg(5) = 0.000440155948912292d0
            xpg(6) = 0.439288094894702550d0
            ypg(6) = 0.151894168909281691d0
            zpg(6) = 0.408567436062806424d0
            hpg(6) = 0.000631622825189968d0
            xpg(7) = 0.913146504567299230d0
            ypg(7) = 0.041353910033118560d0
            zpg(7) = 0.000632597009363113d0
            hpg(7) = 0.000175615901913196d0
            xpg(8) = 0.446310062042201321d0
            ypg(8) = 0.070818937855309348d0
            zpg(8) = 0.035706382390508774d0
            hpg(8) = 0.001418211236604260d0
            xpg(9) = 0.019473084124204713d0
            ypg(9) = 0.839177622280094508d0
            zpg(9) = 0.000228859006372356d0
            hpg(9) = 0.000168488110534687d0
            xpg(10) = 0.604617696699536115d0
            ypg(10) = 0.011942860906145131d0
            zpg(10) = 0.021615346047834469d0
            hpg(10) = 0.000499484240779006d0
            xpg(11) = 0.103751232146433087d0
            ypg(11) = 0.251812433620163382d0
            zpg(11) = 0.097355067510902810d0
            hpg(11) = 0.002765057516234731d0
            xpg(12) = 0.397163690579457339d0
            ypg(12) = 0.009721853696853559d0
            zpg(12) = 0.282569308879231102d0
            hpg(12) = 0.001251012621410417d0
            xpg(13) = 0.021529722487956813d0
            ypg(13) = 0.447816527463372264d0
            zpg(13) = 0.089893759113875664d0
            hpg(13) = 0.001588772686352459d0
            xpg(14) = 0.062690746676197662d0
            ypg(14) = 0.077552226371750557d0
            zpg(14) = 0.071311856466716265d0
            hpg(14) = 0.001028106125828816d0
            xpg(15) = 0.017317948635535781d0
            ypg(15) = 0.154195180853494569d0
            zpg(15) = 0.450466536767989240d0
            hpg(15) = 0.001449578789332851d0
            xpg(16) = 0.021863738731679730d0
            ypg(16) = 0.576844529923029415d0
            zpg(16) = 0.000691397701713393d0
            hpg(16) = 0.000305062014122628d0
            xpg(17) = 0.320235860865633659d0
            ypg(17) = 0.664556748863341406d0
            zpg(17) = 0.012883071166229717d0
            hpg(17) = 0.000206174084211390d0
            xpg(18) = 0.068680048860323333d0
            ypg(18) = 0.639089845949383071d0
            zpg(18) = 0.262857341422232948d0
            hpg(18) = 0.001470140631674588d0
            xpg(19) = 0.318884401570296605d0
            ypg(19) = 0.256327914882753353d0
            zpg(19) = 0.022054131594258476d0
            hpg(19) = 0.001738599322466227d0
            xpg(20) = 0.022589578976510324d0
            ypg(20) = 0.937720406047415612d0
            zpg(20) = 0.019743464689749439d0
            hpg(20) = 0.000163875621418107d0
            xpg(21) = 0.078311055152884541d0
            ypg(21) = 0.307470785843125467d0
            zpg(21) = 0.355095669290221356d0
            hpg(21) = 0.003490066544514464d0
            xpg(22) = 0.271645168850756597d0
            ypg(22) = 0.622189179053759944d0
            zpg(22) = 0.031408891841820921d0
            hpg(22) = 0.001660834003262054d0
            xpg(23) = 0.739707454238827899d0
            ypg(23) = 0.012272593459678340d0
            zpg(23) = 0.215741319147802413d0
            hpg(23) = 0.000466546262845791d0
            xpg(24) = 0.533996038737061268d0
            ypg(24) = 0.028145735078334812d0
            zpg(24) = 0.413294010593600842d0
            hpg(24) = 0.000989956354613176d0
            xpg(25) = 0.000682424254970806d0
            ypg(25) = 0.113764029274711770d0
            zpg(25) = 0.716671532362942034d0
            hpg(25) = 0.000510237424289878d0
            xpg(26) = 0.570191562660287613d0
            ypg(26) = 0.030660588306956372d0
            zpg(26) = 0.124688535176493010d0
            hpg(26) = 0.001933852192034303d0
            xpg(27) = 0.029959666144410446d0
            ypg(27) = 0.246449993220771771d0
            zpg(27) = 0.235683795737423030d0
            hpg(27) = 0.002256869685188125d0
            xpg(28) = 0.272722335362774767d0
            ypg(28) = 0.104423576137331206d0
            zpg(28) = 0.611192997818055432d0
            hpg(28) = 0.000712419665750717d0
            xpg(29) = 0.252804749851328283d0
            ypg(29) = 0.249825956198123443d0
            zpg(29) = 0.248834616285346028d0
            hpg(29) = 0.004872776209705338d0
            xpg(30) = 0.649779843839125881d0
            ypg(30) = 0.096228187859255542d0
            zpg(30) = 0.026984371825727670d0
            hpg(30) = 0.001360362003306534d0
            xpg(31) = 0.056994734232269696d0
            ypg(31) = 0.436057367698598275d0
            zpg(31) = 0.487080703163378818d0
            hpg(31) = 0.001157363673656273d0
            xpg(32) = 0.147008073183899589d0
            ypg(32) = 0.109105306642488351d0
            zpg(32) = 0.681588965314910176d0
            hpg(32) = 0.002127065824961539d0
            xpg(33) = 0.140754013402578585d0
            ypg(33) = 0.007503496946606807d0
            zpg(33) = 0.789315652838747194d0
            hpg(33) = 0.000490104001848109d0
            xpg(34) = 0.246742414843959916d0
            ypg(34) = 0.125150502094054022d0
            zpg(34) = 0.016170640241171650d0
            hpg(34) = 0.001041230341938821d0
            xpg(35) = 0.770269827307140079d0
            ypg(35) = 0.000587377937740670d0
            zpg(35) = 0.093622757280270978d0
            hpg(35) = 0.000350736108291456d0
            xpg(36) = 0.222033328885992751d0
            ypg(36) = 0.064522241397773091d0
            zpg(36) = 0.406709813020482124d0
            hpg(36) = 0.003192469921872888d0
            xpg(37) = 0.023164763903695016d0
            ypg(37) = 0.678197045349286112d0
            zpg(37) = 0.055252495037972575d0
            hpg(37) = 0.001119233439143887d0
            xpg(38) = 0.115381162167599145d0
            ypg(38) = 0.698584529376278651d0
            zpg(38) = 0.018225621329625292d0
            hpg(38) = 0.001145124067972906d0
            xpg(39) = 0.346036261952086300d0
            ypg(39) = 0.020042645135158001d0
            zpg(39) = 0.013112646874512993d0
            hpg(39) = 0.000467336597986420d0
            xpg(40) = 0.088257612665764397d0
            ypg(40) = 0.221890743194826060d0
            zpg(40) = 0.683472856377851737d0
            hpg(40) = 0.000673109011605284d0
            xpg(41) = 0.023284151841602800d0
            ypg(41) = 0.113187653254642652d0
            zpg(41) = 0.835251185930520350d0
            hpg(41) = 0.000589098418756660d0
            xpg(42) = 0.634690392213963504d0
            ypg(42) = 0.217644584420059457d0
            zpg(42) = 0.131249760598345170d0
            hpg(42) = 0.001024975822834534d0
            xpg(43) = 0.745454973014218010d0
            ypg(43) = 0.064598031133060777d0
            zpg(43) = 0.189497970015329175d0
            hpg(43) = 0.000345322966457884d0
            xpg(44) = 0.418266580660706465d0
            ypg(44) = 0.160702961951829284d0
            zpg(44) = 0.096503370901911198d0
            hpg(44) = 0.003278433182654423d0
            xpg(45) = 0.417908686866454497d0
            ypg(45) = 0.406757465470303153d0
            zpg(45) = 0.132321097887058639d0
            hpg(45) = 0.001902208936896543d0
            xpg(46) = 0.197444446625025359d0
            ypg(46) = 0.472071220352845518d0
            zpg(46) = 0.314343068928740679d0
            hpg(46) = 0.001536254193504916d0
            xpg(47) = 0.094166584616699269d0
            ypg(47) = 0.157667410765404058d0
            zpg(47) = 0.021577257901887025d0
            hpg(47) = 0.000978174204766414d0
            xpg(48) = 0.284876489754350499d0
            ypg(48) = 0.467196438356977584d0
            zpg(48) = 0.013722740935102271d0
            hpg(48) = 0.001575032336884164d0
            xpg(49) = 0.794178780590717955d0
            ypg(49) = 0.026170877920180374d0
            zpg(49) = 0.021701086045975721d0
            hpg(49) = 0.000526662234984192d0
            xpg(50) = 0.195590291625625201d0
            ypg(50) = 0.026002474109005220d0
            zpg(50) = 0.235501302560273368d0
            hpg(50) = 0.002099879003872708d0
            xpg(51) = 0.013155058602447971d0
            ypg(51) = 0.199518651418892426d0
            zpg(51) = 0.109327556584519913d0
            hpg(51) = 0.000971833542260749d0
            xpg(52) = 0.011632083588120143d0
            ypg(52) = 0.144660921622416180d0
            zpg(52) = 0.016634291609880988d0
            hpg(52) = 0.000273656874296067d0
            xpg(53) = 0.000771066784710497d0
            ypg(53) = 0.410565942866637932d0
            zpg(53) = 0.274975395518310338d0
            hpg(53) = 0.000864901959574101d0
            xpg(54) = 0.192590984714032432d0
            ypg(54) = 0.187945418012733421d0
            zpg(54) = 0.456322016615619711d0
            hpg(54) = 0.00420302888905984d0
            xpg(55) = 0.045909577413876234d0
            ypg(55) = 0.221113817593067456d0
            zpg(55) = 0.638528455018590329d0
            hpg(55) = 0.001839782770111300d0
            xpg(56) = 0.220971693343654317d0
            ypg(56) = 0.265633128660484985d0
            zpg(56) = 0.486529127435705494d0
            hpg(56) = 0.002116832819819149d0
            xpg(57) = 0.357470777869401210d0
            ypg(57) = 0.277141801979810569d0
            zpg(57) = 0.294759700258399266d0
            hpg(57) = 0.003741027203502297d0
            xpg(58) = 0.183466194257216202d0
            ypg(58) = 0.435430944306639295d0
            zpg(58) = 0.259614795528458698d0
            hpg(58) = 0.004154118162159714d0
            xpg(59) = 0.778678089436258600d0
            ypg(59) = 0.068245562439167864d0
            zpg(59) = 0.090290643840915190d0
            hpg(59) = 0.001438652113813070d0
            xpg(60) = 0.083503684954383061d0
            ypg(60) = 0.013028074713012028d0
            zpg(60) = 0.648107440400148949d0
            hpg(60) = 0.000880133303839030d0
            xpg(61) = 0.411181164699107537d0
            ypg(61) = 0.015358906074005189d0
            zpg(61) = 0.438262647584127045d0
            hpg(61) = 0.001304118812907386d0
            xpg(62) = 0.662659511952267236d0
            ypg(62) = 0.307848134174069460d0
            zpg(62) = 0.018360366146401999d0
            hpg(62) = 0.000372015864893222d0
            xpg(63) = 0.136357487041202633d0
            ypg(63) = 0.031193581339513479d0
            zpg(63) = 0.000196913645331014d0
            hpg(63) = 0.000261354708760932d0
            xpg(64) = 0.031890170593405439d0
            ypg(64) = 0.841530991646365015d0
            zpg(64) = 0.126369124232186355d0
            hpg(64) = 0.000259610895172115d0
            xpg(65) = 0.099651058826745358d0
            ypg(65) = 0.702090188467081304d0
            zpg(65) = 0.122127435345434447d0
            hpg(65) = 0.001877377428352171d0
            xpg(66) = 0.401276754723172710d0
            ypg(66) = 0.100136091003559254d0
            zpg(66) = 0.299151304841005946d0
            hpg(66) = 0.003736953199721342d0
            xpg(67) = 0.111034842723138529d0
            ypg(67) = 0.144266851566993967d0
            zpg(67) = 0.343031822242500720d0
            hpg(67) = 0.003083568160866131d0
            xpg(68) = 0.113280336559172134d0
            ypg(68) = 0.498852395172576347d0
            zpg(68) = 0.033103213917688019d0
            hpg(68) = 0.00203747254993593d0
            xpg(69) = 0.185125128533730573d0
            ypg(69) = 0.025257925676531014d0
            zpg(69) = 0.071959988989683937d0
            hpg(69) = 0.001245138355901754d0
            xpg(70) = 0.168118084191225036d0
            ypg(70) = 0.733419717861681348d0
            zpg(70) = 0.097073594275626304d0
            hpg(70) = 0.000423623226304132d0
            xpg(71) = 0.032848323900184274d0
            ypg(71) = 0.587282168449124307d0
            zpg(71) = 0.175258410473463196d0
            hpg(71) = 0.001856759575857452d0
            xpg(72) = 0.019024889324311123d0
            ypg(72) = 0.279840161794067478d0
            zpg(72) = 0.494125547771558400d0
            hpg(72) = 0.001381182369839022d0
            xpg(73) = 0.221903197223332774d0
            ypg(73) = 0.033367086971297833d0
            zpg(73) = 0.588974819318790321d0
            hpg(73) = 0.002059698489709578d0
            xpg(74) = 0.122554881882387404d0
            ypg(74) = 0.042215898336318085d0
            zpg(74) = 0.833274383322793019d0
            hpg(74) = 0.000318737904171449d0
            xpg(75) = 0.337345868730355053d0
            ypg(75) = 0.079566853310368323d0
            zpg(75) = 0.197500819724511054d0
            hpg(75) = 0.002866129834801582d0
            xpg(76) = 0.000428115271007478d0
            ypg(76) = 0.733945599399890392d0
            zpg(76) = 0.194863559194742178d0
            hpg(76) = 0.000443567902085300d0
            xpg(77) = 0.004161074458799275d0
            ypg(77) = 0.303222184937155532d0
            zpg(77) = 0.666802824462311695d0
            hpg(77) = 0.000351975452532082d0
            xpg(78) = 0.056070459013322461d0
            ypg(78) = 0.005594038297952845d0
            zpg(78) = 0.222340652984313376d0
            hpg(78) = 0.000538513567447074d0
            xpg(79) = 0.023133144513095003d0
            ypg(79) = 0.017752321317989753d0
            zpg(79) = 0.930721535322083907d0
            hpg(79) = 0.000214542442085025d0
            xpg(80) = 0.916092917700588663d0
            ypg(80) = 0.009566613286067694d0
            zpg(80) = 0.055095284958799263d0
            hpg(80) = 0.000199867864688934d0
            xpg(81) = 0.376554224393516091d0
            ypg(81) = 0.015641600921766532d0
            zpg(81) = 0.106988587253400618d0
            hpg(81) = 0.001219967136328668d0
            xpg(82) = 0.108846319536463861d0
            ypg(82) = 0.337519368118120109d0
            zpg(82) = 0.466990736443575035d0
            hpg(82) = 0.002341134863577650d0
            xpg(83) = 0.012737478002952097d0
            ypg(83) = 0.033809784867966477d0
            zpg(83) = 0.596041903282640726d0
            hpg(83) = 0.000659390957044575d0
            xpg(84) = 0.016628241186543641d0
            ypg(84) = 0.080165300158286856d0
            zpg(84) = 0.308670402308715901d0
            hpg(84) = 0.001067235672475957d0
            xpg(85) = 0.028092160358597578d0
            ypg(85) = 0.333976931830175418d0
            zpg(85) = 0.021862188638790174d0
            hpg(85) = 0.000834942651839463d0
            xpg(86) = 0.419130889959823830d0
            ypg(86) = 0.326709036151804164d0
            zpg(86) = 0.253411625169333294d0
            hpg(86) = 0.000734534192805942d0
            xpg(87) = 0.004983324642759546d0
            ypg(87) = 0.446842916140259186d0
            zpg(87) = 0.454715932790208060d0
            hpg(87) = 0.000622654558374897d0
            xpg(88) = 0.059590516552103423d0
            ypg(88) = 0.001465584001853574d0
            zpg(88) = 0.078142972866940428d0
            hpg(88) = 0.000212595738407722d0
            xpg(89) = 0.076229019036361516d0
            ypg(89) = 0.107911620286574939d0
            zpg(89) = 0.578833789104501775d0
            hpg(89) = 0.002760506457723844d0
            xpg(90) = 0.567591351216452310d0
            ypg(90) = 0.127124794601758273d0
            zpg(90) = 0.258032475061934613d0
            hpg(90) = 0.002406324845892238d0
            xpg(91) = 0.205251536375455909d0
            ypg(91) = 0.000632138723830099d0
            zpg(91) = 0.444574604850063003d0
            hpg(91) = 0.000603006125143326d0
            xpg(92) = 0.023348156129191414d0
            ypg(92) = 0.028327999446992159d0
            zpg(92) = 0.022620449825243579d0
            hpg(92) = 0.000240839850977726d0
            xpg(93) = 0.239670663790656278d0
            ypg(93) = 0.316488151078752905d0
            zpg(93) = 0.095843348019856314d0
            hpg(93) = 0.003479439478305232d0
            xpg(94) = 0.840668066350032022d0
            ypg(94) = 0.119040873507988227d0
            zpg(94) = 0.037160060238494001d0
            hpg(94) = 0.000300241190103332d0
            xpg(95) = 0.026226046466612823d0
            ypg(95) = 0.830181622842096733d0
            zpg(95) = 0.064433041051417722d0
            hpg(95) = 0.000725817895092358d0
            xpg(96) = 0.729712082501644023d0
            ypg(96) = 0.174623213665996243d0
            zpg(96) = 0.017613802103044813d0
            hpg(96) = 0.000993765605639344d0
            xpg(97) = 0.455272251501793014d0
            ypg(97) = 0.487856163325192272d0
            zpg(97) = 0.056462617770123451d0
            hpg(97) = 0.000416167243546169d0
            xpg(98) = 0.122877885888670430d0
            ypg(98) = 0.823304084967301288d0
            zpg(98) = 0.021134275184768248d0
            hpg(98) = 0.000628231918166048d0
            xpg(99) = 0.355676065200314857d0
            ypg(99) = 0.110152329140885843d0
            zpg(99) = 0.458505151620200847d0
            hpg(99) = 0.00283042214225853d0
            xpg(100) = 0.148508639530215077d0
            ypg(100) = 0.320360496974352051d0
            zpg(100) = 0.015152988695219824d0
            hpg(100) = 0.001186631235326901d0
            xpg(101) = 0.255183769645292993d0
            ypg(101) = 0.565918718121523920d0
            zpg(101) = 0.143416258879199085d0
            hpg(101) = 0.001914745513488388d0
            xpg(102) = 0.000301733246281871d0
            ypg(102) = 0.616901402631183677d0
            zpg(102) = 0.370585975390405622d0
            hpg(102) = 0.000189742898209514d0
            xpg(103) = 0.604811935980153497d0
            ypg(103) = 0.118167075844067646d0
            zpg(103) = 0.121231848536972998d0
            hpg(103) = 0.002159556523590855d0
            xpg(104) = 0.181235764899320168d0
            ypg(104) = 0.167288068367082484d0
            zpg(104) = 0.207540385263766456d0
            hpg(104) = 0.002993131608156337d0
            xpg(105) = 0.455659877284173754d0
            ypg(105) = 0.223448777664227983d0
            zpg(105) = 0.164400198109932657d0
            hpg(105) = 0.002956738972987333d0
            xpg(106) = 0.316808274763732926d0
            ypg(106) = 0.020391700648517711d0
            zpg(106) = 0.634412113273361319d0
            hpg(106) = 0.000798755419948512d0
            xpg(107) = 0.588572587815202787d0
            ypg(107) = 0.262519550720361131d0
            zpg(107) = 0.084212335374154119d0
            hpg(107) = 0.001619161127969705d0
            xpg(108) = 0.079757088127363543d0
            ypg(108) = 0.086923654159650584d0
            zpg(108) = 0.185712282784471180d0
            hpg(108) = 0.002144159961322959d0
            xpg(109) = 0.045821694743356120d0
            ypg(109) = 0.508728056330881613d0
            zpg(109) = 0.318863637071560957d0
            hpg(109) = 0.002086565345722996d0
            xpg(110) = 0.096491345345384113d0
            ypg(110) = 0.396892709045254671d0
            zpg(110) = 0.177482917714969264d0
            hpg(110) = 0.003582326366862225d0
            xpg(111) = 0.051817353834661560d0
            ypg(111) = 0.045668299129327605d0
            zpg(111) = 0.788572097154935747d0
            hpg(111) = 0.001018437803279553d0
            xpg(112) = 0.157483283913500632d0
            ypg(112) = 0.547958051535508966d0
            zpg(112) = 0.106128579676776613d0
            hpg(112) = 0.002901805502947609d0
            xpg(113) = 0.599792349445298793d0
            ypg(113) = 0.029575905483807844d0
            zpg(113) = 0.245154027848701984d0
            hpg(113) = 0.001722747636558756d0
            xpg(114) = 0.075887324796519839d0
            ypg(114) = 0.031465327566215261d0
            zpg(114) = 0.422506306539806616d0
            hpg(114) = 0.001784162545407309d0
            xpg(115) = 0.481394060019316360d0
            ypg(115) = 0.433055832551902531d0
            zpg(115) = 0.018100651663879534d0
            hpg(115) = 0.001217403964030057d0
            xpg(116) = 0.000357423742141919d0
            ypg(116) = 0.043267314223620196d0
            zpg(116) = 0.141139029348213005d0
            hpg(116) = 0.000320815536278386d0
            xpg(117) = 0.005493262085055545d0
            ypg(117) = 0.002956537533261605d0
            zpg(117) = 0.399678978713579598d0
            hpg(117) = 0.000155700413290461d0

        else if (fapg .eq. 'FPG144') then
! --------- Order 13 : High-order cubature rules for tetrahedra - 10.1002/nme.6313
            xpg(1) = 0.402984982233501364d0
            ypg(1) = 0.055921339432015265d0
            zpg(1) = 0.271219477793538688d0
            hpg(1) = 0.002289248926936554d0
            xpg(2) = 0.004853146143288019d0
            ypg(2) = 0.387033650890423571d0
            zpg(2) = 0.026728360846843798d0
            hpg(2) = 0.000303784992521207d0
            xpg(3) = 0.753957633331647895d0
            ypg(3) = 0.122065601359220950d0
            zpg(3) = 0.005667902372503289d0
            hpg(3) = 0.000476106605963919d0
            xpg(4) = 0.314715187476880793d0
            ypg(4) = 0.023663950466977122d0
            zpg(4) = 0.445700133372992643d0
            hpg(4) = 0.001586536862032418d0
            xpg(5) = 0.654130572185355956d0
            ypg(5) = 0.036524824827036006d0
            zpg(5) = 0.238894607773945171d0
            hpg(5) = 0.001124759477689653d0
            xpg(6) = 0.925310762159826580d0
            ypg(6) = 0.024408337832415457d0
            zpg(6) = 0.023543214558058812d0
            hpg(6) = 0.000254741193220830d0
            xpg(7) = 0.003814050373789552d0
            ypg(7) = 0.419992480052590582d0
            zpg(7) = 0.491165374858810625d0
            hpg(7) = 0.000430389125673390d0
            xpg(8) = 0.127541682781534031d0
            ypg(8) = 0.838898505728828530d0
            zpg(8) = 0.026073037335438894d0
            hpg(8) = 0.000247946025209262d0
            xpg(9) = 0.243365700459289892d0
            ypg(9) = 0.083315832796637217d0
            zpg(9) = 0.604554802435564681d0
            hpg(9) = 0.001543178118220842d0
            xpg(10) = 0.236583611940113317d0
            ypg(10) = 0.043971703620042066d0
            zpg(10) = 0.235806752472383750d0
            hpg(10) = 0.002187510358845000d0
            xpg(11) = 0.099085758499841081d0
            ypg(11) = 0.087803878579364372d0
            zpg(11) = 0.217788330978253766d0
            hpg(11) = 0.001737414431223952d0
            xpg(12) = 0.022209171145425297d0
            ypg(12) = 0.083514758866084460d0
            zpg(12) = 0.078059360456318597d0
            hpg(12) = 0.000689722817811123d0
            xpg(13) = 0.027735802513771241d0
            ypg(13) = 0.240623708560142147d0
            zpg(13) = 0.260299544171853057d0
            hpg(13) = 0.001698322036919328d0
            xpg(14) = 0.026064040858455498d0
            ypg(14) = 0.617069248495987800d0
            zpg(14) = 0.122312993867153189d0
            hpg(14) = 0.001351011404283771d0
            xpg(15) = 0.548377613350833787d0
            ypg(15) = 0.009975161650579685d0
            zpg(15) = 0.225906472332466421d0
            hpg(15) = 0.000567353009670014d0
            xpg(16) = 0.255185712901532492d0
            ypg(16) = 0.632327949540898999d0
            zpg(16) = 0.094499158004900266d0
            hpg(16) = 0.000928936659616203d0
            xpg(17) = 0.060770737098752822d0
            ypg(17) = 0.104187186041621738d0
            zpg(17) = 0.583190384961391733d0
            hpg(17) = 0.001850327528118551d0
            xpg(18) = 0.006887091028037954d0
            ypg(18) = 0.649454777541901479d0
            zpg(18) = 0.028372753366617474d0
            hpg(18) = 0.000350437621592761d0
            xpg(19) = 0.092377515357307454d0
            ypg(19) = 0.566688252212853680d0
            zpg(19) = 0.328477269355757370d0
            hpg(19) = 0.000543661943084833d0
            xpg(20) = 0.103822057768969504d0
            ypg(20) = 0.758692729983167073d0
            zpg(20) = 0.117267932583013573d0
            hpg(20) = 0.000643439256081267d0
            xpg(21) = 0.819381801515008975d0
            ypg(21) = 0.022292662579453507d0
            zpg(21) = 0.021693443903744023d0
            hpg(21) = 0.000464348182741589d0
            xpg(22) = 0.529009305528457029d0
            ypg(22) = 0.174358315803303489d0
            zpg(22) = 0.199666251728784034d0
            hpg(22) = 0.002480773583147739d0
            xpg(23) = 0.293458629693931953d0
            ypg(23) = 0.350070116142154106d0
            zpg(23) = 0.001795274011579795d0
            hpg(23) = 0.000575288018256703d0
            xpg(24) = 0.564988214514982947d0
            ypg(24) = 0.268592829213224079d0
            zpg(24) = 0.143785582764165247d0
            hpg(24) = 0.001234506052990718d0
            xpg(25) = 0.258184009376851105d0
            ypg(25) = 0.013914813987781348d0
            zpg(25) = 0.016018430803131482d0
            hpg(25) = 0.000268814367123680d0
            xpg(26) = 0.733339108976289046d0
            ypg(26) = 0.118940960177329500d0
            zpg(26) = 0.128691081422400734d0
            hpg(26) = 0.000798052840225632d0
            xpg(27) = 0.014881549172686446d0
            ypg(27) = 0.208087595329487181d0
            zpg(27) = 0.561560628877974740d0
            hpg(27) = 0.000948007785729163d0
            xpg(28) = 0.656344836524241954d0
            ypg(28) = 0.018962970895684189d0
            zpg(28) = 0.109046009383029542d0
            hpg(28) = 0.000910548937209906d0
            xpg(29) = 0.020370770140704450d0
            ypg(29) = 0.821646083133138127d0
            zpg(29) = 0.025123040496851188d0
            hpg(29) = 0.000477190012473442d0
            xpg(30) = 0.251926573322684032d0
            ypg(30) = 0.117802106175616010d0
            zpg(30) = 0.098248541159628664d0
            hpg(30) = 0.002445034614002043d0
            xpg(31) = 0.180162775566469209d0
            ypg(31) = 0.254898675164947931d0
            zpg(31) = 0.034323550806945220d0
            hpg(31) = 0.001498474097960685d0
            xpg(32) = 0.117893134955565719d0
            ypg(32) = 0.005805106374632742d0
            zpg(32) = 0.271159688521072745d0
            hpg(32) = 0.000640844323280677d0
            xpg(33) = 0.289341746012952299d0
            ypg(33) = 0.673349414317319746d0
            zpg(33) = 0.013232156981078750d0
            hpg(33) = 0.000412937028767928d0
            xpg(34) = 0.099365533978727170d0
            ypg(34) = 0.022176555764886803d0
            zpg(34) = 0.868322670887145857d0
            hpg(34) = 0.000194684467286981d0
            xpg(35) = 0.093336230918615892d0
            ypg(35) = 0.718307463015777390d0
            zpg(35) = 0.094078031429025482d0
            hpg(35) = 0.001311801844869701d0
            xpg(36) = 0.426758116514513269d0
            ypg(36) = 0.074086566225422657d0
            zpg(36) = 0.021516672509995639d0
            hpg(36) = 0.000988561230910533d0
            xpg(37) = 0.384139598789182279d0
            ypg(37) = 0.201876853314553156d0
            zpg(37) = 0.375361616838647587d0
            hpg(37) = 0.001731769459096245d0
            xpg(38) = 0.476570602766622508d0
            ypg(38) = 0.014787746748208942d0
            zpg(38) = 0.105787439748236554d0
            hpg(38) = 0.000796685409567047d0
            xpg(39) = 0.321365098957279913d0
            ypg(39) = 0.176343635108770512d0
            zpg(39) = 0.323033602262217536d0
            hpg(39) = 0.003771688424302016d0
            xpg(40) = 0.102435062173092363d0
            ypg(40) = 0.113434584355445556d0
            zpg(40) = 0.020123872549235117d0
            hpg(40) = 0.000815389770186190d0
            xpg(41) = 0.135037756242416323d0
            ypg(41) = 0.022016316254346787d0
            zpg(41) = 0.632636797647245032d0
            hpg(41) = 0.001292374652157503d0
            xpg(42) = 0.427155942740121640d0
            ypg(42) = 0.075575862070078761d0
            zpg(42) = 0.400637335226278795d0
            hpg(42) = 0.002394342611920835d0
            xpg(43) = 0.007423352058992606d0
            ypg(43) = 0.505580018305893625d0
            zpg(43) = 0.281421382627391863d0
            hpg(43) = 0.000793193418429923d0
            xpg(44) = 0.841573507009603731d0
            ypg(44) = 0.026808955765909372d0
            zpg(44) = 0.122772413631337963d0
            hpg(44) = 0.000268126474085205d0
            xpg(45) = 0.133314107311733479d0
            ypg(45) = 0.618498202346326114d0
            zpg(45) = 0.219403799104052771d0
            hpg(45) = 0.000991260778377488d0
            xpg(46) = 0.646377375910297796d0
            ypg(46) = 0.021588288299128453d0
            zpg(46) = 0.021836959259042318d0
            hpg(46) = 0.000590590799621874d0
            xpg(47) = 0.668258479817918975d0
            ypg(47) = 0.018169858602626334d0
            zpg(47) = 0.300141838725955202d0
            hpg(47) = 0.000357078802410538d0
            xpg(48) = 0.358934851555677009d0
            ypg(48) = 0.328788775603153552d0
            zpg(48) = 0.234612826768038255d0
            hpg(48) = 0.002732051028333130d0
            xpg(49) = 0.026032852075347600d0
            ypg(49) = 0.015385379193766444d0
            zpg(49) = 0.834470603891131600d0
            hpg(49) = 0.000346759171262602d0
            xpg(50) = 0.392235288840848633d0
            ypg(50) = 0.197542827301242210d0
            zpg(50) = 0.032372169951113028d0
            hpg(50) = 0.001934242916733430d0
            xpg(51) = 0.041110269018704013d0
            ypg(51) = 0.463559879573053404d0
            zpg(51) = 0.074429729567821355d0
            hpg(51) = 0.001524247597881742d0
            xpg(52) = 0.642959234195544140d0
            ypg(52) = 0.306996351051980014d0
            zpg(52) = 0.003031555307071350d0
            hpg(52) = 0.000352273630428293d0
            xpg(53) = 0.023711289540394213d0
            ypg(53) = 0.019933723900515017d0
            zpg(53) = 0.168691966307334958d0
            hpg(53) = 0.000445890338174757d0
            xpg(54) = 0.021353035315088654d0
            ypg(54) = 0.024070750964263175d0
            zpg(54) = 0.328664926728553080d0
            hpg(54) = 0.000566849623873895d0
            xpg(55) = 0.024335203528370929d0
            ypg(55) = 0.132605718115887142d0
            zpg(55) = 0.821734574996847754d0
            hpg(55) = 0.000486074420589401d0
            xpg(56) = 0.250084820890875814d0
            ypg(56) = 0.280993040860044879d0
            zpg(56) = 0.096281055906287962d0
            hpg(56) = 0.002814138418032146d0
            xpg(57) = 0.116133650848346780d0
            ypg(57) = 0.790645117642715923d0
            zpg(57) = 0.016733855340768187d0
            hpg(57) = 0.000630597678354571d0
            xpg(58) = 0.277358194171600939d0
            ypg(58) = 0.016340876406998575d0
            zpg(58) = 0.105688039922541745d0
            hpg(58) = 0.001012708394937109d0
            xpg(59) = 0.144722979054398163d0
            ypg(59) = 0.197633525126866042d0
            zpg(59) = 0.238233833340414637d0
            hpg(59) = 0.003278383176463416d0
            xpg(60) = 0.021479636783768247d0
            ypg(60) = 0.865880812317971689d0
            zpg(60) = 0.091522453183015660d0
            hpg(60) = 0.000303726011329897d0
            xpg(61) = 0.041352307796126416d0
            ypg(61) = 0.507312794765697479d0
            zpg(61) = 0.004708987583563917d0
            hpg(61) = 0.000323189494346036d0
            xpg(62) = 0.125184681712264324d0
            ypg(62) = 0.102142976871279299d0
            zpg(62) = 0.753658261075793100d0
            hpg(62) = 0.000712838205707703d0
            xpg(63) = 0.198887264727459630d0
            ypg(63) = 0.109571766120728380d0
            zpg(63) = 0.518054512015828644d0
            hpg(63) = 0.002401628518231711d0
            xpg(64) = 0.087145851931892808d0
            ypg(64) = 0.079248918765838934d0
            zpg(64) = 0.735543693350126650d0
            hpg(64) = 0.001252625769653002d0
            xpg(65) = 0.067673227190050811d0
            ypg(65) = 0.407168622017727006d0
            zpg(65) = 0.231466262379970431d0
            hpg(65) = 0.002942770244053756d0
            xpg(66) = 0.557240655968317913d0
            ypg(66) = 0.079314622683004829d0
            zpg(66) = 0.193907963885415261d0
            hpg(66) = 0.001837346664473145d0
            xpg(67) = 0.620402904892843208d0
            ypg(67) = 0.226205978429814246d0
            zpg(67) = 0.058553036363721456d0
            hpg(67) = 0.001791939733516674d0
            xpg(68) = 0.002398627388752580d0
            ypg(68) = 0.103575480307320558d0
            zpg(68) = 0.438645411804342848d0
            hpg(68) = 0.000528451297332165d0
            xpg(69) = 0.030133970220982284d0
            ypg(69) = 0.015999745174342064d0
            zpg(69) = 0.688374237970464162d0
            hpg(69) = 0.000458244103170407d0
            xpg(70) = 0.192774902327159691d0
            ypg(70) = 0.085630862627857275d0
            zpg(70) = 0.391167823978110883d0
            hpg(70) = 0.002856126433811271d0
            xpg(71) = 0.054074457199930708d0
            ypg(71) = 0.143373501032070525d0
            zpg(71) = 0.406265298123187443d0
            hpg(71) = 0.00171961331991048d0
            xpg(72) = 0.117358128684145899d0
            ypg(72) = 0.027160581698950446d0
            zpg(72) = 0.110594246078647169d0
            hpg(72) = 0.001076982794115639d0
            xpg(73) = 0.392286516401394412d0
            ypg(73) = 0.316345570093889067d0
            zpg(73) = 0.121385747834984490d0
            hpg(73) = 0.002402307691151500d0
            xpg(74) = 0.128557014198620670d0
            ypg(74) = 0.012363374710923069d0
            zpg(74) = 0.789186126910463496d0
            hpg(74) = 0.000491562317797321d0
            xpg(75) = 0.015869290737260767d0
            ypg(75) = 0.181658622840390635d0
            zpg(75) = 0.016249544899984792d0
            hpg(75) = 0.000318512330345262d0
            xpg(76) = 0.018644578551882941d0
            ypg(76) = 0.242616772316755198d0
            zpg(76) = 0.658101091100629987d0
            hpg(76) = 0.000759286556782974d0
            xpg(77) = 0.020003068069117591d0
            ypg(77) = 0.026491982786442053d0
            zpg(77) = 0.927211894993869758d0
            hpg(77) = 0.000225338516625615d0
            xpg(78) = 0.315221184417950454d0
            ypg(78) = 0.416441473581937280d0
            zpg(78) = 0.042443497933202443d0
            hpg(78) = 0.002135317888181538d0
            xpg(79) = 0.011447245328689248d0
            ypg(79) = 0.535613871061888003d0
            zpg(79) = 0.438327924764136205d0
            hpg(79) = 0.000291728153311234d0
            xpg(80) = 0.021625967230933438d0
            ypg(80) = 0.044254183440393098d0
            zpg(80) = 0.009256170098485279d0
            hpg(80) = 0.000146661872625004d0
            xpg(81) = 0.017676077824465876d0
            ypg(81) = 0.012574062625544197d0
            zpg(81) = 0.513242988713718529d0
            hpg(81) = 0.000299629038022967d0
            xpg(82) = 0.221949586639396489d0
            ypg(82) = 0.236014440187164431d0
            zpg(82) = 0.462149285715051204d0
            hpg(82) = 0.002363477360364317d0
            xpg(83) = 0.023276666837269615d0
            ypg(83) = 0.313195724922053597d0
            zpg(83) = 0.650495142146380992d0
            hpg(83) = 0.000421980838287059d0
            xpg(84) = 0.419210602486394487d0
            ypg(84) = 0.473014369981887487d0
            zpg(84) = 0.008271076205242583d0
            hpg(84) = 0.000679367221053781d0
            xpg(85) = 0.017134656780601341d0
            ypg(85) = 0.115052716338237265d0
            zpg(85) = 0.220772575137674385d0
            hpg(85) = 0.000963332581728092d0
            xpg(86) = 0.217226361056167203d0
            ypg(86) = 0.341408635799275705d0
            zpg(86) = 0.431652044801481373d0
            hpg(86) = 0.000949412961613673d0
            xpg(87) = 0.093269287523023410d0
            ypg(87) = 0.301988531624472158d0
            zpg(87) = 0.126025801385567315d0
            hpg(87) = 0.002315944894241232d0
            xpg(88) = 0.022668456071576786d0
            ypg(88) = 0.734566939691310050d0
            zpg(88) = 0.230023019380377593d0
            hpg(88) = 0.000350051089224689d0
            xpg(89) = 0.444379801443798379d0
            ypg(89) = 0.416247277852935554d0
            zpg(89) = 0.070532517126745026d0
            hpg(89) = 0.001762206803039315d0
            xpg(90) = 0.013162199306739355d0
            ypg(90) = 0.094448785984450962d0
            zpg(90) = 0.770788823962612202d0
            hpg(90) = 0.000567178942869185d0
            xpg(91) = 0.461121410824342413d0
            ypg(91) = 0.015428921605442240d0
            zpg(91) = 0.494650042211685128d0
            hpg(91) = 0.000552771153880178d0
            xpg(92) = 0.447030820647721662d0
            ypg(92) = 0.073859881488471783d0
            zpg(92) = 0.467874066123184101d0
            hpg(92) = 0.000560783361407974d0
            xpg(93) = 0.069719412228011755d0
            ypg(93) = 0.291649193825648616d0
            zpg(93) = 0.018923953609598760d0
            hpg(93) = 0.000921070641826912d0
            xpg(94) = 0.024838065539034195d0
            ypg(94) = 0.930820347861527172d0
            zpg(94) = 0.017705756642195161d0
            hpg(94) = 0.000193252576579459d0
            xpg(95) = 0.582253305811213768d0
            ypg(95) = 0.116468766015040594d0
            zpg(95) = 0.281915998269299368d0
            hpg(95) = 0.001084566925110549d0
            xpg(96) = 0.096431813358047804d0
            ypg(96) = 0.544511554928687437d0
            zpg(96) = 0.228428126919142420d0
            hpg(96) = 0.002380148484425648d0
            xpg(97) = 0.516229095681903012d0
            ypg(97) = 0.005268202700851967d0
            zpg(97) = 0.357853632242750171d0
            hpg(97) = 0.000581498057242608d0
            xpg(98) = 0.078607474439449994d0
            ypg(98) = 0.064673230527195192d0
            zpg(98) = 0.354504953752755034d0
            hpg(98) = 0.001068870515551344d0
            xpg(99) = 0.169474336773575234d0
            ypg(99) = 0.388904311457788639d0
            zpg(99) = 0.343767538744775625d0
            hpg(99) = 0.002688930384078173d0
            xpg(100) = 0.398244963074953139d0
            ypg(100) = 0.447582813298245309d0
            zpg(100) = 0.140051550890120228d0
            hpg(100) = 0.000853978511064380d0
            xpg(101) = 0.161359142629069864d0
            ypg(101) = 0.444746802500024988d0
            zpg(101) = 0.111425041693392120d0
            hpg(101) = 0.002227109480764456d0
            xpg(102) = 0.121565968423608569d0
            ypg(102) = 0.241882635470512247d0
            zpg(102) = 0.393527069207926512d0
            hpg(102) = 0.003552751837418484d0
            xpg(103) = 0.248920233145770919d0
            ypg(103) = 0.472309003805996113d0
            zpg(103) = 0.254599103746278789d0
            hpg(103) = 0.001431237107969192d0
            xpg(104) = 0.185762221331800507d0
            ypg(104) = 0.007818863945579304d0
            zpg(104) = 0.428793648278197791d0
            hpg(104) = 0.000774021806506664d0
            xpg(105) = 0.590030904858065220d0
            ypg(105) = 0.112530612275971389d0
            zpg(105) = 0.002770978521653996d0
            hpg(105) = 0.000412647551459509d0
            xpg(106) = 0.210099490511332964d0
            ypg(106) = 0.591128903555224448d0
            zpg(106) = 0.011896478778350804d0
            hpg(106) = 0.000831670580450223d0
            xpg(107) = 0.510872133826947918d0
            ypg(107) = 0.273947598644604694d0
            zpg(107) = 0.017781295177692646d0
            hpg(107) = 0.001336940295164042d0
            xpg(108) = 0.077119723738413508d0
            ypg(108) = 0.425002437923882617d0
            zpg(108) = 0.465753568992802827d0
            hpg(108) = 0.001347005005772860d0
            xpg(109) = 0.423809162068021346d0
            ypg(109) = 0.282363921682569333d0
            zpg(109) = 0.290044951449305707d0
            hpg(109) = 0.000588501422621103d0
            xpg(110) = 0.076085914100950936d0
            ypg(110) = 0.028632175468987410d0
            zpg(110) = 0.490277411613131289d0
            hpg(110) = 0.001017931028021013d0
            xpg(111) = 0.671673968685108780d0
            ypg(111) = 0.275719975402988254d0
            zpg(111) = 0.042262322063745097d0
            hpg(111) = 0.000447814461434250d0
            xpg(112) = 0.441863010418282516d0
            ypg(112) = 0.190929774807917228d0
            zpg(112) = 0.123984346835658792d0
            hpg(112) = 0.002597785869990159d0
            xpg(113) = 0.220723920933619562d0
            ypg(113) = 0.330818252129761670d0
            zpg(113) = 0.226249672118236863d0
            hpg(113) = 0.003746428030520740d0
            xpg(114) = 0.071614772278386641d0
            ypg(114) = 0.670760123898611885d0
            zpg(114) = 0.015786898512857793d0
            hpg(114) = 0.000721207502910774d0
            xpg(115) = 0.746618529093581791d0
            ypg(115) = 0.079922396051959886d0
            zpg(115) = 0.084016104679697541d0
            hpg(115) = 0.001279445938151351d0
            xpg(116) = 0.248711547764230765d0
            ypg(116) = 0.623530194691585561d0
            zpg(116) = 0.040663399131803008d0
            hpg(116) = 0.001070760507118054d0
            xpg(117) = 0.001471160462491308d0
            ypg(117) = 0.407387415634412373d0
            zpg(117) = 0.168978672956560352d0
            hpg(117) = 0.000600260386243327d0
            xpg(118) = 0.448309961464529554d0
            ypg(118) = 0.011249939859750542d0
            zpg(118) = 0.021345120905877425d0
            hpg(118) = 0.000351286760444479d0
            xpg(119) = 0.090754084754863278d0
            ypg(119) = 0.212686556478483023d0
            zpg(119) = 0.587596080355013830d0
            hpg(119) = 0.002006698611522589d0
            xpg(120) = 0.249807756482585458d0
            ypg(120) = 0.168717392641114103d0
            zpg(120) = 0.003044334322247282d0
            hpg(120) = 0.000488016463067908d0
            xpg(121) = 0.265089000310802026d0
            ypg(121) = 0.024460909200245579d0
            zpg(121) = 0.696952823961462798d0
            hpg(121) = 0.000452385668804994d0
            xpg(122) = 0.829638105523504155d0
            ypg(122) = 0.126030700634485860d0
            zpg(122) = 0.023207801609325012d0
            hpg(122) = 0.000447459356292210d0
            xpg(123) = 0.116986283554964202d0
            ypg(123) = 0.237476617912222390d0
            zpg(123) = 0.624538191516259283d0
            hpg(123) = 0.001060508107699338d0
            xpg(124) = 0.296847100633510166d0
            ypg(124) = 0.153039421627904288d0
            zpg(124) = 0.207431584504410566d0
            hpg(124) = 0.002899014729444639d0
            xpg(125) = 0.003182881770709352d0
            ypg(125) = 0.064033899632597511d0
            zpg(125) = 0.640825089075152894d0
            hpg(125) = 0.000310472009533925d0
            xpg(126) = 0.280440884540757436d0
            ypg(126) = 0.144187386881356641d0
            zpg(126) = 0.561941164655063318d0
            hpg(126) = 0.000889128430234838d0
            xpg(127) = 0.798179325749207821d0
            ypg(127) = 0.006781013385787405d0
            zpg(127) = 0.120007547789255354d0
            hpg(127) = 0.000379122126729971d0
            xpg(128) = 0.017173572059380934d0
            ypg(128) = 0.295705952299060232d0
            zpg(128) = 0.391811361425112138d0
            hpg(128) = 0.001255333413781853d0
            xpg(129) = 0.107129796593810098d0
            ypg(129) = 0.143406692064228142d0
            zpg(129) = 0.099281643449499917d0
            hpg(129) = 0.001727220605454338d0
            xpg(130) = 0.477865790351342578d0
            ypg(130) = 0.483936626223712475d0
            zpg(130) = 0.026975745802658161d0
            hpg(130) = 0.000467305129969196d0
            xpg(131) = 0.014601395984236139d0
            ypg(131) = 0.758278010901685215d0
            zpg(131) = 0.128750800246131804d0
            hpg(131) = 0.000615256171293597d0
            xpg(132) = 0.360266602050313639d0
            ypg(132) = 0.003320148291486242d0
            zpg(132) = 0.255367907216559278d0
            hpg(132) = 0.000639558998508166d0
            xpg(133) = 0.020839276773486417d0
            ypg(133) = 0.005364673957560781d0
            zpg(133) = 0.050727247638595285d0
            hpg(133) = 0.000117938436590962d0
            xpg(134) = 0.244031020980726204d0
            ypg(134) = 0.067523459458363127d0
            zpg(134) = 0.027630108161736950d0
            hpg(134) = 0.000837246113088516d0
            xpg(135) = 0.237903985814082447d0
            ypg(135) = 0.516949165930893039d0
            zpg(135) = 0.133472755072929717d0
            hpg(135) = 0.002737828174808809d0
            xpg(136) = 0.286903740327955155d0
            ypg(136) = 0.013977185724847387d0
            zpg(136) = 0.601794808958686617d0
            hpg(136) = 0.000785968411389362d0
            xpg(137) = 0.143106793911813867d0
            ypg(137) = 0.446205554078804303d0
            zpg(137) = 0.022948701329797554d0
            hpg(137) = 0.001422784558633348d0
            xpg(138) = 0.605413922421708082d0
            ypg(138) = 0.106196466093533932d0
            zpg(138) = 0.056907856627077952d0
            hpg(138) = 0.001774676463116549d0
            xpg(139) = 0.109083823415085586d0
            ypg(139) = 0.019113780778195800d0
            zpg(139) = 0.022500369138579267d0
            hpg(139) = 0.000345365336709017d0
            xpg(140) = 0.021778345009341842d0
            ypg(140) = 0.238362273040131571d0
            zpg(140) = 0.101146951654026797d0
            hpg(140) = 0.001144441720498244d0
            xpg(141) = 0.435104709856069646d0
            ypg(141) = 0.070976890181265985d0
            zpg(141) = 0.113513157977447983d0
            hpg(141) = 0.002002590695961464d0
            xpg(142) = 0.049243350221796227d0
            ypg(142) = 0.384102092848594264d0
            zpg(142) = 0.427059405703110289d0
            hpg(142) = 0.002027201241116856d0
            xpg(143) = 0.025721183835621143d0
            ypg(143) = 0.616684601606356465d0
            zpg(143) = 0.283487993375818630d0
            hpg(143) = 0.001080491883127658d0
            xpg(144) = 0.115216449883332198d0
            ypg(144) = 0.592752720304875462d0
            zpg(144) = 0.073620904664454733d0
            hpg(144) = 0.001501785201332358d0

        else if (fapg .eq. 'FPG204') then
! --------- Order 14 : High-order cubature rules for tetrahedra - 10.1002/nme.6313
            xpg(1) = 0.329815151784619308d0
            ypg(1) = 0.329815151784619308d0
            zpg(1) = 0.329815151784619308d0
            hpg(1) = 0.000574238841720564d0
            xpg(2) = 0.010554544646142073d0
            ypg(2) = 0.329815151784619308d0
            zpg(2) = 0.329815151784619308d0
            hpg(2) = 0.000574238841720564d0
            xpg(3) = 0.329815151784619308d0
            ypg(3) = 0.010554544646142073d0
            zpg(3) = 0.329815151784619308d0
            hpg(3) = 0.000574238841720564d0
            xpg(4) = 0.329815151784619308d0
            ypg(4) = 0.329815151784619308d0
            zpg(4) = 0.010554544646142073d0
            hpg(4) = 0.000574238841720564d0
            xpg(5) = 0.057538282689759202d0
            ypg(5) = 0.057538282689759202d0
            zpg(5) = 0.057538282689759202d0
            hpg(5) = 0.000436875878645082d0
            xpg(6) = 0.827385151930722392d0
            ypg(6) = 0.057538282689759202d0
            zpg(6) = 0.057538282689759202d0
            hpg(6) = 0.000436875878645082d0
            xpg(7) = 0.057538282689759202d0
            ypg(7) = 0.827385151930722392d0
            zpg(7) = 0.057538282689759202d0
            hpg(7) = 0.000436875878645082d0
            xpg(8) = 0.057538282689759202d0
            ypg(8) = 0.057538282689759202d0
            zpg(8) = 0.827385151930722392d0
            hpg(8) = 0.000436875878645082d0
            xpg(9) = 0.005856168613783504d0
            ypg(9) = 0.005856168613783504d0
            zpg(9) = 0.005856168613783504d0
            hpg(9) = 1.845638736489171d-05
            xpg(10) = 0.982431494158649487d0
            ypg(10) = 0.005856168613783504d0
            zpg(10) = 0.005856168613783504d0
            hpg(10) = 1.845638736489171d-05
            xpg(11) = 0.005856168613783504d0
            ypg(11) = 0.982431494158649487d0
            zpg(11) = 0.005856168613783504d0
            hpg(11) = 1.845638736489171d-05
            xpg(12) = 0.005856168613783504d0
            ypg(12) = 0.005856168613783504d0
            zpg(12) = 0.982431494158649487d0
            hpg(12) = 1.845638736489171d-05
            xpg(13) = 0.160555475847956640d0
            ypg(13) = 0.160555475847956640d0
            zpg(13) = 0.160555475847956640d0
            hpg(13) = 0.001671430895412300d0
            xpg(14) = 0.518333572456130079d0
            ypg(14) = 0.160555475847956640d0
            zpg(14) = 0.160555475847956640d0
            hpg(14) = 0.001671430895412300d0
            xpg(15) = 0.160555475847956640d0
            ypg(15) = 0.518333572456130079d0
            zpg(15) = 0.160555475847956640d0
            hpg(15) = 0.001671430895412300d0
            xpg(16) = 0.160555475847956640d0
            ypg(16) = 0.160555475847956640d0
            zpg(16) = 0.518333572456130079d0
            hpg(16) = 0.001671430895412300d0
            xpg(17) = 0.098739646074049086d0
            ypg(17) = 0.098739646074049086d0
            zpg(17) = 0.098739646074049086d0
            hpg(17) = 0.001131734742968235d0
            xpg(18) = 0.703781061777852739d0
            ypg(18) = 0.098739646074049086d0
            zpg(18) = 0.098739646074049086d0
            hpg(18) = 0.001131734742968235d0
            xpg(19) = 0.098739646074049086d0
            ypg(19) = 0.703781061777852739d0
            zpg(19) = 0.098739646074049086d0
            hpg(19) = 0.001131734742968235d0
            xpg(20) = 0.098739646074049086d0
            ypg(20) = 0.098739646074049086d0
            zpg(20) = 0.703781061777852739d0
            hpg(20) = 0.001131734742968235d0
            xpg(21) = 0.208053196159726508d0
            ypg(21) = 0.208053196159726508d0
            zpg(21) = 0.208053196159726508d0
            hpg(21) = 0.001186358056983214d0
            xpg(22) = 0.375840411520820475d0
            ypg(22) = 0.208053196159726508d0
            zpg(22) = 0.208053196159726508d0
            hpg(22) = 0.001186358056983214d0
            xpg(23) = 0.208053196159726508d0
            ypg(23) = 0.375840411520820475d0
            zpg(23) = 0.208053196159726508d0
            hpg(23) = 0.001186358056983214d0
            xpg(24) = 0.208053196159726508d0
            ypg(24) = 0.208053196159726508d0
            zpg(24) = 0.375840411520820475d0
            hpg(24) = 0.001186358056983214d0
            xpg(25) = 0.009752028812223523d0
            ypg(25) = 0.490247971187776476d0
            zpg(25) = 0.490247971187776476d0
            hpg(25) = 0.000176675382427341d0
            xpg(26) = 0.490247971187776476d0
            ypg(26) = 0.009752028812223523d0
            zpg(26) = 0.490247971187776476d0
            hpg(26) = 0.000176675382427341d0
            xpg(27) = 0.009752028812223523d0
            ypg(27) = 0.009752028812223523d0
            zpg(27) = 0.490247971187776476d0
            hpg(27) = 0.000176675382427341d0
            xpg(28) = 0.490247971187776476d0
            ypg(28) = 0.490247971187776476d0
            zpg(28) = 0.009752028812223523d0
            hpg(28) = 0.000176675382427341d0
            xpg(29) = 0.009752028812223523d0
            ypg(29) = 0.490247971187776476d0
            zpg(29) = 0.009752028812223523d0
            hpg(29) = 0.000176675382427341d0
            xpg(30) = 0.490247971187776476d0
            ypg(30) = 0.009752028812223523d0
            zpg(30) = 0.009752028812223523d0
            hpg(30) = 0.000176675382427341d0
            xpg(31) = 0.397060599884467409d0
            ypg(31) = 0.102939400115532590d0
            zpg(31) = 0.102939400115532590d0
            hpg(31) = 0.002278305165830148d0
            xpg(32) = 0.102939400115532590d0
            ypg(32) = 0.397060599884467409d0
            zpg(32) = 0.102939400115532590d0
            hpg(32) = 0.002278305165830148d0
            xpg(33) = 0.397060599884467409d0
            ypg(33) = 0.397060599884467409d0
            zpg(33) = 0.102939400115532590d0
            hpg(33) = 0.002278305165830148d0
            xpg(34) = 0.102939400115532590d0
            ypg(34) = 0.102939400115532590d0
            zpg(34) = 0.397060599884467409d0
            hpg(34) = 0.002278305165830148d0
            xpg(35) = 0.397060599884467409d0
            ypg(35) = 0.102939400115532590d0
            zpg(35) = 0.397060599884467409d0
            hpg(35) = 0.002278305165830148d0
            xpg(36) = 0.102939400115532590d0
            ypg(36) = 0.397060599884467409d0
            zpg(36) = 0.397060599884467409d0
            hpg(36) = 0.002278305165830148d0
            xpg(37) = 0.462086536538784864d0
            ypg(37) = 0.037913463461215135d0
            zpg(37) = 0.037913463461215135d0
            hpg(37) = 0.000448456849112692d0
            xpg(38) = 0.037913463461215135d0
            ypg(38) = 0.462086536538784864d0
            zpg(38) = 0.037913463461215135d0
            hpg(38) = 0.000448456849112692d0
            xpg(39) = 0.462086536538784864d0
            ypg(39) = 0.462086536538784864d0
            zpg(39) = 0.037913463461215135d0
            hpg(39) = 0.000448456849112692d0
            xpg(40) = 0.037913463461215135d0
            ypg(40) = 0.037913463461215135d0
            zpg(40) = 0.462086536538784864d0
            hpg(40) = 0.000448456849112692d0
            xpg(41) = 0.462086536538784864d0
            ypg(41) = 0.037913463461215135d0
            zpg(41) = 0.462086536538784864d0
            hpg(41) = 0.000448456849112692d0
            xpg(42) = 0.037913463461215135d0
            ypg(42) = 0.462086536538784864d0
            zpg(42) = 0.462086536538784864d0
            hpg(42) = 0.000448456849112692d0
            xpg(43) = 0.181426423816139620d0
            ypg(43) = 0.318573576183860379d0
            zpg(43) = 0.318573576183860379d0
            hpg(43) = 0.001694100204307123d0
            xpg(44) = 0.318573576183860379d0
            ypg(44) = 0.181426423816139620d0
            zpg(44) = 0.318573576183860379d0
            hpg(44) = 0.001694100204307123d0
            xpg(45) = 0.181426423816139620d0
            ypg(45) = 0.181426423816139620d0
            zpg(45) = 0.318573576183860379d0
            hpg(45) = 0.001694100204307123d0
            xpg(46) = 0.318573576183860379d0
            ypg(46) = 0.318573576183860379d0
            zpg(46) = 0.181426423816139620d0
            hpg(46) = 0.001694100204307123d0
            xpg(47) = 0.181426423816139620d0
            ypg(47) = 0.318573576183860379d0
            zpg(47) = 0.181426423816139620d0
            hpg(47) = 0.001694100204307123d0
            xpg(48) = 0.318573576183860379d0
            ypg(48) = 0.181426423816139620d0
            zpg(48) = 0.181426423816139620d0
            hpg(48) = 0.001694100204307123d0
            xpg(49) = 0.418938778837666962d0
            ypg(49) = 0.249909218863499831d0
            zpg(49) = 0.249909218863499831d0
            hpg(49) = 0.002257399939621233d0
            xpg(50) = 0.249909218863499831d0
            ypg(50) = 0.418938778837666962d0
            zpg(50) = 0.249909218863499831d0
            hpg(50) = 0.002257399939621233d0
            xpg(51) = 0.249909218863499831d0
            ypg(51) = 0.249909218863499831d0
            zpg(51) = 0.418938778837666962d0
            hpg(51) = 0.002257399939621233d0
            xpg(52) = 0.081242783435333373d0
            ypg(52) = 0.249909218863499831d0
            zpg(52) = 0.249909218863499831d0
            hpg(52) = 0.002257399939621233d0
            xpg(53) = 0.081242783435333373d0
            ypg(53) = 0.418938778837666962d0
            zpg(53) = 0.249909218863499831d0
            hpg(53) = 0.002257399939621233d0
            xpg(54) = 0.081242783435333373d0
            ypg(54) = 0.249909218863499831d0
            zpg(54) = 0.418938778837666962d0
            hpg(54) = 0.002257399939621233d0
            xpg(55) = 0.249909218863499831d0
            ypg(55) = 0.081242783435333373d0
            zpg(55) = 0.249909218863499831d0
            hpg(55) = 0.002257399939621233d0
            xpg(56) = 0.418938778837666962d0
            ypg(56) = 0.081242783435333373d0
            zpg(56) = 0.249909218863499831d0
            hpg(56) = 0.002257399939621233d0
            xpg(57) = 0.249909218863499831d0
            ypg(57) = 0.081242783435333373d0
            zpg(57) = 0.418938778837666962d0
            hpg(57) = 0.002257399939621233d0
            xpg(58) = 0.249909218863499831d0
            ypg(58) = 0.249909218863499831d0
            zpg(58) = 0.081242783435333373d0
            hpg(58) = 0.002257399939621233d0
            xpg(59) = 0.418938778837666962d0
            ypg(59) = 0.249909218863499831d0
            zpg(59) = 0.081242783435333373d0
            hpg(59) = 0.002257399939621233d0
            xpg(60) = 0.249909218863499831d0
            ypg(60) = 0.418938778837666962d0
            zpg(60) = 0.081242783435333373d0
            hpg(60) = 0.002257399939621233d0
            xpg(61) = 0.562876680249783716d0
            ypg(61) = 0.213263378061875731d0
            zpg(61) = 0.213263378061875731d0
            hpg(61) = 0.000770556674611631d0
            xpg(62) = 0.213263378061875731d0
            ypg(62) = 0.562876680249783716d0
            zpg(62) = 0.213263378061875731d0
            hpg(62) = 0.000770556674611631d0
            xpg(63) = 0.213263378061875731d0
            ypg(63) = 0.213263378061875731d0
            zpg(63) = 0.562876680249783716d0
            hpg(63) = 0.000770556674611631d0
            xpg(64) = 0.010596563626464820d0
            ypg(64) = 0.213263378061875731d0
            zpg(64) = 0.213263378061875731d0
            hpg(64) = 0.000770556674611631d0
            xpg(65) = 0.010596563626464820d0
            ypg(65) = 0.562876680249783716d0
            zpg(65) = 0.213263378061875731d0
            hpg(65) = 0.000770556674611631d0
            xpg(66) = 0.010596563626464820d0
            ypg(66) = 0.213263378061875731d0
            zpg(66) = 0.562876680249783716d0
            hpg(66) = 0.000770556674611631d0
            xpg(67) = 0.213263378061875731d0
            ypg(67) = 0.010596563626464820d0
            zpg(67) = 0.213263378061875731d0
            hpg(67) = 0.000770556674611631d0
            xpg(68) = 0.562876680249783716d0
            ypg(68) = 0.010596563626464820d0
            zpg(68) = 0.213263378061875731d0
            hpg(68) = 0.000770556674611631d0
            xpg(69) = 0.213263378061875731d0
            ypg(69) = 0.010596563626464820d0
            zpg(69) = 0.562876680249783716d0
            hpg(69) = 0.000770556674611631d0
            xpg(70) = 0.213263378061875731d0
            ypg(70) = 0.213263378061875731d0
            zpg(70) = 0.010596563626464820d0
            hpg(70) = 0.000770556674611631d0
            xpg(71) = 0.562876680249783716d0
            ypg(71) = 0.213263378061875731d0
            zpg(71) = 0.010596563626464820d0
            hpg(71) = 0.000770556674611631d0
            xpg(72) = 0.213263378061875731d0
            ypg(72) = 0.562876680249783716d0
            zpg(72) = 0.010596563626464820d0
            hpg(72) = 0.000770556674611631d0
            xpg(73) = 0.895468972714016204d0
            ypg(73) = 0.049181284940159052d0
            zpg(73) = 0.049181284940159052d0
            hpg(73) = 0.000138858307541221d0
            xpg(74) = 0.049181284940159052d0
            ypg(74) = 0.895468972714016204d0
            zpg(74) = 0.049181284940159052d0
            hpg(74) = 0.000138858307541221d0
            xpg(75) = 0.049181284940159052d0
            ypg(75) = 0.049181284940159052d0
            zpg(75) = 0.895468972714016204d0
            hpg(75) = 0.000138858307541221d0
            xpg(76) = 0.006168457405665689d0
            ypg(76) = 0.049181284940159052d0
            zpg(76) = 0.049181284940159052d0
            hpg(76) = 0.000138858307541221d0
            xpg(77) = 0.006168457405665689d0
            ypg(77) = 0.895468972714016204d0
            zpg(77) = 0.049181284940159052d0
            hpg(77) = 0.000138858307541221d0
            xpg(78) = 0.006168457405665689d0
            ypg(78) = 0.049181284940159052d0
            zpg(78) = 0.895468972714016204d0
            hpg(78) = 0.000138858307541221d0
            xpg(79) = 0.049181284940159052d0
            ypg(79) = 0.006168457405665689d0
            zpg(79) = 0.049181284940159052d0
            hpg(79) = 0.000138858307541221d0
            xpg(80) = 0.895468972714016204d0
            ypg(80) = 0.006168457405665689d0
            zpg(80) = 0.049181284940159052d0
            hpg(80) = 0.000138858307541221d0
            xpg(81) = 0.049181284940159052d0
            ypg(81) = 0.006168457405665689d0
            zpg(81) = 0.895468972714016204d0
            hpg(81) = 0.000138858307541221d0
            xpg(82) = 0.049181284940159052d0
            ypg(82) = 0.049181284940159052d0
            zpg(82) = 0.006168457405665689d0
            hpg(82) = 0.000138858307541221d0
            xpg(83) = 0.895468972714016204d0
            ypg(83) = 0.049181284940159052d0
            zpg(83) = 0.006168457405665689d0
            hpg(83) = 0.000138858307541221d0
            xpg(84) = 0.049181284940159052d0
            ypg(84) = 0.895468972714016204d0
            zpg(84) = 0.006168457405665689d0
            hpg(84) = 0.000138858307541221d0
            xpg(85) = 0.020663578929673465d0
            ypg(85) = 0.392862617970060138d0
            zpg(85) = 0.392862617970060138d0
            hpg(85) = 0.001206716126809541d0
            xpg(86) = 0.392862617970060138d0
            ypg(86) = 0.020663578929673465d0
            zpg(86) = 0.392862617970060138d0
            hpg(86) = 0.001206716126809541d0
            xpg(87) = 0.392862617970060138d0
            ypg(87) = 0.392862617970060138d0
            zpg(87) = 0.020663578929673465d0
            hpg(87) = 0.001206716126809541d0
            xpg(88) = 0.193611185130206257d0
            ypg(88) = 0.392862617970060138d0
            zpg(88) = 0.392862617970060138d0
            hpg(88) = 0.001206716126809541d0
            xpg(89) = 0.193611185130206257d0
            ypg(89) = 0.020663578929673465d0
            zpg(89) = 0.392862617970060138d0
            hpg(89) = 0.001206716126809541d0
            xpg(90) = 0.193611185130206257d0
            ypg(90) = 0.392862617970060138d0
            zpg(90) = 0.020663578929673465d0
            hpg(90) = 0.001206716126809541d0
            xpg(91) = 0.392862617970060138d0
            ypg(91) = 0.193611185130206257d0
            zpg(91) = 0.392862617970060138d0
            hpg(91) = 0.001206716126809541d0
            xpg(92) = 0.020663578929673465d0
            ypg(92) = 0.193611185130206257d0
            zpg(92) = 0.392862617970060138d0
            hpg(92) = 0.001206716126809541d0
            xpg(93) = 0.392862617970060138d0
            ypg(93) = 0.193611185130206257d0
            zpg(93) = 0.020663578929673465d0
            hpg(93) = 0.001206716126809541d0
            xpg(94) = 0.392862617970060138d0
            ypg(94) = 0.392862617970060138d0
            zpg(94) = 0.193611185130206257d0
            hpg(94) = 0.001206716126809541d0
            xpg(95) = 0.020663578929673465d0
            ypg(95) = 0.392862617970060138d0
            zpg(95) = 0.193611185130206257d0
            hpg(95) = 0.001206716126809541d0
            xpg(96) = 0.392862617970060138d0
            ypg(96) = 0.020663578929673465d0
            zpg(96) = 0.193611185130206257d0
            hpg(96) = 0.001206716126809541d0
            xpg(97) = 0.126902202407495768d0
            ypg(97) = 0.012617605532570709d0
            zpg(97) = 0.012617605532570709d0
            hpg(97) = 0.000145863463499634d0
            xpg(98) = 0.012617605532570709d0
            ypg(98) = 0.126902202407495768d0
            zpg(98) = 0.012617605532570709d0
            hpg(98) = 0.000145863463499634d0
            xpg(99) = 0.012617605532570709d0
            ypg(99) = 0.012617605532570709d0
            zpg(99) = 0.126902202407495768d0
            hpg(99) = 0.000145863463499634d0
            xpg(100) = 0.847862586527362812d0
            ypg(100) = 0.012617605532570709d0
            zpg(100) = 0.012617605532570709d0
            hpg(100) = 0.000145863463499634d0
            xpg(101) = 0.847862586527362812d0
            ypg(101) = 0.126902202407495768d0
            zpg(101) = 0.012617605532570709d0
            hpg(101) = 0.000145863463499634d0
            xpg(102) = 0.847862586527362812d0
            ypg(102) = 0.012617605532570709d0
            zpg(102) = 0.126902202407495768d0
            hpg(102) = 0.000145863463499634d0
            xpg(103) = 0.012617605532570709d0
            ypg(103) = 0.847862586527362812d0
            zpg(103) = 0.012617605532570709d0
            hpg(103) = 0.000145863463499634d0
            xpg(104) = 0.126902202407495768d0
            ypg(104) = 0.847862586527362812d0
            zpg(104) = 0.012617605532570709d0
            hpg(104) = 0.000145863463499634d0
            xpg(105) = 0.012617605532570709d0
            ypg(105) = 0.847862586527362812d0
            zpg(105) = 0.126902202407495768d0
            hpg(105) = 0.000145863463499634d0
            xpg(106) = 0.012617605532570709d0
            ypg(106) = 0.012617605532570709d0
            zpg(106) = 0.847862586527362812d0
            hpg(106) = 0.000145863463499634d0
            xpg(107) = 0.126902202407495768d0
            ypg(107) = 0.012617605532570709d0
            zpg(107) = 0.847862586527362812d0
            hpg(107) = 0.000145863463499634d0
            xpg(108) = 0.012617605532570709d0
            ypg(108) = 0.126902202407495768d0
            zpg(108) = 0.847862586527362812d0
            hpg(108) = 0.000145863463499634d0
            xpg(109) = 0.062140003117621531d0
            ypg(109) = 0.121618784348698988d0
            zpg(109) = 0.576740617455094747d0
            hpg(109) = 0.001195422753016406d0
            xpg(110) = 0.062140003117621531d0
            ypg(110) = 0.576740617455094747d0
            zpg(110) = 0.121618784348698988d0
            hpg(110) = 0.001195422753016406d0
            xpg(111) = 0.121618784348698988d0
            ypg(111) = 0.062140003117621531d0
            zpg(111) = 0.576740617455094747d0
            hpg(111) = 0.001195422753016406d0
            xpg(112) = 0.576740617455094747d0
            ypg(112) = 0.062140003117621531d0
            zpg(112) = 0.121618784348698988d0
            hpg(112) = 0.001195422753016406d0
            xpg(113) = 0.121618784348698988d0
            ypg(113) = 0.576740617455094747d0
            zpg(113) = 0.062140003117621531d0
            hpg(113) = 0.001195422753016406d0
            xpg(114) = 0.576740617455094747d0
            ypg(114) = 0.121618784348698988d0
            zpg(114) = 0.062140003117621531d0
            hpg(114) = 0.001195422753016406d0
            xpg(115) = 0.239500595078584732d0
            ypg(115) = 0.121618784348698988d0
            zpg(115) = 0.576740617455094747d0
            hpg(115) = 0.001195422753016406d0
            xpg(116) = 0.239500595078584732d0
            ypg(116) = 0.576740617455094747d0
            zpg(116) = 0.121618784348698988d0
            hpg(116) = 0.001195422753016406d0
            xpg(117) = 0.239500595078584732d0
            ypg(117) = 0.062140003117621531d0
            zpg(117) = 0.576740617455094747d0
            hpg(117) = 0.001195422753016406d0
            xpg(118) = 0.239500595078584732d0
            ypg(118) = 0.062140003117621531d0
            zpg(118) = 0.121618784348698988d0
            hpg(118) = 0.001195422753016406d0
            xpg(119) = 0.239500595078584732d0
            ypg(119) = 0.576740617455094747d0
            zpg(119) = 0.062140003117621531d0
            hpg(119) = 0.001195422753016406d0
            xpg(120) = 0.239500595078584732d0
            ypg(120) = 0.121618784348698988d0
            zpg(120) = 0.062140003117621531d0
            hpg(120) = 0.001195422753016406d0
            xpg(121) = 0.121618784348698988d0
            ypg(121) = 0.239500595078584732d0
            zpg(121) = 0.576740617455094747d0
            hpg(121) = 0.001195422753016406d0
            xpg(122) = 0.576740617455094747d0
            ypg(122) = 0.239500595078584732d0
            zpg(122) = 0.121618784348698988d0
            hpg(122) = 0.001195422753016406d0
            xpg(123) = 0.062140003117621531d0
            ypg(123) = 0.239500595078584732d0
            zpg(123) = 0.576740617455094747d0
            hpg(123) = 0.001195422753016406d0
            xpg(124) = 0.062140003117621531d0
            ypg(124) = 0.239500595078584732d0
            zpg(124) = 0.121618784348698988d0
            hpg(124) = 0.001195422753016406d0
            xpg(125) = 0.576740617455094747d0
            ypg(125) = 0.239500595078584732d0
            zpg(125) = 0.062140003117621531d0
            hpg(125) = 0.001195422753016406d0
            xpg(126) = 0.121618784348698988d0
            ypg(126) = 0.239500595078584732d0
            zpg(126) = 0.062140003117621531d0
            hpg(126) = 0.001195422753016406d0
            xpg(127) = 0.121618784348698988d0
            ypg(127) = 0.576740617455094747d0
            zpg(127) = 0.239500595078584732d0
            hpg(127) = 0.001195422753016406d0
            xpg(128) = 0.576740617455094747d0
            ypg(128) = 0.121618784348698988d0
            zpg(128) = 0.239500595078584732d0
            hpg(128) = 0.001195422753016406d0
            xpg(129) = 0.062140003117621531d0
            ypg(129) = 0.576740617455094747d0
            zpg(129) = 0.239500595078584732d0
            hpg(129) = 0.001195422753016406d0
            xpg(130) = 0.062140003117621531d0
            ypg(130) = 0.121618784348698988d0
            zpg(130) = 0.239500595078584732d0
            hpg(130) = 0.001195422753016406d0
            xpg(131) = 0.576740617455094747d0
            ypg(131) = 0.062140003117621531d0
            zpg(131) = 0.239500595078584732d0
            hpg(131) = 0.001195422753016406d0
            xpg(132) = 0.121618784348698988d0
            ypg(132) = 0.062140003117621531d0
            zpg(132) = 0.239500595078584732d0
            hpg(132) = 0.001195422753016406d0
            xpg(133) = 0.551794132118424634d0
            ypg(133) = 0.348902296047024577d0
            zpg(133) = 0.015522119223207871d0
            hpg(133) = 0.000714476057131796d0
            xpg(134) = 0.551794132118424634d0
            ypg(134) = 0.015522119223207871d0
            zpg(134) = 0.348902296047024577d0
            hpg(134) = 0.000714476057131796d0
            xpg(135) = 0.348902296047024577d0
            ypg(135) = 0.551794132118424634d0
            zpg(135) = 0.015522119223207871d0
            hpg(135) = 0.000714476057131796d0
            xpg(136) = 0.015522119223207871d0
            ypg(136) = 0.551794132118424634d0
            zpg(136) = 0.348902296047024577d0
            hpg(136) = 0.000714476057131796d0
            xpg(137) = 0.348902296047024577d0
            ypg(137) = 0.015522119223207871d0
            zpg(137) = 0.551794132118424634d0
            hpg(137) = 0.000714476057131796d0
            xpg(138) = 0.015522119223207871d0
            ypg(138) = 0.348902296047024577d0
            zpg(138) = 0.551794132118424634d0
            hpg(138) = 0.000714476057131796d0
            xpg(139) = 0.083781452611342916d0
            ypg(139) = 0.348902296047024577d0
            zpg(139) = 0.015522119223207871d0
            hpg(139) = 0.000714476057131796d0
            xpg(140) = 0.083781452611342916d0
            ypg(140) = 0.015522119223207871d0
            zpg(140) = 0.348902296047024577d0
            hpg(140) = 0.000714476057131796d0
            xpg(141) = 0.083781452611342916d0
            ypg(141) = 0.551794132118424634d0
            zpg(141) = 0.015522119223207871d0
            hpg(141) = 0.000714476057131796d0
            xpg(142) = 0.083781452611342916d0
            ypg(142) = 0.551794132118424634d0
            zpg(142) = 0.348902296047024577d0
            hpg(142) = 0.000714476057131796d0
            xpg(143) = 0.083781452611342916d0
            ypg(143) = 0.015522119223207871d0
            zpg(143) = 0.551794132118424634d0
            hpg(143) = 0.000714476057131796d0
            xpg(144) = 0.083781452611342916d0
            ypg(144) = 0.348902296047024577d0
            zpg(144) = 0.551794132118424634d0
            hpg(144) = 0.000714476057131796d0
            xpg(145) = 0.348902296047024577d0
            ypg(145) = 0.083781452611342916d0
            zpg(145) = 0.015522119223207871d0
            hpg(145) = 0.000714476057131796d0
            xpg(146) = 0.015522119223207871d0
            ypg(146) = 0.083781452611342916d0
            zpg(146) = 0.348902296047024577d0
            hpg(146) = 0.000714476057131796d0
            xpg(147) = 0.551794132118424634d0
            ypg(147) = 0.083781452611342916d0
            zpg(147) = 0.015522119223207871d0
            hpg(147) = 0.000714476057131796d0
            xpg(148) = 0.551794132118424634d0
            ypg(148) = 0.083781452611342916d0
            zpg(148) = 0.348902296047024577d0
            hpg(148) = 0.000714476057131796d0
            xpg(149) = 0.015522119223207871d0
            ypg(149) = 0.083781452611342916d0
            zpg(149) = 0.551794132118424634d0
            hpg(149) = 0.000714476057131796d0
            xpg(150) = 0.348902296047024577d0
            ypg(150) = 0.083781452611342916d0
            zpg(150) = 0.551794132118424634d0
            hpg(150) = 0.000714476057131796d0
            xpg(151) = 0.348902296047024577d0
            ypg(151) = 0.015522119223207871d0
            zpg(151) = 0.083781452611342916d0
            hpg(151) = 0.000714476057131796d0
            xpg(152) = 0.015522119223207871d0
            ypg(152) = 0.348902296047024577d0
            zpg(152) = 0.083781452611342916d0
            hpg(152) = 0.000714476057131796d0
            xpg(153) = 0.551794132118424634d0
            ypg(153) = 0.015522119223207871d0
            zpg(153) = 0.083781452611342916d0
            hpg(153) = 0.000714476057131796d0
            xpg(154) = 0.551794132118424634d0
            ypg(154) = 0.348902296047024577d0
            zpg(154) = 0.083781452611342916d0
            hpg(154) = 0.000714476057131796d0
            xpg(155) = 0.015522119223207871d0
            ypg(155) = 0.551794132118424634d0
            zpg(155) = 0.083781452611342916d0
            hpg(155) = 0.000714476057131796d0
            xpg(156) = 0.348902296047024577d0
            ypg(156) = 0.551794132118424634d0
            zpg(156) = 0.083781452611342916d0
            hpg(156) = 0.000714476057131796d0
            xpg(157) = 0.017841640647654563d0
            ypg(157) = 0.682238320168289338d0
            zpg(157) = 0.015479297313172471d0
            hpg(157) = 0.000162625460691781d0
            xpg(158) = 0.017841640647654563d0
            ypg(158) = 0.015479297313172471d0
            zpg(158) = 0.682238320168289338d0
            hpg(158) = 0.000162625460691781d0
            xpg(159) = 0.682238320168289338d0
            ypg(159) = 0.017841640647654563d0
            zpg(159) = 0.015479297313172471d0
            hpg(159) = 0.000162625460691781d0
            xpg(160) = 0.015479297313172471d0
            ypg(160) = 0.017841640647654563d0
            zpg(160) = 0.682238320168289338d0
            hpg(160) = 0.000162625460691781d0
            xpg(161) = 0.682238320168289338d0
            ypg(161) = 0.015479297313172471d0
            zpg(161) = 0.017841640647654563d0
            hpg(161) = 0.000162625460691781d0
            xpg(162) = 0.015479297313172471d0
            ypg(162) = 0.682238320168289338d0
            zpg(162) = 0.017841640647654563d0
            hpg(162) = 0.000162625460691781d0
            xpg(163) = 0.284440741870883626d0
            ypg(163) = 0.682238320168289338d0
            zpg(163) = 0.015479297313172471d0
            hpg(163) = 0.000162625460691781d0
            xpg(164) = 0.284440741870883626d0
            ypg(164) = 0.015479297313172471d0
            zpg(164) = 0.682238320168289338d0
            hpg(164) = 0.000162625460691781d0
            xpg(165) = 0.284440741870883626d0
            ypg(165) = 0.017841640647654563d0
            zpg(165) = 0.015479297313172471d0
            hpg(165) = 0.000162625460691781d0
            xpg(166) = 0.284440741870883626d0
            ypg(166) = 0.017841640647654563d0
            zpg(166) = 0.682238320168289338d0
            hpg(166) = 0.000162625460691781d0
            xpg(167) = 0.284440741870883626d0
            ypg(167) = 0.015479297313172471d0
            zpg(167) = 0.017841640647654563d0
            hpg(167) = 0.000162625460691781d0
            xpg(168) = 0.284440741870883626d0
            ypg(168) = 0.682238320168289338d0
            zpg(168) = 0.017841640647654563d0
            hpg(168) = 0.000162625460691781d0
            xpg(169) = 0.682238320168289338d0
            ypg(169) = 0.284440741870883626d0
            zpg(169) = 0.015479297313172471d0
            hpg(169) = 0.000162625460691781d0
            xpg(170) = 0.015479297313172471d0
            ypg(170) = 0.284440741870883626d0
            zpg(170) = 0.682238320168289338d0
            hpg(170) = 0.000162625460691781d0
            xpg(171) = 0.017841640647654563d0
            ypg(171) = 0.284440741870883626d0
            zpg(171) = 0.015479297313172471d0
            hpg(171) = 0.000162625460691781d0
            xpg(172) = 0.017841640647654563d0
            ypg(172) = 0.284440741870883626d0
            zpg(172) = 0.682238320168289338d0
            hpg(172) = 0.000162625460691781d0
            xpg(173) = 0.015479297313172471d0
            ypg(173) = 0.284440741870883626d0
            zpg(173) = 0.017841640647654563d0
            hpg(173) = 0.000162625460691781d0
            xpg(174) = 0.682238320168289338d0
            ypg(174) = 0.284440741870883626d0
            zpg(174) = 0.017841640647654563d0
            hpg(174) = 0.000162625460691781d0
            xpg(175) = 0.682238320168289338d0
            ypg(175) = 0.015479297313172471d0
            zpg(175) = 0.284440741870883626d0
            hpg(175) = 0.000162625460691781d0
            xpg(176) = 0.015479297313172471d0
            ypg(176) = 0.682238320168289338d0
            zpg(176) = 0.284440741870883626d0
            hpg(176) = 0.000162625460691781d0
            xpg(177) = 0.017841640647654563d0
            ypg(177) = 0.015479297313172471d0
            zpg(177) = 0.284440741870883626d0
            hpg(177) = 0.000162625460691781d0
            xpg(178) = 0.017841640647654563d0
            ypg(178) = 0.682238320168289338d0
            zpg(178) = 0.284440741870883626d0
            hpg(178) = 0.000162625460691781d0
            xpg(179) = 0.015479297313172471d0
            ypg(179) = 0.017841640647654563d0
            zpg(179) = 0.284440741870883626d0
            hpg(179) = 0.000162625460691781d0
            xpg(180) = 0.682238320168289338d0
            ypg(180) = 0.017841640647654563d0
            zpg(180) = 0.284440741870883626d0
            hpg(180) = 0.000162625460691781d0
            xpg(181) = 0.169539550700162253d0
            ypg(181) = 0.733998298658581071d0
            zpg(181) = 0.017776602770142781d0
            hpg(181) = 0.000626322716627787d0
            xpg(182) = 0.169539550700162253d0
            ypg(182) = 0.017776602770142781d0
            zpg(182) = 0.733998298658581071d0
            hpg(182) = 0.000626322716627787d0
            xpg(183) = 0.733998298658581071d0
            ypg(183) = 0.169539550700162253d0
            zpg(183) = 0.017776602770142781d0
            hpg(183) = 0.000626322716627787d0
            xpg(184) = 0.017776602770142781d0
            ypg(184) = 0.169539550700162253d0
            zpg(184) = 0.733998298658581071d0
            hpg(184) = 0.000626322716627787d0
            xpg(185) = 0.733998298658581071d0
            ypg(185) = 0.017776602770142781d0
            zpg(185) = 0.169539550700162253d0
            hpg(185) = 0.000626322716627787d0
            xpg(186) = 0.017776602770142781d0
            ypg(186) = 0.733998298658581071d0
            zpg(186) = 0.169539550700162253d0
            hpg(186) = 0.000626322716627787d0
            xpg(187) = 0.078685547871113893d0
            ypg(187) = 0.733998298658581071d0
            zpg(187) = 0.017776602770142781d0
            hpg(187) = 0.000626322716627787d0
            xpg(188) = 0.078685547871113893d0
            ypg(188) = 0.017776602770142781d0
            zpg(188) = 0.733998298658581071d0
            hpg(188) = 0.000626322716627787d0
            xpg(189) = 0.078685547871113893d0
            ypg(189) = 0.169539550700162253d0
            zpg(189) = 0.017776602770142781d0
            hpg(189) = 0.000626322716627787d0
            xpg(190) = 0.078685547871113893d0
            ypg(190) = 0.169539550700162253d0
            zpg(190) = 0.733998298658581071d0
            hpg(190) = 0.000626322716627787d0
            xpg(191) = 0.078685547871113893d0
            ypg(191) = 0.017776602770142781d0
            zpg(191) = 0.169539550700162253d0
            hpg(191) = 0.000626322716627787d0
            xpg(192) = 0.078685547871113893d0
            ypg(192) = 0.733998298658581071d0
            zpg(192) = 0.169539550700162253d0
            hpg(192) = 0.000626322716627787d0
            xpg(193) = 0.733998298658581071d0
            ypg(193) = 0.078685547871113893d0
            zpg(193) = 0.017776602770142781d0
            hpg(193) = 0.000626322716627787d0
            xpg(194) = 0.017776602770142781d0
            ypg(194) = 0.078685547871113893d0
            zpg(194) = 0.733998298658581071d0
            hpg(194) = 0.000626322716627787d0
            xpg(195) = 0.169539550700162253d0
            ypg(195) = 0.078685547871113893d0
            zpg(195) = 0.017776602770142781d0
            hpg(195) = 0.000626322716627787d0
            xpg(196) = 0.169539550700162253d0
            ypg(196) = 0.078685547871113893d0
            zpg(196) = 0.733998298658581071d0
            hpg(196) = 0.000626322716627787d0
            xpg(197) = 0.017776602770142781d0
            ypg(197) = 0.078685547871113893d0
            zpg(197) = 0.169539550700162253d0
            hpg(197) = 0.000626322716627787d0
            xpg(198) = 0.733998298658581071d0
            ypg(198) = 0.078685547871113893d0
            zpg(198) = 0.169539550700162253d0
            hpg(198) = 0.000626322716627787d0
            xpg(199) = 0.733998298658581071d0
            ypg(199) = 0.017776602770142781d0
            zpg(199) = 0.078685547871113893d0
            hpg(199) = 0.000626322716627787d0
            xpg(200) = 0.017776602770142781d0
            ypg(200) = 0.733998298658581071d0
            zpg(200) = 0.078685547871113893d0
            hpg(200) = 0.000626322716627787d0
            xpg(201) = 0.169539550700162253d0
            ypg(201) = 0.017776602770142781d0
            zpg(201) = 0.078685547871113893d0
            hpg(201) = 0.000626322716627787d0
            xpg(202) = 0.169539550700162253d0
            ypg(202) = 0.733998298658581071d0
            zpg(202) = 0.078685547871113893d0
            hpg(202) = 0.000626322716627787d0
            xpg(203) = 0.017776602770142781d0
            ypg(203) = 0.169539550700162253d0
            zpg(203) = 0.078685547871113893d0
            hpg(203) = 0.000626322716627787d0
            xpg(204) = 0.733998298658581071d0
            ypg(204) = 0.169539550700162253d0
            zpg(204) = 0.078685547871113893d0
            hpg(204) = 0.000626322716627787d0

        else if (fapg .eq. 'FPG4NOS') then
! --------- POUR LES POINTS DE GAUSS
            a1 = (5.d0-rac5)/20.d0
            b1 = (5.d0+3.d0*rac5)/20.d0
            h5 = un/24.d0
            npi = 0
            do i = 1, 4
                npi = npi+1
                xpg(npi) = a1
                ypg(npi) = a1
                zpg(npi) = a1
                hpg(npi) = h5
            end do
            zpg(2) = b1
            ypg(3) = b1
            xpg(4) = b1
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+4) = vol/nnos
                xpg(iNode+4) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+4) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+4) = xno(ndim*(iNode-1)+3)
            end do
        end if

    else if (elrefa .eq. 'PY5' .or. elrefa .eq. 'P13' .or. elrefa .eq. 'P19') then
        if (fapg .eq. 'FPG1B') then
! --------- ORDRE 1
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
            xpg(1) = 0.d0
            ypg(1) = 0.d0
            zpg(1) = 0.25d0
            hpg(1) = 0.6666666666666667d0
        else if (fapg .eq. 'FPG5') then
! --------- ORDRE 2
            p1 = deux/15.d0
            h1 = 0.1531754163448146d0
            h2 = 0.6372983346207416d0
            xpg(1) = undemi
            xpg(2) = zero
            xpg(3) = -undemi
            xpg(4) = zero
            xpg(5) = zero
            ypg(1) = zero
            ypg(2) = undemi
            ypg(3) = zero
            ypg(4) = -undemi
            ypg(5) = zero
            zpg(1) = h1
            zpg(2) = h1
            zpg(3) = h1
            zpg(4) = h1
            zpg(5) = h2
            hpg(1) = p1
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p1
            hpg(5) = p1

        else if (fapg .eq. 'FPG6') then
!    Order 3 (REX 33086)
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.

            xpg(1) = 0.d0
            ypg(1) = 0.d0
            zpg(1) = 0.5714285703860683d0
            hpg(1) = 0.1681372559485071d0

            xpg(2) = 0.d0
            ypg(2) = 0.d0
            zpg(2) = 5.67585d-09
            hpg(2) = 0.07500000404404333d0

            xpg(3) = -0.5610836110587396d0
            ypg(3) = 0.d0
            zpg(3) = 0.1666666666666667d0
            hpg(3) = 0.1058823516685291d0

            xpg(4) = 0.d0
            ypg(4) = -0.5610836110587396d0
            zpg(4) = 0.1666666666666667d0
            hpg(4) = 0.1058823516685291d0

            xpg(5) = 0.d0
            ypg(5) = 0.5610836110587396d0
            zpg(5) = 0.1666666666666667d0
            hpg(5) = 0.1058823516685291d0

            xpg(6) = 0.5610836110587396d0
            ypg(6) = 0.d0
            zpg(6) = 0.1666666666666667d0
            hpg(6) = 0.1058823516685291d0

        else if (fapg .eq. 'FPG10') then
! Order 4 (REX 20813)
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.

            xpg(1) = 0.d0
            ypg(1) = 0.d0
            zpg(1) = 0.6772327888861374d0
            hpg(1) = 0.07582792211376127d0

            xpg(2) = 0.d0
            ypg(2) = 0.d0
            zpg(2) = 0.1251369531087465d0
            hpg(2) = 0.1379222683930349d0

            xpg(3) = -0.3252907781991163d0
            ypg(3) = -0.3252907781991163d0
            zpg(3) = 0.3223841495782137d0
            hpg(3) = 0.07088305859288367d0

            xpg(4) = -0.3252907781991163d0
            ypg(4) = 0.3252907781991163d0
            zpg(4) = 0.3223841495782137d0
            hpg(4) = 0.07088305859288367d0

            xpg(5) = 0.3252907781991163d0
            ypg(5) = 0.3252907781991163d0
            zpg(5) = 0.3223841495782137d0
            hpg(5) = 0.07088305859288367d0

            xpg(6) = 0.3252907781991163d0
            ypg(6) = -0.3252907781991163d0
            zpg(6) = 0.3223841495782137d0
            hpg(6) = 0.07088305859288367d0

            xpg(7) = -0.65796699712169d0
            ypg(7) = 0.d0
            zpg(7) = 0.0392482838988154d0
            hpg(7) = 0.04234606044708394d0

            xpg(8) = 0.d0
            ypg(8) = -0.65796699712169d0
            zpg(8) = 0.0392482838988154d0
            hpg(8) = 0.04234606044708394d0

            xpg(9) = 0.d0
            ypg(9) = 0.65796699712169d0
            zpg(9) = 0.0392482838988154d0
            hpg(9) = 0.04234606044708394d0

            xpg(10) = 0.65796699712169d0
            ypg(10) = 0.d0
            zpg(10) = 0.0392482838988154d0
            hpg(10) = 0.04234606044708394d0

        else if (fapg .eq. 'FPG15') then
            ! Order 5
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
            xpg(1) = 0.000000000000000d0
            ypg(1) = 0.000000000000000d0
            zpg(1) = 0.7298578807825067d0
            hpg(1) = 0.04562357993942674d0

            xpg(2) = 0.000000000000000d0
            ypg(2) = 0.000000000000000d0
            zpg(2) = 0.300401020813769d0
            hpg(2) = 0.112931409661816d0

            xpg(3) = 0.000000000000000d0
            ypg(3) = 0.000000000000000d0
            zpg(3) = 6.4917722d-09
            hpg(3) = 0.03913635721904967d0

            xpg(4) = -0.3532630157731623d0
            ypg(4) = -0.3532630157731623d0
            zpg(4) = 0.125d0
            hpg(4) = 0.05096086209874681d0

            xpg(5) = -0.3532630157731623d0
            ypg(5) = 0.3532630157731623d0
            zpg(5) = 0.125d0
            hpg(5) = 0.05096086209874681d0

            xpg(6) = 0.3532630157731623d0
            ypg(6) = 0.3532630157731623d0
            zpg(6) = 0.125d0
            hpg(6) = 0.05096086209874681d0

            xpg(7) = 0.3532630157731623d0
            ypg(7) = -0.3532630157731623d0
            zpg(7) = 0.125d0
            hpg(7) = 0.05096086209874681d0

            xpg(8) = -0.7051171227788277d0
            ypg(8) = 0.000000000000000d0
            zpg(8) = 0.061111907062023d0
            hpg(8) = 0.02644726771976367d0

            xpg(9) = 0.000000000000000d0
            ypg(9) = -0.7051171227788277d0
            zpg(9) = 0.061111907062023d0
            hpg(9) = 0.02644726771976367d0

            xpg(10) = 0.000000000000000d0
            ypg(10) = 0.7051171227788277d0
            zpg(10) = 0.061111907062023d0
            hpg(10) = 0.02644726771976367d0

            xpg(11) = 0.7051171227788277d0
            ypg(11) = 0.000000000000000d0
            zpg(11) = 0.061111907062023d0
            hpg(11) = 0.02644726771976367d0

            xpg(12) = -0.432882864103541d0
            ypg(12) = 0.000000000000000d0
            zpg(12) = 0.4236013371197248d0
            hpg(12) = 0.03983570014308314d0

            xpg(13) = 0.000000000000000d0
            ypg(13) = -0.432882864103541d0
            zpg(13) = 0.4236013371197248d0
            hpg(13) = 0.03983570014308314d0

            xpg(14) = 0.000000000000000d0
            ypg(14) = 0.432882864103541d0
            zpg(14) = 0.4236013371197248d0
            hpg(14) = 0.03983570014308314d0

            xpg(15) = 0.432882864103541d0
            ypg(15) = 0.000000000000000d0
            zpg(15) = 0.4236013371197248d0
            hpg(15) = 0.03983570014308314d0

        else if (fapg .eq. 'FPG24') then
            ! Order 6
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
            xpg(1) = 0.000000000000000d0
            ypg(1) = 0.000000000000000d0
            zpg(1) = 0.8076457976939595d0
            hpg(1) = 0.01697526244176133d0

            xpg(2) = 0.000000000000000d0
            ypg(2) = 0.000000000000000d0
            zpg(2) = 0.0017638088528196d0
            hpg(2) = 0.0107023421167942d0

            xpg(3) = 0.000000000000000d0
            ypg(3) = 0.000000000000000d0
            zpg(3) = 0.1382628064637306d0
            hpg(3) = 0.0797197029683492d0

            xpg(4) = 0.000000000000000d0
            ypg(4) = 0.000000000000000d0
            zpg(4) = 0.4214239119356371d0
            hpg(4) = 0.0687071134661012d0

            xpg(5) = -0.4172976755573542d0
            ypg(5) = -0.4172976755573542d0
            zpg(5) = 0.097447341025462d0
            hpg(5) = 0.02463372529088633d0

            xpg(6) = -0.4172976755573542d0
            ypg(6) = 0.4172976755573542d0
            zpg(6) = 0.097447341025462d0
            hpg(6) = 0.02463372529088633d0

            xpg(7) = 0.4172976755573542d0
            ypg(7) = 0.4172976755573542d0
            zpg(7) = 0.097447341025462d0
            hpg(7) = 0.02463372529088633d0

            xpg(8) = 0.4172976755573542d0
            ypg(8) = -0.4172976755573542d0
            zpg(8) = 0.097447341025462d0
            hpg(8) = 0.02463372529088633d0

            xpg(9) = -0.2169627046883496d0
            ypg(9) = -0.2169627046883496d0
            zpg(9) = 0.5660745906233009d0
            hpg(9) = 0.02105838632544886d0

            xpg(10) = -0.2169627046883496d0
            ypg(10) = 0.2169627046883496d0
            zpg(10) = 0.5660745906233009d0
            hpg(10) = 0.02105838632544886d0

            xpg(11) = 0.2169627046883496d0
            ypg(11) = 0.2169627046883496d0
            zpg(11) = 0.5660745906233009d0
            hpg(11) = 0.02105838632544886d0

            xpg(12) = 0.2169627046883496d0
            ypg(12) = -0.2169627046883496d0
            zpg(12) = 0.5660745906233009d0
            hpg(12) = 0.02105838632544886d0

            xpg(13) = -0.5656808544256755d0
            ypg(13) = 0.000000000000000d0
            zpg(13) = 0.0294777308457207d0
            hpg(13) = 0.0248000862596322d0

            xpg(14) = 0.000000000000000d0
            ypg(14) = -0.5656808544256755d0
            zpg(14) = 0.0294777308457207d0
            hpg(14) = 0.0248000862596322d0

            xpg(15) = 0.000000000000000d0
            ypg(15) = 0.5656808544256755d0
            zpg(15) = 0.0294777308457207d0
            hpg(15) = 0.0248000862596322d0

            xpg(16) = 0.5656808544256755d0
            ypg(16) = 0.000000000000000d0
            zpg(16) = 0.0294777308457207d0
            hpg(16) = 0.0248000862596322d0

            xpg(17) = -0.4980790917807059d0
            ypg(17) = 0.000000000000000d0
            zpg(17) = 0.2649158632121295d0
            hpg(17) = 0.04925492311795127d0

            xpg(18) = 0.000000000000000d0
            ypg(18) = -0.4980790917807059d0
            zpg(18) = 0.2649158632121295d0
            hpg(18) = 0.04925492311795127d0

            xpg(19) = 0.000000000000000d0
            ypg(19) = 0.4980790917807059d0
            zpg(19) = 0.2649158632121295d0
            hpg(19) = 0.04925492311795127d0

            xpg(20) = 0.4980790917807059d0
            ypg(20) = 0.000000000000000d0
            zpg(20) = 0.2649158632121295d0
            hpg(20) = 0.04925492311795127d0

            xpg(21) = -0.9508994872144825d0
            ypg(21) = 0.000000000000000d0
            zpg(21) = 0.048249070631936d0
            hpg(21) = 0.0028934404244966d0

            xpg(22) = 0.000000000000000d0
            ypg(22) = -0.9508994872144825d0
            zpg(22) = 0.048249070631936d0
            hpg(22) = 0.0028934404244966d0

            xpg(23) = 0.000000000000000d0
            ypg(23) = 0.9508994872144825d0
            zpg(23) = 0.048249070631936d0
            hpg(23) = 0.0028934404244966d0

            xpg(24) = 0.9508994872144825d0
            ypg(24) = 0.000000000000000d0
            zpg(24) = 0.048249070631936d0
            hpg(24) = 0.0028934404244966d0
        else if (fapg .eq. 'FPG31') then
            ! Order 7
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
            xpg(1) = 0.000000000000000d0
            ypg(1) = 0.000000000000000d0
            zpg(1) = 4.54802082028d-05
            hpg(1) = 0.01666705135613133d0

            xpg(2) = 0.000000000000000d0
            ypg(2) = 0.000000000000000d0
            zpg(2) = 0.3935828767576542d0
            hpg(2) = 0.06702012179992874d0

            xpg(3) = 0.000000000000000d0
            ypg(3) = 0.000000000000000d0
            zpg(3) = 0.8386827828168515d0
            hpg(3) = 0.01047144961801153d0

            xpg(4) = -0.4320493798938573d0
            ypg(4) = -0.4320493798938573d0
            zpg(4) = 0.06666666666666669d0
            hpg(4) = 0.01779450619533527d0

            xpg(5) = -0.4320493798938573d0
            ypg(5) = 0.4320493798938573d0
            zpg(5) = 0.06666666666666669d0
            hpg(5) = 0.01779450619533527d0

            xpg(6) = 0.4320493798938573d0
            ypg(6) = 0.4320493798938573d0
            zpg(6) = 0.06666666666666669d0
            hpg(6) = 0.01779450619533527d0

            xpg(7) = 0.4320493798938573d0
            ypg(7) = -0.4320493798938573d0
            zpg(7) = 0.06666666666666669d0
            hpg(7) = 0.01779450619533527d0

            xpg(8) = -0.3086066999241838d0
            ypg(8) = -0.3086066999241838d0
            zpg(8) = 0.3333333333333334d0
            hpg(8) = 0.019140625d0

            xpg(9) = -0.3086066999241838d0
            ypg(9) = 0.3086066999241838d0
            zpg(9) = 0.3333333333333334d0
            hpg(9) = 0.019140625d0

            xpg(10) = 0.3086066999241838d0
            ypg(10) = 0.3086066999241838d0
            zpg(10) = 0.3333333333333334d0
            hpg(10) = 0.019140625d0

            xpg(11) = 0.3086066999241838d0
            ypg(11) = -0.3086066999241838d0
            zpg(11) = 0.3333333333333334d0
            hpg(11) = 0.019140625d0

            xpg(12) = -0.3541523808681161d0
            ypg(12) = 0.000000000000000d0
            zpg(12) = 0.1292864095395495d0
            hpg(12) = 0.04397345809998487d0

            xpg(13) = 0.000000000000000d0
            ypg(13) = -0.3541523808681161d0
            zpg(13) = 0.1292864095395495d0
            hpg(13) = 0.04397345809998487d0

            xpg(14) = 0.000000000000000d0
            ypg(14) = 0.3541523808681161d0
            zpg(14) = 0.1292864095395495d0
            hpg(14) = 0.04397345809998487d0

            xpg(15) = 0.3541523808681161d0
            ypg(15) = 0.000000000000000d0
            zpg(15) = 0.1292864095395495d0
            hpg(15) = 0.04397345809998487d0

            xpg(16) = -0.8027258501000878d0
            ypg(16) = 0.000000000000000d0
            zpg(16) = 0.08013853017977809d0
            hpg(16) = 0.009857425304157133d0

            xpg(17) = 0.000000000000000d0
            ypg(17) = -0.8027258501000878d0
            zpg(17) = 0.08013853017977809d0
            hpg(17) = 0.009857425304157133d0

            xpg(18) = 0.000000000000000d0
            ypg(18) = 0.8027258501000878d0
            zpg(18) = 0.08013853017977809d0
            hpg(18) = 0.009857425304157133d0

            xpg(19) = 0.8027258501000878d0
            ypg(19) = 0.000000000000000d0
            zpg(19) = 0.08013853017977809d0
            hpg(19) = 0.009857425304157133d0

            xpg(20) = -0.2541353468618572d0
            ypg(20) = 0.000000000000000d0
            zpg(20) = 0.6055199301105999d0
            hpg(20) = 0.01967940646281593d0

            xpg(21) = 0.000000000000000d0
            ypg(21) = -0.2541353468618572d0
            zpg(21) = 0.6055199301105999d0
            hpg(21) = 0.01967940646281593d0

            xpg(22) = 0.000000000000000d0
            ypg(22) = 0.2541353468618572d0
            zpg(22) = 0.6055199301105999d0
            hpg(22) = 0.01967940646281593d0

            xpg(23) = 0.2541353468618572d0
            ypg(23) = 0.000000000000000d0
            zpg(23) = 0.6055199301105999d0
            hpg(23) = 0.01967940646281593d0

            xpg(24) = -0.6143051077207853d0
            ypg(24) = 0.000000000000000d0
            zpg(24) = 8.97583d-11
            hpg(24) = 0.008869229979749401d0

            xpg(25) = 0.000000000000000d0
            ypg(25) = -0.6143051077207853d0
            zpg(25) = 8.97583d-11
            hpg(25) = 0.008869229979749401d0

            xpg(26) = 0.000000000000000d0
            ypg(26) = 0.6143051077207853d0
            zpg(26) = 8.97583d-11
            hpg(26) = 0.008869229979749401d0

            xpg(27) = 0.6143051077207853d0
            ypg(27) = 0.000000000000000d0
            zpg(27) = 8.97583d-11
            hpg(27) = 0.008869229979749401d0

            xpg(28) = -0.5248326982543755d0
            ypg(28) = 0.000000000000000d0
            zpg(28) = 0.2905430754945976d0
            hpg(28) = 0.0238123599311062d0

            xpg(29) = 0.000000000000000d0
            ypg(29) = -0.5248326982543755d0
            zpg(29) = 0.2905430754945976d0
            hpg(29) = 0.0238123599311062d0

            xpg(30) = 0.000000000000000d0
            ypg(30) = 0.5248326982543755d0
            zpg(30) = 0.2905430754945976d0
            hpg(30) = 0.0238123599311062d0

            xpg(31) = 0.5248326982543755d0
            ypg(31) = 0.000000000000000d0
            zpg(31) = 0.2905430754945976d0
            hpg(31) = 0.0238123599311062d0
        else if (fapg .eq. 'FPG47') then
            ! Order 8
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
            xpg(1) = 0.d0
            ypg(1) = 0.d0
            zpg(1) = 0.584805434108126d0
            hpg(1) = 0.03296066343718627d0

            xpg(2) = 0.d0
            ypg(2) = 0.d0
            zpg(2) = 0.2714503954191895d0
            hpg(2) = 0.05823645962743066d0

            xpg(3) = 0.d0
            ypg(3) = 0.d0
            zpg(3) = 0.071440430215454d0
            hpg(3) = 0.02564038499790467d0

            xpg(4) = -0.3980480443871291d0
            ypg(4) = -0.3980480443871291d0
            zpg(4) = 7.323639d-9
            hpg(4) = 0.006250706266751667d0

            xpg(5) = -0.3980480443871291d0
            ypg(5) = 0.3980480443871291d0
            zpg(5) = 7.323639d-9
            hpg(5) = 0.006250706266751667d0

            xpg(6) = 0.3980480443871291d0
            ypg(6) = 0.3980480443871291d0
            zpg(6) = 7.323639d-9
            hpg(6) = 0.006250706266751667d0

            xpg(7) = 0.3980480443871291d0
            ypg(7) = -0.3980480443871291d0
            zpg(7) = 7.323639d-9
            hpg(7) = 0.006250706266751667d0

            xpg(8) = -0.2824444793709199d0
            ypg(8) = -0.2824444793709199d0
            zpg(8) = 0.4221837006010586d0
            hpg(8) = 0.01087401281095773d0

            xpg(9) = -0.2824444793709199d0
            ypg(9) = 0.2824444793709199d0
            zpg(9) = 0.4221837006010586d0
            hpg(9) = 0.01087401281095773d0

            xpg(10) = 0.2824444793709199d0
            ypg(10) = 0.2824444793709199d0
            zpg(10) = 0.4221837006010586d0
            hpg(10) = 0.01087401281095773d0

            xpg(11) = 0.2824444793709199d0
            ypg(11) = -0.2824444793709199d0
            zpg(11) = 0.4221837006010586d0
            hpg(11) = 0.01087401281095773d0

            xpg(12) = -0.3765113467744738d0
            ypg(12) = -0.3765113467744738d0
            zpg(12) = 0.1477639972965497d0
            hpg(12) = 0.0214931480473806d0

            xpg(13) = -0.3765113467744738d0
            ypg(13) = 0.3765113467744738d0
            zpg(13) = 0.1477639972965497d0
            hpg(13) = 0.0214931480473806d0

            xpg(14) = 0.3765113467744738d0
            ypg(14) = 0.3765113467744738d0
            zpg(14) = 0.1477639972965497d0
            hpg(14) = 0.0214931480473806d0

            xpg(15) = 0.3765113467744738d0
            ypg(15) = -0.3765113467744738d0
            zpg(15) = 0.1477639972965497d0
            hpg(15) = 0.0214931480473806d0

            xpg(16) = -0.246577942462287d0
            ypg(16) = 0.d0
            zpg(16) = 0.6734528479254847d0
            hpg(16) = 0.008423399747923534d0

            xpg(17) = 0.d0
            ypg(17) = -0.246577942462287d0
            zpg(17) = 0.6734528479254847d0
            hpg(17) = 0.008423399747923534d0

            xpg(18) = 0.d0
            ypg(18) = 0.246577942462287d0
            zpg(18) = 0.6734528479254847d0
            hpg(18) = 0.008423399747923534d0

            xpg(19) = 0.246577942462287d0
            ypg(19) = 0.d0
            zpg(19) = 0.6734528479254847d0
            hpg(19) = 0.008423399747923534d0

            xpg(20) = -0.7338895325467361d0
            ypg(20) = 0.d0
            zpg(20) = 0.2455708050362506d0
            hpg(20) = 0.004143696684076733d0

            xpg(21) = 0.d0
            ypg(21) = -0.7338895325467361d0
            zpg(21) = 0.2455708050362506d0
            hpg(21) = 0.004143696684076733d0

            xpg(22) = 0.d0
            ypg(22) = 0.7338895325467361d0
            zpg(22) = 0.2455708050362506d0
            hpg(22) = 0.004143696684076733d0

            xpg(23) = 0.7338895325467361d0
            ypg(23) = 0.d0
            zpg(23) = 0.2455708050362506d0
            hpg(23) = 0.004143696684076733d0

            xpg(24) = -8.26826715128d-5
            ypg(24) = 0.d0
            zpg(24) = 0.8582503626505816d0
            hpg(24) = 0.0017421685557662d0

            xpg(25) = 0.d0
            ypg(25) = -8.26826715128d-5
            zpg(25) = 0.8582503626505816d0
            hpg(25) = 0.0017421685557662d0

            xpg(26) = 0.d0
            ypg(26) = 8.26826715128d-5
            zpg(26) = 0.8582503626505816d0
            hpg(26) = 0.0017421685557662d0

            xpg(27) = 8.26826715128d-5
            ypg(27) = 0.d0
            zpg(27) = 0.8582503626505816d0
            hpg(27) = 0.0017421685557662d0

            xpg(28) = -0.3989965570952277d0
            ypg(28) = 0.d0
            zpg(28) = 0.0474032316194875d0
            hpg(28) = 0.019549978796694d0

            xpg(29) = 0.d0
            ypg(29) = -0.3989965570952277d0
            zpg(29) = 0.0474032316194875d0
            hpg(29) = 0.019549978796694d0

            xpg(30) = 0.d0
            ypg(30) = 0.3989965570952277d0
            zpg(30) = 0.0474032316194875d0
            hpg(30) = 0.019549978796694d0

            xpg(31) = 0.3989965570952277d0
            ypg(31) = 0.d0
            zpg(31) = 0.0474032316194875d0
            hpg(31) = 0.019549978796694d0

            xpg(32) = -0.36531241447959d0
            ypg(32) = 0.d0
            zpg(32) = 0.4021857062205061d0
            hpg(32) = 0.0242395895094498d0

            xpg(33) = 0.d0
            ypg(33) = -0.36531241447959d0
            zpg(33) = 0.4021857062205061d0
            hpg(33) = 0.0242395895094498d0

            xpg(34) = 0.d0
            ypg(34) = 0.36531241447959d0
            zpg(34) = 0.4021857062205061d0
            hpg(34) = 0.0242395895094498d0

            xpg(35) = 0.36531241447959d0
            ypg(35) = 0.d0
            zpg(35) = 0.4021857062205061d0
            hpg(35) = 0.0242395895094498d0

            xpg(36) = -0.4834252440408382d0
            ypg(36) = 0.d0
            zpg(36) = 0.1903148210864091d0
            hpg(36) = 0.02732137396394747d0

            xpg(37) = 0.d0
            ypg(37) = -0.4834252440408382d0
            zpg(37) = 0.1903148210864091d0
            hpg(37) = 0.02732137396394747d0

            xpg(38) = 0.d0
            ypg(38) = 0.4834252440408382d0
            zpg(38) = 0.1903148210864091d0
            hpg(38) = 0.02732137396394747d0

            xpg(39) = 0.4834252440408382d0
            ypg(39) = 0.d0
            zpg(39) = 0.1903148210864091d0
            hpg(39) = 0.02732137396394747d0

            xpg(40) = -0.7631304970460457d0
            ypg(40) = 0.1317256994454848d0
            zpg(40) = 0.0454758042237511d0
            hpg(40) = 0.006709607634044334d0

            xpg(41) = -0.7631304970460457d0
            ypg(41) = -0.1317256994454848d0
            zpg(41) = 0.0454758042237511d0
            hpg(41) = 0.006709607634044334d0

            xpg(42) = 0.1317256994454848d0
            ypg(42) = -0.7631304970460457d0
            zpg(42) = 0.0454758042237511d0
            hpg(42) = 0.006709607634044334d0

            xpg(43) = 0.1317256994454848d0
            ypg(43) = 0.7631304970460457d0
            zpg(43) = 0.0454758042237511d0
            hpg(43) = 0.006709607634044334d0

            xpg(44) = -0.1317256994454848d0
            ypg(44) = 0.7631304970460457d0
            zpg(44) = 0.0454758042237511d0
            hpg(44) = 0.006709607634044334d0

            xpg(45) = -0.1317256994454848d0
            ypg(45) = -0.7631304970460457d0
            zpg(45) = 0.0454758042237511d0
            hpg(45) = 0.006709607634044334d0

            xpg(46) = 0.7631304970460457d0
            ypg(46) = -0.1317256994454848d0
            zpg(46) = 0.0454758042237511d0
            hpg(46) = 0.006709607634044334d0

            xpg(47) = 0.7631304970460457d0
            ypg(47) = 0.1317256994454848d0
            zpg(47) = 0.0454758042237511d0
            hpg(47) = 0.006709607634044334d0
!
!
        else if (fapg .eq. 'FPG62') then
            ! Order 9
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
            xpg(1) = 0.d0
            ypg(1) = 0.d0
            zpg(1) = 0.0371459309315055d0
            hpg(1) = 0.02136524609022833d0

            xpg(2) = 0.d0
            ypg(2) = 0.d0
            zpg(2) = 0.6942002631703261d0
            hpg(2) = 0.01586214729671873d0

            xpg(3) = -0.4629100479304213d0
            ypg(3) = -0.4629100479304213d0
            zpg(3) = 4.2251285d-9
            hpg(3) = 0.002782248201996267d0

            xpg(4) = -0.4629100479304213d0
            ypg(4) = 0.4629100479304213d0
            zpg(4) = 4.2251285d-9
            hpg(4) = 0.002782248201996267d0

            xpg(5) = 0.4629100479304213d0
            ypg(5) = 0.4629100479304213d0
            zpg(5) = 4.2251285d-9
            hpg(5) = 0.002782248201996267d0

            xpg(6) = 0.4629100479304213d0
            ypg(6) = -0.4629100479304213d0
            zpg(6) = 4.2251285d-9
            hpg(6) = 0.002782248201996267d0

            xpg(7) = -0.3380484960581819d0
            ypg(7) = -0.3380484960581819d0
            zpg(7) = 0.2697317845200571d0
            hpg(7) = 0.01486252757307867d0

            xpg(8) = -0.3380484960581819d0
            ypg(8) = 0.3380484960581819d0
            zpg(8) = 0.2697317845200571d0
            hpg(8) = 0.01486252757307867d0

            xpg(9) = 0.3380484960581819d0
            ypg(9) = 0.3380484960581819d0
            zpg(9) = 0.2697317845200571d0
            hpg(9) = 0.01486252757307867d0

            xpg(10) = 0.3380484960581819d0
            ypg(10) = -0.3380484960581819d0
            zpg(10) = 0.2697317845200571d0
            hpg(10) = 0.01486252757307867d0

            xpg(11) = -0.3081660889254626d0
            ypg(11) = -0.3081660889254626d0
            zpg(11) = 0.0833333333333333d0
            hpg(11) = 0.0194186818519688d0

            xpg(12) = -0.3081660889254626d0
            ypg(12) = 0.3081660889254626d0
            zpg(12) = 0.0833333333333333d0
            hpg(12) = 0.0194186818519688d0

            xpg(13) = 0.3081660889254626d0
            ypg(13) = 0.3081660889254626d0
            zpg(13) = 0.0833333333333333d0
            hpg(13) = 0.0194186818519688d0

            xpg(14) = 0.3081660889254626d0
            ypg(14) = -0.3081660889254626d0
            zpg(14) = 0.0833333333333333d0
            hpg(14) = 0.0194186818519688d0

            xpg(15) = -0.2159855846425902d0
            ypg(15) = -0.2159855846425902d0
            zpg(15) = 0.5334178104457834d0
            hpg(15) = 0.007694537698867333d0

            xpg(16) = -0.2159855846425902d0
            ypg(16) = 0.2159855846425902d0
            zpg(16) = 0.5334178104457834d0
            hpg(16) = 0.007694537698867333d0

            xpg(17) = 0.2159855846425902d0
            ypg(17) = 0.2159855846425902d0
            zpg(17) = 0.5334178104457834d0
            hpg(17) = 0.007694537698867333d0

            xpg(18) = 0.2159855846425902d0
            ypg(18) = -0.2159855846425902d0
            zpg(18) = 0.5334178104457834d0
            hpg(18) = 0.007694537698867333d0

            xpg(19) = -0.2418018136307699d0
            ypg(19) = 0.d0
            zpg(19) = 0.6619880087153507d0
            hpg(19) = 0.006733716551573132d0

            xpg(20) = 0.d0
            ypg(20) = -0.2418018136307699d0
            zpg(20) = 0.6619880087153507d0
            hpg(20) = 0.006733716551573132d0

            xpg(21) = 0.d0
            ypg(21) = 0.2418018136307699d0
            zpg(21) = 0.6619880087153507d0
            hpg(21) = 0.006733716551573132d0

            xpg(22) = 0.2418018136307699d0
            ypg(22) = 0.d0
            zpg(22) = 0.6619880087153507d0
            hpg(22) = 0.006733716551573132d0

            xpg(23) = -0.8440791629031899d0
            ypg(23) = 0.d0
            zpg(23) = 0.1522047098054636d0
            hpg(23) = 0.0013726521342804d0

            xpg(24) = 0.d0
            ypg(24) = -0.8440791629031899d0
            zpg(24) = 0.1522047098054636d0
            hpg(24) = 0.0013726521342804d0

            xpg(25) = 0.d0
            ypg(25) = 0.8440791629031899d0
            zpg(25) = 0.1522047098054636d0
            hpg(25) = 0.0013726521342804d0

            xpg(26) = 0.8440791629031899d0
            ypg(26) = 0.d0
            zpg(26) = 0.1522047098054636d0
            hpg(26) = 0.0013726521342804d0

            xpg(27) = -0.2279669387608253d0
            ypg(27) = 0.d0
            zpg(27) = 0.4222655365894455d0
            hpg(27) = 0.0228304591771452d0

            xpg(28) = 0.d0
            ypg(28) = -0.2279669387608253d0
            zpg(28) = 0.4222655365894455d0
            hpg(28) = 0.0228304591771452d0

            xpg(29) = 0.d0
            ypg(29) = 0.2279669387608253d0
            zpg(29) = 0.4222655365894455d0
            hpg(29) = 0.0228304591771452d0

            xpg(30) = 0.2279669387608253d0
            ypg(30) = 0.d0
            zpg(30) = 0.4222655365894455d0
            hpg(30) = 0.0228304591771452d0

            xpg(31) = -0.5027849771350101d0
            ypg(31) = 0.d0
            zpg(31) = 0.3967643698988902d0
            hpg(31) = 0.008606275761264799d0

            xpg(32) = 0.d0
            ypg(32) = -0.5027849771350101d0
            zpg(32) = 0.3967643698988902d0
            hpg(32) = 0.008606275761264799d0

            xpg(33) = 0.d0
            ypg(33) = 0.5027849771350101d0
            zpg(33) = 0.3967643698988902d0
            hpg(33) = 0.008606275761264799d0

            xpg(34) = 0.5027849771350101d0
            ypg(34) = 0.d0
            zpg(34) = 0.3967643698988902d0
            hpg(34) = 0.008606275761264799d0

            xpg(35) = -0.2605010749834311d0
            ypg(35) = 0.d0
            zpg(35) = 0.2008207219371707d0
            hpg(35) = 0.02474152076522233d0

            xpg(36) = 0.d0
            ypg(36) = -0.2605010749834311d0
            zpg(36) = 0.2008207219371707d0
            hpg(36) = 0.02474152076522233d0

            xpg(37) = 0.d0
            ypg(37) = 0.2605010749834311d0
            zpg(37) = 0.2008207219371707d0
            hpg(37) = 0.02474152076522233d0

            xpg(38) = 0.2605010749834311d0
            ypg(38) = 0.d0
            zpg(38) = 0.2008207219371707d0
            hpg(38) = 0.02474152076522233d0

            xpg(39) = -0.09269587308672431d0
            ypg(39) = 0.d0
            zpg(39) = 0.8661453707796378d0
            hpg(39) = 0.0014791321191674d0

            xpg(40) = 0.d0
            ypg(40) = -0.09269587308672431d0
            zpg(40) = 0.8661453707796378d0
            hpg(40) = 0.0014791321191674d0

            xpg(41) = 0.d0
            ypg(41) = 0.09269587308672431d0
            zpg(41) = 0.8661453707796378d0
            hpg(41) = 0.0014791321191674d0

            xpg(42) = 0.09269587308672431d0
            ypg(42) = 0.d0
            zpg(42) = 0.8661453707796378d0
            hpg(42) = 0.0014791321191674d0

            xpg(43) = -0.4832161680706445d0
            ypg(43) = 0.d0
            zpg(43) = 0.0195040124924115d0
            hpg(43) = 0.0106374727574634d0

            xpg(44) = 0.d0
            ypg(44) = -0.4832161680706445d0
            zpg(44) = 0.0195040124924115d0
            hpg(44) = 0.0106374727574634d0

            xpg(45) = 0.d0
            ypg(45) = 0.4832161680706445d0
            zpg(45) = 0.0195040124924115d0
            hpg(45) = 0.0106374727574634d0

            xpg(46) = 0.4832161680706445d0
            ypg(46) = 0.d0
            zpg(46) = 0.0195040124924115d0
            hpg(46) = 0.0106374727574634d0

            xpg(47) = -0.5671855056082324d0
            ypg(47) = 0.d0
            zpg(47) = 0.2131686017065601d0
            hpg(47) = 0.01537821208488507d0

            xpg(48) = 0.d0
            ypg(48) = -0.5671855056082324d0
            zpg(48) = 0.2131686017065601d0
            hpg(48) = 0.01537821208488507d0

            xpg(49) = 0.d0
            ypg(49) = 0.5671855056082324d0
            zpg(49) = 0.2131686017065601d0
            hpg(49) = 0.01537821208488507d0

            xpg(50) = 0.5671855056082324d0
            ypg(50) = 0.d0
            zpg(50) = 0.2131686017065601d0
            hpg(50) = 0.01537821208488507d0

            xpg(51) = -0.8463949915138761d0
            ypg(51) = 0.d0
            zpg(51) = 0.0153477972480097d0
            hpg(51) = 0.0030662670484654d0

            xpg(52) = 0.d0
            ypg(52) = -0.8463949915138761d0
            zpg(52) = 0.0153477972480097d0
            hpg(52) = 0.0030662670484654d0

            xpg(53) = 0.d0
            ypg(53) = 0.8463949915138761d0
            zpg(53) = 0.0153477972480097d0
            hpg(53) = 0.0030662670484654d0

            xpg(54) = 0.8463949915138761d0
            ypg(54) = 0.d0
            zpg(54) = 0.0153477972480097d0
            hpg(54) = 0.0030662670484654d0

            xpg(55) = -0.6556842905419458d0
            ypg(55) = -0.1914935775757995d0
            zpg(55) = 0.0833333333333333d0
            hpg(55) = 0.008878057297275867d0

            xpg(56) = -0.6556842905419458d0
            ypg(56) = 0.1914935775757995d0
            zpg(56) = 0.0833333333333333d0
            hpg(56) = 0.008878057297275867d0

            xpg(57) = -0.1914935775757995d0
            ypg(57) = -0.6556842905419458d0
            zpg(57) = 0.0833333333333333d0
            hpg(57) = 0.008878057297275867d0

            xpg(58) = -0.1914935775757995d0
            ypg(58) = 0.6556842905419458d0
            zpg(58) = 0.0833333333333333d0
            hpg(58) = 0.008878057297275867d0

            xpg(59) = 0.1914935775757995d0
            ypg(59) = 0.6556842905419458d0
            zpg(59) = 0.0833333333333333d0
            hpg(59) = 0.008878057297275867d0

            xpg(60) = 0.1914935775757995d0
            ypg(60) = -0.6556842905419458d0
            zpg(60) = 0.0833333333333333d0
            hpg(60) = 0.008878057297275867d0

            xpg(61) = 0.6556842905419458d0
            ypg(61) = 0.1914935775757995d0
            zpg(61) = 0.0833333333333333d0
            hpg(61) = 0.008878057297275867d0

            xpg(62) = 0.6556842905419458d0
            ypg(62) = -0.1914935775757995d0
            zpg(62) = 0.0833333333333333d0
            hpg(62) = 0.008878057297275867d0
!
        else if (fapg .eq. 'FPG83') then
            ! Order 10
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!
            xpg(1) = 0.d0
            ypg(1) = 0.d0
            zpg(1) = 0.9167957817791272d0
            hpg(1) = 0.001600817841680933d0

            xpg(2) = 0.d0
            ypg(2) = 0.d0
            zpg(2) = 0.3511368063403526d0
            hpg(2) = 0.0312772852731758d0

            xpg(3) = 0.d0
            ypg(3) = 0.d0
            zpg(3) = 0.7419135679633d0
            hpg(3) = 0.006564440499019467d0

            xpg(4) = -0.4605650780073202d0
            ypg(4) = -0.4605650780073202d0
            zpg(4) = 0.06316159727991449d0
            hpg(4) = 0.0037849282119542d0

            xpg(5) = -0.4605650780073202d0
            ypg(5) = 0.4605650780073202d0
            zpg(5) = 0.06316159727991449d0
            hpg(5) = 0.0037849282119542d0

            xpg(6) = 0.4605650780073202d0
            ypg(6) = 0.4605650780073202d0
            zpg(6) = 0.06316159727991449d0
            hpg(6) = 0.0037849282119542d0

            xpg(7) = 0.4605650780073202d0
            ypg(7) = -0.4605650780073202d0
            zpg(7) = 0.06316159727991449d0
            hpg(7) = 0.0037849282119542d0

            xpg(8) = -0.2273474422126268d0
            ypg(8) = -0.2273474422126268d0
            zpg(8) = 0.1767304111617324d0
            hpg(8) = 0.0274808398278762d0

            xpg(9) = -0.2273474422126268d0
            ypg(9) = 0.2273474422126268d0
            zpg(9) = 0.1767304111617324d0
            hpg(9) = 0.0274808398278762d0

            xpg(10) = 0.2273474422126268d0
            ypg(10) = 0.2273474422126268d0
            zpg(10) = 0.1767304111617324d0
            hpg(10) = 0.0274808398278762d0

            xpg(11) = 0.2273474422126268d0
            ypg(11) = -0.2273474422126268d0
            zpg(11) = 0.1767304111617324d0
            hpg(11) = 0.0274808398278762d0

            xpg(12) = -0.1747734922483533d0
            ypg(12) = -0.1747734922483533d0
            zpg(12) = 0.602093810707914d0
            hpg(12) = 0.006139168845354467d0

            xpg(13) = -0.1747734922483533d0
            ypg(13) = 0.1747734922483533d0
            zpg(13) = 0.602093810707914d0
            hpg(13) = 0.006139168845354467d0

            xpg(14) = 0.1747734922483533d0
            ypg(14) = 0.1747734922483533d0
            zpg(14) = 0.602093810707914d0
            hpg(14) = 0.006139168845354467d0

            xpg(15) = 0.1747734922483533d0
            ypg(15) = -0.1747734922483533d0
            zpg(15) = 0.602093810707914d0
            hpg(15) = 0.006139168845354467d0

            xpg(16) = -0.4036909359989291d0
            ypg(16) = 0.d0
            zpg(16) = 0.5270494697681531d0
            hpg(16) = 0.004712264110046067d0

            xpg(17) = 0.d0
            ypg(17) = -0.4036909359989291d0
            zpg(17) = 0.5270494697681531d0
            hpg(17) = 0.004712264110046067d0

            xpg(18) = 0.d0
            ypg(18) = 0.4036909359989291d0
            zpg(18) = 0.5270494697681531d0
            hpg(18) = 0.004712264110046067d0

            xpg(19) = 0.4036909359989291d0
            ypg(19) = 0.d0
            zpg(19) = 0.5270494697681531d0
            hpg(19) = 0.004712264110046067d0

            xpg(20) = -0.3714608909064742d0
            ypg(20) = 0.d0
            zpg(20) = 0.3491104418985491d0
            hpg(20) = 0.01988712388425d0

            xpg(21) = 0.d0
            ypg(21) = -0.3714608909064742d0
            zpg(21) = 0.3491104418985491d0
            hpg(21) = 0.01988712388425d0

            xpg(22) = 0.d0
            ypg(22) = 0.3714608909064742d0
            zpg(22) = 0.3491104418985491d0
            hpg(22) = 0.01988712388425d0

            xpg(23) = 0.3714608909064742d0
            ypg(23) = 0.d0
            zpg(23) = 0.3491104418985491d0
            hpg(23) = 0.01988712388425d0

            xpg(24) = -0.2995979966434122d0
            ypg(24) = 0.d0
            zpg(24) = 0.0069339728754634d0
            hpg(24) = 0.003913328534544533d0

            xpg(25) = 0.d0
            ypg(25) = -0.2995979966434122d0
            zpg(25) = 0.0069339728754634d0
            hpg(25) = 0.003913328534544533d0

            xpg(26) = 0.d0
            ypg(26) = 0.2995979966434122d0
            zpg(26) = 0.0069339728754634d0
            hpg(26) = 0.003913328534544533d0

            xpg(27) = 0.2995979966434122d0
            ypg(27) = 0.d0
            zpg(27) = 0.0069339728754634d0
            hpg(27) = 0.003913328534544533d0

            xpg(28) = -0.1582262545541785d0
            ypg(28) = 0.d0
            zpg(28) = 0.7752643179331966d0
            hpg(28) = 0.003302910289845933d0

            xpg(29) = 0.d0
            ypg(29) = -0.1582262545541785d0
            zpg(29) = 0.7752643179331966d0
            hpg(29) = 0.003302910289845933d0

            xpg(30) = 0.d0
            ypg(30) = 0.1582262545541785d0
            zpg(30) = 0.7752643179331966d0
            hpg(30) = 0.003302910289845933d0

            xpg(31) = 0.1582262545541785d0
            ypg(31) = 0.d0
            zpg(31) = 0.7752643179331966d0
            hpg(31) = 0.003302910289845933d0

            xpg(32) = -0.1819979349519966d0
            ypg(32) = 0.d0
            zpg(32) = 0.5470650181463391d0
            hpg(32) = 0.01062062484777893d0

            xpg(33) = 0.d0
            ypg(33) = -0.1819979349519966d0
            zpg(33) = 0.5470650181463391d0
            hpg(33) = 0.01062062484777893d0

            xpg(34) = 0.d0
            ypg(34) = 0.1819979349519966d0
            zpg(34) = 0.5470650181463391d0
            hpg(34) = 0.01062062484777893d0

            xpg(35) = 0.1819979349519966d0
            ypg(35) = 0.d0
            zpg(35) = 0.5470650181463391d0
            hpg(35) = 0.01062062484777893d0

            xpg(36) = -0.8455841498012413d0
            ypg(36) = 0.d0
            zpg(36) = 0.07289841147566049d0
            hpg(36) = 0.0033943516567886d0

            xpg(37) = 0.d0
            ypg(37) = -0.8455841498012413d0
            zpg(37) = 0.07289841147566049d0
            hpg(37) = 0.0033943516567886d0

            xpg(38) = 0.d0
            ypg(38) = 0.8455841498012413d0
            zpg(38) = 0.07289841147566049d0
            hpg(38) = 0.0033943516567886d0

            xpg(39) = 0.8455841498012413d0
            ypg(39) = 0.d0
            zpg(39) = 0.07289841147566049d0
            hpg(39) = 0.0033943516567886d0

            xpg(40) = -0.6468694059373429d0
            ypg(40) = 0.d0
            zpg(40) = 0.2836223397977548d0
            hpg(40) = 0.004563862230205867d0

            xpg(41) = 0.d0
            ypg(41) = -0.6468694059373429d0
            zpg(41) = 0.2836223397977548d0
            hpg(41) = 0.004563862230205867d0

            xpg(42) = 0.d0
            ypg(42) = 0.6468694059373429d0
            zpg(42) = 0.2836223397977548d0
            hpg(42) = 0.004563862230205867d0

            xpg(43) = 0.6468694059373429d0
            ypg(43) = 0.d0
            zpg(43) = 0.2836223397977548d0
            hpg(43) = 0.004563862230205867d0

            xpg(44) = -0.5875870294087208d0
            ypg(44) = 0.d0
            zpg(44) = 0.0598130217927265d0
            hpg(44) = 0.01002534692024427d0

            xpg(45) = 0.d0
            ypg(45) = -0.5875870294087208d0
            zpg(45) = 0.0598130217927265d0
            hpg(45) = 0.01002534692024427d0

            xpg(46) = 0.d0
            ypg(46) = 0.5875870294087208d0
            zpg(46) = 0.0598130217927265d0
            hpg(46) = 0.01002534692024427d0

            xpg(47) = 0.5875870294087208d0
            ypg(47) = 0.d0
            zpg(47) = 0.0598130217927265d0
            hpg(47) = 0.01002534692024427d0

            xpg(48) = -0.2237620599838161d0
            ypg(48) = 0.d0
            zpg(48) = 0.0658308799806233d0
            hpg(48) = 0.01041010874919613d0

            xpg(49) = 0.d0
            ypg(49) = -0.2237620599838161d0
            zpg(49) = 0.0658308799806233d0
            hpg(49) = 0.01041010874919613d0

            xpg(50) = 0.d0
            ypg(50) = 0.2237620599838161d0
            zpg(50) = 0.0658308799806233d0
            hpg(50) = 0.01041010874919613d0

            xpg(51) = 0.2237620599838161d0
            ypg(51) = 0.d0
            zpg(51) = 0.0658308799806233d0
            hpg(51) = 0.01041010874919613d0

            xpg(52) = -0.4851415701773136d0
            ypg(52) = -0.2604241382924074d0
            zpg(52) = 0.027250255274674d0
            hpg(52) = 0.005939573377931999d0

            xpg(53) = -0.4851415701773136d0
            ypg(53) = 0.2604241382924074d0
            zpg(53) = 0.027250255274674d0
            hpg(53) = 0.005939573377931999d0

            xpg(54) = -0.2604241382924074d0
            ypg(54) = -0.4851415701773136d0
            zpg(54) = 0.027250255274674d0
            hpg(54) = 0.005939573377931999d0

            xpg(55) = -0.2604241382924074d0
            ypg(55) = 0.4851415701773136d0
            zpg(55) = 0.027250255274674d0
            hpg(55) = 0.005939573377931999d0

            xpg(56) = 0.2604241382924074d0
            ypg(56) = 0.4851415701773136d0
            zpg(56) = 0.027250255274674d0
            hpg(56) = 0.005939573377931999d0

            xpg(57) = 0.2604241382924074d0
            ypg(57) = -0.4851415701773136d0
            zpg(57) = 0.027250255274674d0
            hpg(57) = 0.005939573377931999d0

            xpg(58) = 0.4851415701773136d0
            ypg(58) = 0.2604241382924074d0
            zpg(58) = 0.027250255274674d0
            hpg(58) = 0.005939573377931999d0

            xpg(59) = 0.4851415701773136d0
            ypg(59) = -0.2604241382924074d0
            zpg(59) = 0.027250255274674d0
            hpg(59) = 0.005939573377931999d0

            xpg(60) = -0.5771943577668779d0
            ypg(60) = 0.1715324050350008d0
            zpg(60) = 0.1618000129270508d0
            hpg(60) = 0.01119448097323473d0

            xpg(61) = -0.5771943577668779d0
            ypg(61) = -0.1715324050350008d0
            zpg(61) = 0.1618000129270508d0
            hpg(61) = 0.01119448097323473d0

            xpg(62) = 0.1715324050350008d0
            ypg(62) = -0.5771943577668779d0
            zpg(62) = 0.1618000129270508d0
            hpg(62) = 0.01119448097323473d0

            xpg(63) = 0.1715324050350008d0
            ypg(63) = 0.5771943577668779d0
            zpg(63) = 0.1618000129270508d0
            hpg(63) = 0.01119448097323473d0

            xpg(64) = -0.1715324050350008d0
            ypg(64) = 0.5771943577668779d0
            zpg(64) = 0.1618000129270508d0
            hpg(64) = 0.01119448097323473d0

            xpg(65) = -0.1715324050350008d0
            ypg(65) = -0.5771943577668779d0
            zpg(65) = 0.1618000129270508d0
            hpg(65) = 0.01119448097323473d0

            xpg(66) = 0.5771943577668779d0
            ypg(66) = -0.1715324050350008d0
            zpg(66) = 0.1618000129270508d0
            hpg(66) = 0.01119448097323473d0

            xpg(67) = 0.5771943577668779d0
            ypg(67) = 0.1715324050350008d0
            zpg(67) = 0.1618000129270508d0
            hpg(67) = 0.01119448097323473d0

            xpg(68) = -0.3355997397040416d0
            ypg(68) = 0.264857171422832d0
            zpg(68) = 0.3584705635890507d0
            hpg(68) = 0.005470499403535867d0

            xpg(69) = -0.3355997397040416d0
            ypg(69) = -0.264857171422832d0
            zpg(69) = 0.3584705635890507d0
            hpg(69) = 0.005470499403535867d0

            xpg(70) = 0.264857171422832d0
            ypg(70) = -0.3355997397040416d0
            zpg(70) = 0.3584705635890507d0
            hpg(70) = 0.005470499403535867d0

            xpg(71) = 0.264857171422832d0
            ypg(71) = 0.3355997397040416d0
            zpg(71) = 0.3584705635890507d0
            hpg(71) = 0.005470499403535867d0

            xpg(72) = -0.264857171422832d0
            ypg(72) = 0.3355997397040416d0
            zpg(72) = 0.3584705635890507d0
            hpg(72) = 0.005470499403535867d0

            xpg(73) = -0.264857171422832d0
            ypg(73) = -0.3355997397040416d0
            zpg(73) = 0.3584705635890507d0
            hpg(73) = 0.005470499403535867d0

            xpg(74) = 0.3355997397040416d0
            ypg(74) = -0.264857171422832d0
            zpg(74) = 0.3584705635890507d0
            hpg(74) = 0.005470499403535867d0

            xpg(75) = 0.3355997397040416d0
            ypg(75) = 0.264857171422832d0
            zpg(75) = 0.3584705635890507d0
            hpg(75) = 0.005470499403535867d0

            xpg(76) = -0.7898911255791035d0
            ypg(76) = -0.1523280215005647d0
            zpg(76) = 0.0074406410768847d0
            hpg(76) = 0.0016810325728536d0

            xpg(77) = -0.7898911255791035d0
            ypg(77) = 0.1523280215005647d0
            zpg(77) = 0.0074406410768847d0
            hpg(77) = 0.0016810325728536d0

            xpg(78) = -0.1523280215005647d0
            ypg(78) = -0.7898911255791035d0
            zpg(78) = 0.0074406410768847d0
            hpg(78) = 0.0016810325728536d0

            xpg(79) = -0.1523280215005647d0
            ypg(79) = 0.7898911255791035d0
            zpg(79) = 0.0074406410768847d0
            hpg(79) = 0.0016810325728536d0

            xpg(80) = 0.1523280215005647d0
            ypg(80) = 0.7898911255791035d0
            zpg(80) = 0.0074406410768847d0
            hpg(80) = 0.0016810325728536d0

            xpg(81) = 0.1523280215005647d0
            ypg(81) = -0.7898911255791035d0
            zpg(81) = 0.0074406410768847d0
            hpg(81) = 0.0016810325728536d0

            xpg(82) = 0.7898911255791035d0
            ypg(82) = 0.1523280215005647d0
            zpg(82) = 0.0074406410768847d0
            hpg(82) = 0.0016810325728536d0

            xpg(83) = 0.7898911255791035d0
            ypg(83) = -0.1523280215005647d0
            zpg(83) = 0.0074406410768847d0
            hpg(83) = 0.0016810325728536d0
!
        else if (fapg .eq. 'FPG5NOS') then
! --------- POUR LES POINTS DE GAUSS
            p1 = deux/15.d0
            h1 = 0.1531754163448146d0
            h2 = 0.6372983346207416d0
            xpg(1) = undemi
            xpg(2) = zero
            xpg(3) = -undemi
            xpg(4) = zero
            xpg(5) = zero
            ypg(1) = zero
            ypg(2) = undemi
            ypg(3) = zero
            ypg(4) = -undemi
            ypg(5) = zero
            zpg(1) = h1
            zpg(2) = h1
            zpg(3) = h1
            zpg(4) = h1
            zpg(5) = h2
            hpg(1) = p1
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p1
            hpg(5) = p1
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+5) = vol/nnos
                xpg(iNode+5) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+5) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+5) = xno(ndim*(iNode-1)+3)
            end do
        else
            ASSERT(ASTER_FALSE)

        end if

    else if (elrefa .eq. 'TR3' .or. elrefa .eq. 'TR6' .or. elrefa .eq. 'TR7' &
             .or. elrefa .eq. 'TR1') then
        if (fapg .eq. 'FPG1') then
            xpg(1) = un/3.d0
            ypg(1) = un/3.d0
            hpg(1) = un/deux

        else if (fapg .eq. 'FPG3') then
            xpg(1) = un/6.d0
            ypg(1) = un/6.d0
            xpg(2) = deux/3.d0
            ypg(2) = un/6.d0
            xpg(3) = un/6.d0
            ypg(3) = deux/3.d0
            hpg(1) = un/6.d0
            hpg(2) = un/6.d0
            hpg(3) = un/6.d0

        else if (fapg .eq. 'FPG4') then
            xpg(1) = 0.2d0
            ypg(1) = 0.2d0
            xpg(2) = 0.6d0
            ypg(2) = 0.2d0
            xpg(3) = 0.2d0
            ypg(3) = 0.6d0
            xpg(4) = un/3.d0
            ypg(4) = un/3.d0
            hpg(1) = 25.d0/96.d0
            hpg(2) = 25.d0/96.d0
            hpg(3) = 25.d0/96.d0
            hpg(4) = -27.d0/96.d0

        else if (fapg .eq. 'FPG6') then
            h1 = 0.111690794839005d0
            h2 = 0.054975871827661d0
            a1 = 0.445948490915965d0
            b1 = 0.091576213509771d0
            xpg(3) = b1
            ypg(3) = un-deux*b1
            xpg(1) = b1
            ypg(1) = b1
            xpg(2) = un-deux*b1
            ypg(2) = b1
            xpg(6) = un-deux*a1
            ypg(6) = a1
            xpg(4) = a1
            ypg(4) = un-deux*a1
            xpg(5) = a1
            ypg(5) = a1
            hpg(1) = h2
            hpg(2) = h2
            hpg(3) = h2
            hpg(4) = h1
            hpg(5) = h1
            hpg(6) = h1

        else if (fapg .eq. 'FPG7') then
            p1 = (155.d0+rac15)/2400.d0
            p2 = (155.d0-rac15)/2400.d0
            a2 = (6.d0+rac15)/21.d0
            b2 = (6.d0-rac15)/21.d0
            xpg(1) = un/3.d0
            ypg(1) = un/3.d0
            xpg(2) = a2
            ypg(2) = a2
            xpg(3) = un-deux*a2
            ypg(3) = a2
            xpg(4) = a2
            ypg(4) = un-deux*a2
            xpg(5) = b2
            ypg(5) = b2
            xpg(6) = un-deux*b2
            ypg(6) = b2
            xpg(7) = b2
            ypg(7) = un-deux*b2
            hpg(1) = 9.d0/80.d0
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p1
            hpg(5) = p2
            hpg(6) = p2
            hpg(7) = p2

        else if (fapg .eq. 'FPG12') then
            a1 = 0.063089014491502d0
            b1 = 0.249286745170910d0
            c1 = 0.310352451033785d0
            d1 = 0.053145049844816d0
            xpg(1) = a1
            ypg(1) = a1
            xpg(2) = un-deux*a1
            ypg(2) = a1
            xpg(3) = a1
            ypg(3) = un-deux*a1
            xpg(4) = b1
            ypg(4) = b1
            xpg(5) = un-deux*b1
            ypg(5) = b1
            xpg(6) = b1
            ypg(6) = un-deux*b1
            xpg(7) = c1
            ypg(7) = d1
            xpg(8) = d1
            ypg(8) = c1
            xpg(9) = un-c1-d1
            ypg(9) = c1
            xpg(10) = un-c1-d1
            ypg(10) = d1
            xpg(11) = c1
            ypg(11) = un-c1-d1
            xpg(12) = d1
            ypg(12) = un-c1-d1
            p1 = 0.025422453185103d0
            p2 = 0.058393137863189d0
            p3 = 0.041425537809187d0
            hpg(1) = p1
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p2
            hpg(5) = p2
            hpg(6) = p2
            hpg(7) = p3
            hpg(8) = p3
            hpg(9) = p3
            hpg(10) = p3
            hpg(11) = p3
            hpg(12) = p3

        else if (fapg .eq. 'FPG13') then
!         FORMULE A 13 POINTS : ORDRE 7  (CF BATHE, PAGE 280)
            xpg(1) = 0.0651301029022d0
            ypg(1) = 0.0651301029022d0
            xpg(2) = 0.8697397941956d0
            ypg(2) = 0.0651301029022d0
            xpg(3) = 0.0651301029022d0
            ypg(3) = 0.8697397941956d0
            xpg(4) = 0.3128654960049d0
            ypg(4) = 0.0486903154253d0
            xpg(5) = 0.6384441885698d0
            ypg(5) = 0.3128654960049d0
            xpg(6) = 0.0486903154253d0
            ypg(6) = 0.6384441885698d0
            xpg(7) = 0.6384441885698d0
            ypg(7) = 0.0486903154253d0
            xpg(8) = 0.3128654960049d0
            ypg(8) = 0.6384441885698d0
            xpg(9) = 0.0486903154253d0
            ypg(9) = 0.3128654960049d0
            xpg(10) = 0.2603459660790d0
            ypg(10) = 0.2603459660790d0
            xpg(11) = 0.4793080678419d0
            ypg(11) = 0.2603459660790d0
            xpg(12) = 0.2603459660790d0
            ypg(12) = 0.4793080678419d0
            xpg(13) = 0.3333333333333d0
            ypg(13) = 0.3333333333333d0
            p1 = 0.0533472356088d0/deux
            p2 = 0.0771137608903d0/deux
            p3 = 0.1756152574332d0/deux
            p4 = -0.1495700444677d0/deux
            hpg(1) = p1
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p2
            hpg(5) = p2
            hpg(6) = p2
            hpg(7) = p2
            hpg(8) = p2
            hpg(9) = p2
            hpg(10) = p3
            hpg(11) = p3
            hpg(12) = p3
            hpg(13) = p4

        else if (fapg .eq. 'FPG16') then
            xpg(1) = un/3.d0
            ypg(1) = un/3.d0
            xpg(2) = 0.081414823414554d0
            ypg(2) = 0.459292588292723d0
            xpg(3) = 0.459292588292723d0
            ypg(3) = 0.081414823414554d0
            xpg(4) = 0.459292588292723d0
            ypg(4) = 0.459292588292723d0
            xpg(5) = 0.658861384496480d0
            ypg(5) = 0.170569307751760d0
            xpg(6) = 0.170569307751760d0
            ypg(6) = 0.658861384496480d0
            xpg(7) = 0.170569307751760d0
            ypg(7) = 0.170569307751760d0
            xpg(8) = 0.898905543365938d0
            ypg(8) = 0.050547228317031d0
            xpg(9) = 0.050547228317031d0
            ypg(9) = 0.898905543365938d0
            xpg(10) = 0.050547228317031d0
            ypg(10) = 0.050547228317031d0
            xpg(11) = 0.008394777409958d0
            ypg(11) = 0.728492392955404d0
            xpg(12) = 0.728492392955404d0
            ypg(12) = 0.008394777409958d0
            xpg(13) = 0.263112829634638d0
            ypg(13) = 0.008394777409958d0
            xpg(14) = 0.008394777409958d0
            ypg(14) = 0.263112829634638d0
            xpg(15) = 0.263112829634638d0
            ypg(15) = 0.728492392955404d0
            xpg(16) = 0.728492392955404d0
            ypg(16) = 0.263112829634638d0
            p1 = 0.144315607677787d0/deux
            p2 = 0.095091634267285d0/deux
            p3 = 0.103217370534718d0/deux
            p4 = 0.032458497623198d0/deux
            p5 = 0.027230314174435d0/deux
            hpg(1) = p1
            hpg(2) = p2
            hpg(3) = p2
            hpg(4) = p2
            hpg(5) = p3
            hpg(6) = p3
            hpg(7) = p3
            hpg(8) = p4
            hpg(9) = p4
            hpg(10) = p4
            hpg(11) = p5
            hpg(12) = p5
            hpg(13) = p5
            hpg(14) = p5
            hpg(15) = p5
            hpg(16) = p5

        else if (fapg .eq. 'FPG19') then
! --------- Order 9 : Very high-order symmetric positive-interior
!                      quadrature rules on triangles and tetrahedra
!                      Zelalem Arega Worku  and Jason E. Hicken
!           https://github.com/OptimalDesignLab/SummationByParts.jl
            xpg(1) = 0.0447295133944527098d0
            ypg(1) = 0.9105409732110945803d0
            hpg(1) = 0.0127888378293490156d0
            xpg(2) = 0.0447295133944527098d0
            ypg(2) = 0.0447295133944527098d0
            hpg(2) = 0.0127888378293490156d0
            xpg(3) = 0.9105409732110945803d0
            ypg(3) = 0.0447295133944527098d0
            hpg(3) = 0.0127888378293490156d0
            xpg(4) = 0.1882035356190327302d0
            ypg(4) = 0.6235929287619345395d0
            hpg(4) = 0.0398238694636051265d0
            xpg(5) = 0.1882035356190327302d0
            ypg(5) = 0.1882035356190327302d0
            hpg(5) = 0.0398238694636051265d0
            xpg(6) = 0.6235929287619345395d0
            ypg(6) = 0.1882035356190327302d0
            hpg(6) = 0.0398238694636051265d0
            xpg(7) = 0.4896825191987376277d0
            ypg(7) = 0.0206349616025247444d0
            hpg(7) = 0.0156673501135695352d0
            xpg(8) = 0.4896825191987376277d0
            ypg(8) = 0.4896825191987376277d0
            hpg(8) = 0.0156673501135695352d0
            xpg(9) = 0.0206349616025247444d0
            ypg(9) = 0.4896825191987376277d0
            hpg(9) = 0.0156673501135695352d0
            xpg(10) = 0.4370895914929366372d0
            ypg(10) = 0.1258208170141267254d0
            hpg(10) = 0.0389137705023871396d0
            xpg(11) = 0.4370895914929366372d0
            ypg(11) = 0.4370895914929366372d0
            hpg(11) = 0.0389137705023871396d0
            xpg(12) = 0.1258208170141267254d0
            ypg(12) = 0.4370895914929366372d0
            hpg(12) = 0.0389137705023871396d0
            xpg(13) = 0.0368384120547362836d0
            ypg(13) = 0.2219629891607656956d0
            hpg(13) = 0.0216417696886446886d0
            xpg(14) = 0.7411985987844980207d0
            ypg(14) = 0.2219629891607656956d0
            hpg(14) = 0.0216417696886446886d0
            xpg(15) = 0.7411985987844980207d0
            ypg(15) = 0.0368384120547362836d0
            hpg(15) = 0.0216417696886446886d0
            xpg(16) = 0.0368384120547362836d0
            ypg(16) = 0.7411985987844980207d0
            hpg(16) = 0.0216417696886446886d0
            xpg(17) = 0.2219629891607656956d0
            ypg(17) = 0.7411985987844980207d0
            hpg(17) = 0.0216417696886446886d0
            xpg(18) = 0.2219629891607656956d0
            ypg(18) = 0.0368384120547362836d0
            hpg(18) = 0.0216417696886446886d0
            xpg(19) = 0.3333333333333333333d0
            ypg(19) = 0.3333333333333333333d0
            hpg(19) = 0.0485678981413994169d0

        else if (fapg .eq. 'FPG25') then
! --------- Order 10 : Very high-order symmetric positive-interior
!                      quadrature rules on triangles and tetrahedra
!                      Zelalem Arega Worku  and Jason E. Hicken
!           https://github.com/OptimalDesignLab/SummationByParts.jl
            xpg(1) = 0.0077803194517149288d0
            ypg(1) = 0.9844393610965701424d0
            hpg(1) = 0.0016159489092647598d0
            xpg(2) = 0.0077803194517149288d0
            ypg(2) = 0.0077803194517149288d0
            hpg(2) = 0.0016159489092647598d0
            xpg(3) = 0.9844393610965701424d0
            ypg(3) = 0.0077803194517149288d0
            hpg(3) = 0.0016159489092647598d0
            xpg(4) = 0.1821457744859497251d0
            ypg(4) = 0.6357084510281005498d0
            hpg(4) = 0.0391537997959688796d0
            xpg(5) = 0.1821457744859497251d0
            ypg(5) = 0.1821457744859497251d0
            hpg(5) = 0.0391537997959688796d0
            xpg(6) = 0.6357084510281005498d0
            ypg(6) = 0.1821457744859497251d0
            hpg(6) = 0.0391537997959688796d0
            xpg(7) = 0.4890421285464169384d0
            ypg(7) = 0.0219157429071661231d0
            hpg(7) = 0.0103181844086345647d0
            xpg(8) = 0.4890421285464169384d0
            ypg(8) = 0.4890421285464169384d0
            hpg(8) = 0.0103181844086345647d0
            xpg(9) = 0.0219157429071661231d0
            ypg(9) = 0.4890421285464169384d0
            hpg(9) = 0.0103181844086345647d0
            xpg(10) = 0.4261065018332104709d0
            ypg(10) = 0.1477869963335790581d0
            hpg(10) = 0.0391535362251084783d0
            xpg(11) = 0.4261065018332104709d0
            ypg(11) = 0.4261065018332104709d0
            hpg(11) = 0.0391535362251084783d0
            xpg(12) = 0.1477869963335790581d0
            ypg(12) = 0.4261065018332104709d0
            hpg(12) = 0.0391535362251084783d0
            xpg(13) = 0.0324596699274846068d0
            ypg(13) = 0.1171324896904364880d0
            hpg(13) = 0.0129818721495140074d0
            xpg(14) = 0.8504078403820789051d0
            ypg(14) = 0.1171324896904364880d0
            hpg(14) = 0.0129818721495140074d0
            xpg(15) = 0.8504078403820789051d0
            ypg(15) = 0.0324596699274846068d0
            hpg(15) = 0.0129818721495140074d0
            xpg(16) = 0.0324596699274846068d0
            ypg(16) = 0.8504078403820789051d0
            hpg(16) = 0.0129818721495140074d0
            xpg(17) = 0.1171324896904364880d0
            ypg(17) = 0.8504078403820789051d0
            hpg(17) = 0.0129818721495140074d0
            xpg(18) = 0.1171324896904364880d0
            ypg(18) = 0.0324596699274846068d0
            hpg(18) = 0.0129818721495140074d0
            xpg(19) = 0.0377228854559406535d0
            ypg(19) = 0.3004402842074705852d0
            hpg(19) = 0.0187177221401866412d0
            xpg(20) = 0.6618368303365887612d0
            ypg(20) = 0.3004402842074705852d0
            hpg(20) = 0.0187177221401866412d0
            xpg(21) = 0.6618368303365887612d0
            ypg(21) = 0.0377228854559406535d0
            hpg(21) = 0.0187177221401866412d0
            xpg(22) = 0.0377228854559406535d0
            ypg(22) = 0.6618368303365887612d0
            hpg(22) = 0.0187177221401866412d0
            xpg(23) = 0.3004402842074705852d0
            ypg(23) = 0.6618368303365887612d0
            hpg(23) = 0.0187177221401866412d0
            xpg(24) = 0.3004402842074705852d0
            ypg(24) = 0.0377228854559406535d0
            hpg(24) = 0.0187177221401866412d0
            xpg(25) = 0.3333333333333333333d0
            ypg(25) = 0.3333333333333333333d0
            hpg(25) = 0.0390780262448660605d0

        else if (fapg .eq. 'FPG28') then
! --------- Order 11 : Very high-order symmetric positive-interior
!                      quadrature rules on triangles and tetrahedra
!                      Zelalem Arega Worku  and Jason E. Hicken
!           https://github.com/OptimalDesignLab/SummationByParts.jl
            xpg(1) = 0.0306255245353121603d0
            ypg(1) = 0.9387489509293756793d0
            hpg(1) = 0.0060348721385088178d0
            xpg(2) = 0.0306255245353121603d0
            ypg(2) = 0.0306255245353121603d0
            hpg(2) = 0.0060348721385088178d0
            xpg(3) = 0.9387489509293756793d0
            ypg(3) = 0.0306255245353121603d0
            hpg(3) = 0.0060348721385088178d0
            xpg(4) = 0.1120794384974225991d0
            ypg(4) = 0.7758411230051548016d0
            hpg(4) = 0.0199610731868902660d0
            xpg(5) = 0.1120794384974225991d0
            ypg(5) = 0.1120794384974225991d0
            hpg(5) = 0.0199610731868902660d0
            xpg(6) = 0.7758411230051548016d0
            ypg(6) = 0.1120794384974225991d0
            hpg(6) = 0.0199610731868902660d0
            xpg(7) = 0.2139994793245382346d0
            ypg(7) = 0.5720010413509235306d0
            hpg(7) = 0.0340877415313096921d0
            xpg(8) = 0.2139994793245382346d0
            ypg(8) = 0.2139994793245382346d0
            hpg(8) = 0.0340877415313096921d0
            xpg(9) = 0.5720010413509235306d0
            ypg(9) = 0.2139994793245382346d0
            hpg(9) = 0.0340877415313096921d0
            xpg(10) = 0.4983342186162269636d0
            ypg(10) = 0.0033315627675460727d0
            hpg(10) = 0.0065100046990440848d0
            xpg(11) = 0.4983342186162269636d0
            ypg(11) = 0.4983342186162269636d0
            hpg(11) = 0.0065100046990440848d0
            xpg(12) = 0.0033315627675460727d0
            ypg(12) = 0.4983342186162269636d0
            hpg(12) = 0.0065100046990440848d0
            xpg(13) = 0.4369864234941642685d0
            ypg(13) = 0.1260271530116714628d0
            hpg(13) = 0.0318464179605561488d0
            xpg(14) = 0.4369864234941642685d0
            ypg(14) = 0.4369864234941642685d0
            hpg(14) = 0.0318464179605561488d0
            xpg(15) = 0.1260271530116714628d0
            ypg(15) = 0.4369864234941642685d0
            hpg(15) = 0.0318464179605561488d0
            xpg(16) = 0.0137073703144710521d0
            ypg(16) = 0.1582595826844803475d0
            hpg(16) = 0.0070380771122731376d0
            xpg(17) = 0.8280330470010486002d0
            ypg(17) = 0.1582595826844803475d0
            hpg(17) = 0.0070380771122731376d0
            xpg(18) = 0.8280330470010486002d0
            ypg(18) = 0.0137073703144710521d0
            hpg(18) = 0.0070380771122731376d0
            xpg(19) = 0.0137073703144710521d0
            ypg(19) = 0.8280330470010486002d0
            hpg(19) = 0.0070380771122731376d0
            xpg(20) = 0.1582595826844803475d0
            ypg(20) = 0.8280330470010486002d0
            hpg(20) = 0.0070380771122731376d0
            xpg(21) = 0.1582595826844803475d0
            ypg(21) = 0.0137073703144710521d0
            hpg(21) = 0.0070380771122731376d0
            xpg(22) = 0.0474766004458184237d0
            ypg(22) = 0.3079833895041687134d0
            hpg(22) = 0.0202382836943673032d0
            xpg(23) = 0.6445400100500128628d0
            ypg(23) = 0.3079833895041687134d0
            hpg(23) = 0.0202382836943673032d0
            xpg(24) = 0.6445400100500128628d0
            ypg(24) = 0.0474766004458184237d0
            hpg(24) = 0.0202382836943673032d0
            xpg(25) = 0.0474766004458184237d0
            ypg(25) = 0.6445400100500128628d0
            hpg(25) = 0.0202382836943673032d0
            xpg(26) = 0.3079833895041687134d0
            ypg(26) = 0.6445400100500128628d0
            hpg(26) = 0.0202382836943673032d0
            xpg(27) = 0.3079833895041687134d0
            ypg(27) = 0.0474766004458184237d0
            hpg(27) = 0.0202382836943673032d0
            xpg(28) = 0.3333333333333333333d0
            ypg(28) = 0.3333333333333333333d0
            hpg(28) = 0.0410215066112303255d0

        else if (fapg .eq. 'FPG33') then
! --------- Order 12 : Very high-order symmetric positive-interior
!                      quadrature rules on triangles and tetrahedra
!                      Zelalem Arega Worku  and Jason E. Hicken
!           https://github.com/OptimalDesignLab/SummationByParts.jl
            xpg(1) = 0.0213173504532103702d0
            ypg(1) = 0.9573652990935792594d0
            hpg(1) = 0.0030831305257795086d0
            xpg(2) = 0.0213173504532103702d0
            ypg(2) = 0.0213173504532103702d0
            hpg(2) = 0.0030831305257795086d0
            xpg(3) = 0.9573652990935792594d0
            ypg(3) = 0.0213173504532103702d0
            hpg(3) = 0.0030831305257795086d0
            xpg(4) = 0.1275761455415859246d0
            ypg(4) = 0.7448477089168281507d0
            hpg(4) = 0.0173980564653544714d0
            xpg(5) = 0.1275761455415859246d0
            ypg(5) = 0.1275761455415859246d0
            hpg(5) = 0.0173980564653544714d0
            xpg(6) = 0.7448477089168281507d0
            ypg(6) = 0.1275761455415859246d0
            hpg(6) = 0.0173980564653544714d0
            xpg(7) = 0.2712103850121159223d0
            ypg(7) = 0.4575792299757681553d0
            hpg(7) = 0.0314291121089425501d0
            xpg(8) = 0.2712103850121159223d0
            ypg(8) = 0.2712103850121159223d0
            hpg(8) = 0.0314291121089425501d0
            xpg(9) = 0.4575792299757681553d0
            ypg(9) = 0.2712103850121159223d0
            hpg(9) = 0.0314291121089425501d0
            xpg(10) = 0.4882173897738048825d0
            ypg(10) = 0.0235652204523902348d0
            hpg(10) = 0.0128655332202276677d0
            xpg(11) = 0.4882173897738048825d0
            ypg(11) = 0.4882173897738048825d0
            hpg(11) = 0.0128655332202276677d0
            xpg(12) = 0.0235652204523902348d0
            ypg(12) = 0.4882173897738048825d0
            hpg(12) = 0.0128655332202276677d0
            xpg(13) = 0.4397243922944602729d0
            ypg(13) = 0.1205512154110794540d0
            hpg(13) = 0.0218462722690192010d0
            xpg(14) = 0.4397243922944602729d0
            ypg(14) = 0.4397243922944602729d0
            hpg(14) = 0.0218462722690192010d0
            xpg(15) = 0.1205512154110794540d0
            ypg(15) = 0.4397243922944602729d0
            hpg(15) = 0.0218462722690192010d0
            xpg(16) = 0.0257340505483302281d0
            ypg(16) = 0.1162519159075971412d0
            hpg(16) = 0.0086581155543294461d0
            xpg(17) = 0.8580140335440726305d0
            ypg(17) = 0.1162519159075971412d0
            hpg(17) = 0.0086581155543294461d0
            xpg(18) = 0.8580140335440726305d0
            ypg(18) = 0.0257340505483302281d0
            hpg(18) = 0.0086581155543294461d0
            xpg(19) = 0.0257340505483302281d0
            ypg(19) = 0.8580140335440726305d0
            hpg(19) = 0.0086581155543294461d0
            xpg(20) = 0.1162519159075971412d0
            ypg(20) = 0.8580140335440726305d0
            hpg(20) = 0.0086581155543294461d0
            xpg(21) = 0.1162519159075971412d0
            ypg(21) = 0.0257340505483302281d0
            hpg(21) = 0.0086581155543294461d0
            xpg(22) = 0.0228383322222570296d0
            ypg(22) = 0.2813255809899395482d0
            hpg(22) = 0.0111783866011517228d0
            xpg(23) = 0.6958360867878034221d0
            ypg(23) = 0.2813255809899395482d0
            hpg(23) = 0.0111783866011517228d0
            xpg(24) = 0.6958360867878034221d0
            ypg(24) = 0.0228383322222570296d0
            hpg(24) = 0.0111783866011517228d0
            xpg(25) = 0.0228383322222570296d0
            ypg(25) = 0.6958360867878034221d0
            hpg(25) = 0.0111783866011517228d0
            xpg(26) = 0.2813255809899395482d0
            ypg(26) = 0.6958360867878034221d0
            hpg(26) = 0.0111783866011517228d0
            xpg(27) = 0.2813255809899395482d0
            ypg(27) = 0.0228383322222570296d0
            hpg(27) = 0.0111783866011517228d0
            xpg(28) = 0.1153434945346979991d0
            ypg(28) = 0.275713269685514194d0
            hpg(28) = 0.0201857788831904647d0
            xpg(29) = 0.6089432357797878068d0
            ypg(29) = 0.275713269685514194d0
            hpg(29) = 0.0201857788831904647d0
            xpg(30) = 0.6089432357797878068d0
            ypg(30) = 0.1153434945346979991d0
            hpg(30) = 0.0201857788831904647d0
            xpg(31) = 0.1153434945346979991d0
            ypg(31) = 0.6089432357797878068d0
            hpg(31) = 0.0201857788831904647d0
            xpg(32) = 0.275713269685514194d0
            ypg(32) = 0.6089432357797878068d0
            hpg(32) = 0.0201857788831904647d0
            xpg(33) = 0.275713269685514194d0
            ypg(33) = 0.1153434945346979991d0
            hpg(33) = 0.0201857788831904647d0

        else if (fapg .eq. 'FPG37') then
! --------- Order 13 : Very high-order symmetric positive-interior
!                      quadrature rules on triangles and tetrahedra
!                      Zelalem Arega Worku  and Jason E. Hicken
!           https://github.com/OptimalDesignLab/SummationByParts.jl
            xpg(1) = 0.0249000475950202539d0
            ypg(1) = 0.9501999048099594921d0
            hpg(1) = 0.0040160266401608808d0
            xpg(2) = 0.0249000475950202539d0
            ypg(2) = 0.0249000475950202539d0
            hpg(2) = 0.0040160266401608808d0
            xpg(3) = 0.9501999048099594921d0
            ypg(3) = 0.0249000475950202539d0
            hpg(3) = 0.0040160266401608808d0
            xpg(4) = 0.1135142530658452311d0
            ypg(4) = 0.7729714938683095378d0
            hpg(4) = 0.0154441132150080990d0
            xpg(5) = 0.1135142530658452311d0
            ypg(5) = 0.1135142530658452311d0
            hpg(5) = 0.0154441132150080990d0
            xpg(6) = 0.7729714938683095378d0
            ypg(6) = 0.1135142530658452311d0
            hpg(6) = 0.0154441132150080990d0
            xpg(7) = 0.2312256581889772392d0
            ypg(7) = 0.5375486836220455216d0
            hpg(7) = 0.0229987917394664872d0
            xpg(8) = 0.2312256581889772392d0
            ypg(8) = 0.2312256581889772392d0
            hpg(8) = 0.0229987917394664872d0
            xpg(9) = 0.5375486836220455216d0
            ypg(9) = 0.2312256581889772392d0
            hpg(9) = 0.0229987917394664872d0
            xpg(10) = 0.4962502629183093431d0
            ypg(10) = 0.0074994741633813136d0
            hpg(10) = 0.0049011586458473182d0
            xpg(11) = 0.4962502629183093431d0
            ypg(11) = 0.4962502629183093431d0
            hpg(11) = 0.0049011586458473182d0
            xpg(12) = 0.0074994741633813136d0
            ypg(12) = 0.4962502629183093431d0
            hpg(12) = 0.0049011586458473182d0
            xpg(13) = 0.4697418433578275545d0
            ypg(13) = 0.0605163132843448909d0
            hpg(13) = 0.0164291923372711773d0
            xpg(14) = 0.4697418433578275545d0
            ypg(14) = 0.4697418433578275545d0
            hpg(14) = 0.0164291923372711773d0
            xpg(15) = 0.0605163132843448909d0
            ypg(15) = 0.4697418433578275545d0
            hpg(15) = 0.0164291923372711773d0
            xpg(16) = 0.4145466153608832050d0
            ypg(16) = 0.1709067692782335899d0
            hpg(16) = 0.0234848668864528373d0
            xpg(17) = 0.4145466153608832050d0
            ypg(17) = 0.4145466153608832050d0
            hpg(17) = 0.0234848668864528373d0
            xpg(18) = 0.1709067692782335899d0
            ypg(18) = 0.4145466153608832050d0
            hpg(18) = 0.0234848668864528373d0
            xpg(19) = 0.0219525540418703654d0
            ypg(19) = 0.1268226967738466181d0
            hpg(19) = 0.0076932418768962641d0
            xpg(20) = 0.8512247491842830164d0
            ypg(20) = 0.1268226967738466181d0
            hpg(20) = 0.0076932418768962641d0
            xpg(21) = 0.8512247491842830164d0
            ypg(21) = 0.0219525540418703654d0
            hpg(21) = 0.0076932418768962641d0
            xpg(22) = 0.0219525540418703654d0
            ypg(22) = 0.8512247491842830164d0
            hpg(22) = 0.0076932418768962641d0
            xpg(23) = 0.1268226967738466181d0
            ypg(23) = 0.8512247491842830164d0
            hpg(23) = 0.0076932418768962641d0
            xpg(24) = 0.1268226967738466181d0
            ypg(24) = 0.0219525540418703654d0
            hpg(24) = 0.0076932418768962641d0
            xpg(25) = 0.0190342631016713422d0
            ypg(25) = 0.2920961426868683708d0
            hpg(25) = 0.0090819190938150692d0
            xpg(26) = 0.688869594211460287d0
            ypg(26) = 0.2920961426868683708d0
            hpg(26) = 0.0090819190938150692d0
            xpg(27) = 0.688869594211460287d0
            ypg(27) = 0.0190342631016713422d0
            hpg(27) = 0.0090819190938150692d0
            xpg(28) = 0.0190342631016713422d0
            ypg(28) = 0.688869594211460287d0
            hpg(28) = 0.0090819190938150692d0
            xpg(29) = 0.2920961426868683708d0
            ypg(29) = 0.688869594211460287d0
            hpg(29) = 0.0090819190938150692d0
            xpg(30) = 0.2920961426868683708d0
            ypg(30) = 0.0190342631016713422d0
            hpg(30) = 0.0090819190938150692d0
            xpg(31) = 0.0978946109941108393d0
            ypg(31) = 0.2666310672358204842d0
            hpg(31) = 0.0186181629969218637d0
            xpg(32) = 0.6354743217700686764d0
            ypg(32) = 0.2666310672358204842d0
            hpg(32) = 0.0186181629969218637d0
            xpg(33) = 0.6354743217700686764d0
            ypg(33) = 0.0978946109941108393d0
            hpg(33) = 0.0186181629969218637d0
            xpg(34) = 0.0978946109941108393d0
            ypg(34) = 0.6354743217700686764d0
            hpg(34) = 0.0186181629969218637d0
            xpg(35) = 0.2666310672358204842d0
            ypg(35) = 0.6354743217700686764d0
            hpg(35) = 0.0186181629969218637d0
            xpg(36) = 0.2666310672358204842d0
            ypg(36) = 0.0978946109941108393d0
            hpg(36) = 0.0186181629969218637d0
            xpg(37) = 0.3333333333333333333d0
            ypg(37) = 0.3333333333333333333d0
            hpg(37) = 0.0258176078015804164d0

        else if (fapg .eq. 'FPG42') then
! --------- Order 14 : Very high-order symmetric positive-interior
!                      quadrature rules on triangles and tetrahedra
!                      Zelalem Arega Worku  and Jason E. Hicken
!           https://github.com/OptimalDesignLab/SummationByParts.jl
            xpg(1) = 0.0193909612487010481d0
            ypg(1) = 0.9612180775025979037d0
            hpg(1) = 0.0024617018012000408d0
            xpg(2) = 0.0193909612487010481d0
            ypg(2) = 0.0193909612487010481d0
            hpg(2) = 0.0024617018012000408d0
            xpg(3) = 0.9612180775025979037d0
            ypg(3) = 0.0193909612487010481d0
            hpg(3) = 0.0024617018012000408d0
            xpg(4) = 0.0617998830908726012d0
            ypg(4) = 0.8764002338182547974d0
            hpg(4) = 0.0072168498348883338d0
            xpg(5) = 0.0617998830908726012d0
            ypg(5) = 0.0617998830908726012d0
            hpg(5) = 0.0072168498348883338d0
            xpg(6) = 0.8764002338182547974d0
            ypg(6) = 0.0617998830908726012d0
            hpg(6) = 0.0072168498348883338d0
            xpg(7) = 0.1772055324125434369d0
            ypg(7) = 0.6455889351749131261d0
            hpg(7) = 0.0210812943684965087d0
            xpg(8) = 0.1772055324125434369d0
            ypg(8) = 0.1772055324125434369d0
            hpg(8) = 0.0210812943684965087d0
            xpg(9) = 0.6455889351749131261d0
            ypg(9) = 0.1772055324125434369d0
            hpg(9) = 0.0210812943684965087d0
            xpg(10) = 0.2734775283088386597d0
            ypg(10) = 0.4530449433823226805d0
            hpg(10) = 0.0258870522536457931d0
            xpg(11) = 0.2734775283088386597d0
            ypg(11) = 0.2734775283088386597d0
            hpg(11) = 0.0258870522536457931d0
            xpg(12) = 0.4530449433823226805d0
            ypg(12) = 0.2734775283088386597d0
            hpg(12) = 0.0258870522536457931d0
            xpg(13) = 0.4889639103621786386d0
            ypg(13) = 0.0220721792756427226d0
            hpg(13) = 0.0109417906847144453d0
            xpg(14) = 0.4889639103621786386d0
            ypg(14) = 0.4889639103621786386d0
            hpg(14) = 0.0109417906847144453d0
            xpg(15) = 0.0220721792756427226d0
            ypg(15) = 0.4889639103621786386d0
            hpg(15) = 0.0109417906847144453d0
            xpg(16) = 0.4176447193404539225d0
            ypg(16) = 0.1647105613190921549d0
            hpg(16) = 0.0163941767720626753d0
            xpg(17) = 0.4176447193404539225d0
            ypg(17) = 0.4176447193404539225d0
            hpg(17) = 0.0163941767720626753d0
            xpg(18) = 0.1647105613190921549d0
            ypg(18) = 0.4176447193404539225d0
            hpg(18) = 0.0163941767720626753d0
            xpg(19) = 0.0012683309328720250d0
            ypg(19) = 0.1189744976969568453d0
            hpg(19) = 0.0025051144192503358d0
            xpg(20) = 0.8797571713701711295d0
            ypg(20) = 0.1189744976969568453d0
            hpg(20) = 0.0025051144192503358d0
            xpg(21) = 0.8797571713701711295d0
            ypg(21) = 0.0012683309328720250d0
            hpg(21) = 0.0025051144192503358d0
            xpg(22) = 0.0012683309328720250d0
            ypg(22) = 0.8797571713701711295d0
            hpg(22) = 0.0025051144192503358d0
            xpg(23) = 0.1189744976969568453d0
            ypg(23) = 0.8797571713701711295d0
            hpg(23) = 0.0025051144192503358d0
            xpg(24) = 0.1189744976969568453d0
            ypg(24) = 0.0012683309328720250d0
            hpg(24) = 0.0025051144192503358d0
            xpg(25) = 0.0146469500556544096d0
            ypg(25) = 0.298372882136257753d0
            hpg(25) = 0.0072181540567669202d0
            xpg(26) = 0.6869801678080878374d0
            ypg(26) = 0.298372882136257753d0
            hpg(26) = 0.0072181540567669202d0
            xpg(27) = 0.6869801678080878374d0
            ypg(27) = 0.0146469500556544096d0
            hpg(27) = 0.0072181540567669202d0
            xpg(28) = 0.0146469500556544096d0
            ypg(28) = 0.6869801678080878374d0
            hpg(28) = 0.0072181540567669202d0
            xpg(29) = 0.298372882136257753d0
            ypg(29) = 0.6869801678080878374d0
            hpg(29) = 0.0072181540567669202d0
            xpg(30) = 0.298372882136257753d0
            ypg(30) = 0.0146469500556544096d0
            hpg(30) = 0.0072181540567669202d0
            xpg(31) = 0.0571247574036479390d0
            ypg(31) = 0.1722666878213555783d0
            hpg(31) = 0.0123328766062818369d0
            xpg(32) = 0.7706085547749964826d0
            ypg(32) = 0.1722666878213555783d0
            hpg(32) = 0.0123328766062818369d0
            xpg(33) = 0.7706085547749964826d0
            ypg(33) = 0.0571247574036479390d0
            hpg(33) = 0.0123328766062818369d0
            xpg(34) = 0.0571247574036479390d0
            ypg(34) = 0.7706085547749964826d0
            hpg(34) = 0.0123328766062818369d0
            xpg(35) = 0.1722666878213555783d0
            ypg(35) = 0.7706085547749964826d0
            hpg(35) = 0.0123328766062818369d0
            xpg(36) = 0.1722666878213555783d0
            ypg(36) = 0.0571247574036479390d0
            hpg(36) = 0.0123328766062818369d0
            xpg(37) = 0.0929162493569718247d0
            ypg(37) = 0.3368614597963450017d0
            hpg(37) = 0.0192857553935303416d0
            xpg(38) = 0.5702222908466831735d0
            ypg(38) = 0.3368614597963450017d0
            hpg(38) = 0.0192857553935303416d0
            xpg(39) = 0.5702222908466831735d0
            ypg(39) = 0.0929162493569718247d0
            hpg(39) = 0.0192857553935303416d0
            xpg(40) = 0.0929162493569718247d0
            ypg(40) = 0.5702222908466831735d0
            hpg(40) = 0.0192857553935303416d0
            xpg(41) = 0.3368614597963450017d0
            ypg(41) = 0.5702222908466831735d0
            hpg(41) = 0.0192857553935303416d0
            xpg(42) = 0.3368614597963450017d0
            ypg(42) = 0.0929162493569718247d0
            hpg(42) = 0.0192857553935303416d0

        else if (fapg .eq. 'COT3') then
            xpg(1) = undemi
            ypg(1) = undemi
            xpg(2) = zero
            ypg(2) = undemi
            xpg(3) = undemi
            ypg(3) = zero
            hpg(1) = un/6.d0
            hpg(2) = un/6.d0
            hpg(3) = un/6.d0

        else if (fapg .eq. 'SIMP') then
            xpg(1) = zero
            ypg(1) = zero
            xpg(2) = un
            ypg(2) = zero
            xpg(3) = zero
            ypg(3) = un
            xpg(4) = undemi
            ypg(4) = zero
            xpg(5) = undemi
            ypg(5) = undemi
            xpg(6) = zero
            ypg(6) = undemi
            hpg(1) = un/30.d0
            hpg(2) = un/30.d0
            hpg(3) = un/30.d0
            hpg(4) = 4.d0/30.d0
            hpg(5) = 4.d0/30.d0
            hpg(6) = 4.d0/30.d0

        else if (fapg .eq. 'FPG3NOS') then
! --------- POUR LES POINTS DE GAUSS
            xpg(1) = un/6.d0
            ypg(1) = un/6.d0
            xpg(2) = deux/3.d0
            ypg(2) = un/6.d0
            xpg(3) = un/6.d0
            ypg(3) = deux/3.d0
            hpg(1) = un/6.d0
            hpg(2) = un/6.d0
            hpg(3) = un/6.d0
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+3) = vol/nnos
                xpg(iNode+3) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+3) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+3) = xno(ndim*(iNode-1)+3)
            end do

        else
            ASSERT(ASTER_FALSE)

        end if

    else if (elrefa .eq. 'QU4' .or. elrefa .eq. 'QU8' .or. elrefa .eq. 'QU9' &
             .or. elrefa .eq. 'Q12') then
        if (fapg .eq. 'FPG1') then
            xpg(1) = zero
            ypg(1) = zero
            hpg(1) = 4.d0

        else if (fapg .eq. 'FIS2') then
! ------- ELEMENT PARTICULIER DE FISSURE, S'APPUIE SUR UN SEG2
            xpg(1) = -rac_1div3
            ypg(1) = zero
            xpg(2) = rac_1div3
            ypg(2) = zero
            hpg(1) = deux
            hpg(2) = deux

        else if (fapg .eq. 'FPG4') then
            xpg(1) = -rac_1div3
            ypg(1) = -rac_1div3
            xpg(2) = rac_1div3
            ypg(2) = -rac_1div3
            xpg(3) = rac_1div3
            ypg(3) = rac_1div3
            xpg(4) = -rac_1div3
            ypg(4) = rac_1div3
            hpg(1) = un
            hpg(2) = un
            hpg(3) = un
            hpg(4) = un

        else if (fapg .eq. 'FPG9') then
            hpg(1) = 25.d0/81.0d0
            hpg(2) = 25.d0/81.0d0
            hpg(3) = 25.d0/81.0d0
            hpg(4) = 25.d0/81.0d0
            hpg(5) = 40.d0/81.0d0
            hpg(6) = 40.d0/81.0d0
            hpg(7) = 40.d0/81.0d0
            hpg(8) = 40.d0/81.0d0
            hpg(9) = 64.d0/81.0d0
            xpg(1) = -rac_3div5
            ypg(1) = -rac_3div5
            xpg(2) = rac_3div5
            ypg(2) = -rac_3div5
            xpg(3) = rac_3div5
            ypg(3) = rac_3div5
            xpg(4) = -rac_3div5
            ypg(4) = rac_3div5
            xpg(5) = zero
            ypg(5) = -rac_3div5
            xpg(6) = rac_3div5
            ypg(6) = zero
            xpg(7) = zero
            ypg(7) = rac_3div5
            xpg(8) = -rac_3div5
            ypg(8) = zero
            xpg(9) = zero
            ypg(9) = zero

        else if (fapg .eq. 'FPG9COQ') then
            hpg(7) = 25.d0/81.0d0
            hpg(1) = 25.d0/81.0d0
            hpg(3) = 25.d0/81.0d0
            hpg(5) = 25.d0/81.0d0
            hpg(8) = 40.d0/81.0d0
            hpg(2) = 40.d0/81.0d0
            hpg(4) = 40.d0/81.0d0
            hpg(6) = 40.d0/81.0d0
            hpg(9) = 64.d0/81.0d0
            xpg(1) = -rac_3div5
            ypg(1) = -rac_3div5
            xpg(3) = rac_3div5
            ypg(3) = -rac_3div5
            xpg(5) = rac_3div5
            ypg(5) = rac_3div5
            xpg(7) = -rac_3div5
            ypg(7) = rac_3div5
            xpg(2) = zero
            ypg(2) = -rac_3div5
            xpg(4) = rac_3div5
            ypg(4) = zero
            xpg(6) = zero
            ypg(6) = rac_3div5
            xpg(8) = -rac_3div5
            ypg(8) = zero
            xpg(9) = zero
            ypg(9) = zero

        else if (fapg .eq. 'FPG16') then
            h(1) = (18.d0+rac30)/36.d0
            h(2) = h(1)
            h(3) = (18.d0-rac30)/36.d0
            h(4) = h(3)
            a(1) = -gauss4p12
            a(2) = -a(1)
            a(3) = -gauss4p34
            a(4) = -a(3)
            npar = 4
            npi = 0
            do ix = 1, npar
                do iy = 1, npar
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = a(iy)
                    hpg(npi) = h(ix)*h(iy)
                end do
            end do

        else if (fapg .eq. 'FPG25') then
            ! order 9 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            a(1) = 0.d0
            a(2) = 0.53846931010568309103d0
            a(3) = -a(2)
            a(4) = 0.90617984593866399279d0
            a(5) = -a(4)
            h(1) = 0.56888888888888888888d0
            h(2) = 0.47862867049936646804d0
            h(3) = h(2)
            h(4) = 0.23692688505618908751d0
            h(5) = h(4)
            npar = 5
            npi = 0
            do ix = 1, npar
                do iy = 1, npar
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = a(iy)
                    hpg(npi) = h(ix)*h(iy)
                end do
            end do

        else if (fapg .eq. 'FPG36') then
            ! order 11 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            a(1) = 0.661209386466264513661d0
            a(2) = -a(1)
            a(3) = 0.238619186083196908630d0
            a(4) = -a(3)
            a(5) = 0.932469514203152027812d0
            a(6) = -a(5)
            h(1) = 0.360761573048138607569d0
            h(2) = h(1)
            h(3) = 0.467913934572691047389d0
            h(4) = h(3)
            h(5) = 0.171324492379170345040d0
            h(6) = h(5)
            npar = 6
            npi = 0
            do ix = 1, npar
                do iy = 1, npar
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = a(iy)
                    hpg(npi) = h(ix)*h(iy)
                end do
            end do

        else if (fapg .eq. 'FPG49') then
            ! order 13 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            a(1) = 0.d0
            a(2) = 0.405845151377397166906d0
            a(3) = -a(2)
            a(4) = 0.741531185599394439863d0
            a(5) = -a(4)
            a(6) = 0.949107912342758524526d0
            a(7) = -a(6)
            h(1) = 0.417959183673469387755d0
            h(2) = 0.381830050505118944950d0
            h(3) = h(2)
            h(4) = 0.279705391489276667901d0
            h(5) = h(4)
            h(6) = 0.129484966168869693270d0
            h(7) = h(6)
            npar = 7
            npi = 0
            do ix = 1, npar
                do iy = 1, npar
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = a(iy)
                    hpg(npi) = h(ix)*h(iy)
                end do
            end do

        else if (fapg .eq. 'FPG64') then
            ! order 15 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            a(1) = 0.183434642495649804939d0
            a(2) = -a(1)
            a(3) = 0.525532409916328985817d0
            a(4) = -a(3)
            a(5) = 0.796666477413626739591d0
            a(6) = -a(5)
            a(7) = 0.960289856497536231683d0
            a(8) = -a(7)
            h(1) = 0.362683783378361982965d0
            h(2) = h(1)
            h(3) = 0.313706645877887287337d0
            h(4) = h(3)
            h(5) = 0.222381034453374470544d0
            h(6) = h(5)
            h(7) = 0.101228536290376259152d0
            h(8) = h(7)
            npar = 8
            npi = 0
            do ix = 1, npar
                do iy = 1, npar
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = a(iy)
                    hpg(npi) = h(ix)*h(iy)
                end do
            end do

        else if (fapg .eq. 'FPG4NOS') then
! --------- POUR LES POINTS DE GAUSS
            xpg(1) = -rac_1div3
            ypg(1) = -rac_1div3
            xpg(2) = rac_1div3
            ypg(2) = -rac_1div3
            xpg(3) = rac_1div3
            ypg(3) = rac_1div3
            xpg(4) = -rac_1div3
            ypg(4) = rac_1div3
            hpg(1) = un
            hpg(2) = un
            hpg(3) = un
            hpg(4) = un
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+4) = vol/nnos
                xpg(iNode+4) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+4) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+4) = xno(ndim*(iNode-1)+3)
            end do

        else
            ASSERT(ASTER_FALSE)

        end if

    else if (elrefa .eq. 'SE2' .or. elrefa .eq. 'SE3' .or. elrefa .eq. 'SE4') then
        if (fapg .eq. 'FPG1') then
            ! order 1
            xpg(1) = zero
            hpg(1) = deux

        else if (fapg .eq. 'FPG2') then
            ! order 3
            xpg(1) = rac_1div3
            xpg(2) = -xpg(1)
            hpg(1) = un
            hpg(2) = hpg(1)

        else if (fapg .eq. 'FPG3') then
            ! order 5
            xpg(1) = -rac_3div5
            xpg(2) = zero
            xpg(3) = rac_3div5
            hpg(1) = 5.d0/9.d0
            hpg(2) = 8.d0/9.d0
            hpg(3) = 5.d0/9.d0

        else if (fapg .eq. 'FPG4') then
            ! order 7
            xpg(1) = gauss4p12
            xpg(2) = -xpg(1)
            xpg(3) = gauss4p34
            xpg(4) = -xpg(3)
            hpg(1) = (18.d0+rac30)/36.d0
            hpg(2) = hpg(1)
            hpg(3) = (18.d0-rac30)/36.d0
            hpg(4) = hpg(3)

        else if (fapg .eq. 'FPG5') then
            ! order 9 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            xpg(1) = 0.d0
            xpg(2) = 0.53846931010568309103d0
            xpg(3) = -xpg(2)
            xpg(4) = 0.90617984593866399279d0
            xpg(5) = -xpg(4)
            hpg(1) = 0.56888888888888888888d0
            hpg(2) = 0.47862867049936646804d0
            hpg(3) = hpg(2)
            hpg(4) = 0.23692688505618908751d0
            hpg(5) = hpg(4)

        else if (fapg .eq. 'FPG6') then
            ! order 11 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            xpg(1) = 0.661209386466264513661d0
            xpg(2) = -xpg(1)
            xpg(3) = 0.238619186083196908630d0
            xpg(4) = -xpg(3)
            xpg(5) = 0.932469514203152027812d0
            xpg(6) = -xpg(5)
            hpg(1) = 0.360761573048138607569d0
            hpg(2) = hpg(1)
            hpg(3) = 0.467913934572691047389d0
            hpg(4) = hpg(3)
            hpg(5) = 0.171324492379170345040d0
            hpg(6) = hpg(5)

        else if (fapg .eq. 'FPG7') then
            ! order 13 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            xpg(1) = 0.d0
            xpg(2) = 0.405845151377397166906d0
            xpg(3) = -xpg(2)
            xpg(4) = 0.741531185599394439863d0
            xpg(5) = -xpg(4)
            xpg(6) = 0.949107912342758524526d0
            xpg(7) = -xpg(6)
            hpg(1) = 0.417959183673469387755d0
            hpg(2) = 0.381830050505118944950d0
            hpg(3) = hpg(2)
            hpg(4) = 0.279705391489276667901d0
            hpg(5) = hpg(4)
            hpg(6) = 0.129484966168869693270d0
            hpg(7) = hpg(6)

        else if (fapg .eq. 'FPG8') then
            ! order 15 : https://pomax.github.io/bezierinfo/legendre-gauss.html
            xpg(1) = 0.183434642495649804939d0
            xpg(2) = -xpg(1)
            xpg(3) = 0.525532409916328985817d0
            xpg(4) = -xpg(3)
            xpg(5) = 0.796666477413626739591d0
            xpg(6) = -xpg(5)
            xpg(7) = 0.960289856497536231683d0
            xpg(8) = -xpg(7)
            hpg(1) = 0.362683783378361982965d0
            hpg(2) = hpg(1)
            hpg(3) = 0.313706645877887287337d0
            hpg(4) = hpg(3)
            hpg(5) = 0.222381034453374470544d0
            hpg(6) = hpg(5)
            hpg(7) = 0.101228536290376259152d0
            hpg(8) = hpg(7)

        else if (fapg .eq. 'FPG2NOS') then
            xpg(1) = rac_1div3
            xpg(2) = -xpg(1)
            xpg(3) = xno(1)
            xpg(4) = xno(2)
            hpg(1) = un
            hpg(2) = hpg(1)
            hpg(3) = vol/nnos
            hpg(4) = hpg(3)

        else if (fapg .eq. 'FPG3NOS') then
            xpg(1) = -rac_3div5
            xpg(2) = zero
            xpg(3) = rac_3div5
            xpg(4) = xno(1)
            xpg(5) = xno(nnos)
            hpg(1) = 5.d0/9.d0
            hpg(2) = 8.d0/9.d0
            hpg(3) = 5.d0/9.d0
            hpg(4) = vol/nnos
            hpg(5) = hpg(4)

        else if (fapg .eq. 'SIMP') then
            xpg(1) = -un
            xpg(2) = zero
            xpg(3) = un
            hpg(1) = un/3.d0
            hpg(2) = 4.d0/3.d0
            hpg(3) = un/3.d0

        else if (fapg .eq. 'SIMP1') then
            xpg(1) = -un
            xpg(2) = -undemi
            xpg(3) = zero
            xpg(4) = undemi
            xpg(5) = un
            hpg(1) = un/6.d0
            hpg(2) = deux/3.d0
            hpg(3) = un/3.d0
            hpg(4) = deux/3.d0
            hpg(5) = un/6.d0

        else if (fapg .eq. 'COTES') then
            xpg(1) = -un
            xpg(2) = -un/3.d0
            xpg(3) = un/3.d0
            xpg(4) = un
            hpg(1) = un/4.d0
            hpg(2) = 3.d0/4.d0
            hpg(3) = 3.d0/4.d0
            hpg(4) = un/4.d0

        else if (fapg .eq. 'COTES1') then
            xpg(1) = -un
            xpg(2) = -un/deux
            xpg(3) = zero
            xpg(4) = un/deux
            xpg(5) = un
            hpg(1) = 7.d0/45.d0
            hpg(2) = 32.d0/45.d0
            hpg(3) = 12.d0/45.d0
            hpg(4) = 32.d0/45.d0
            hpg(5) = 7.d0/45.d0

        else if (fapg .eq. 'COTES2') then
            xpg(1) = -un
            xpg(2) = -7.d0/9.d0
            xpg(3) = -5.d0/9.d0
            xpg(4) = -un/3.d0
            xpg(5) = -un/9.d0
            xpg(6) = un/9.d0
            xpg(7) = un/3.d0
            xpg(8) = 5.d0/9.d0
            xpg(9) = 7.d0/9.d0
            xpg(10) = un
            hpg(1) = un/12.d0
            hpg(2) = un/4.d0
            hpg(3) = un/4.d0
            hpg(4) = un/6.d0
            hpg(5) = un/4.d0
            hpg(6) = un/4.d0
            hpg(7) = un/6.d0
            hpg(8) = un/4.d0
            hpg(9) = un/4.d0
            hpg(10) = un/12.d0

        else
            ASSERT(ASTER_FALSE)

        end if

    else if (elrefa .eq. 'PO1') then
        hpg(1) = un

    else
        ASSERT(ASTER_FALSE)
    end if
!
170 continue
!
    do i = 1, nbpg
        poipg(i) = hpg(i)
        if (ndim .ge. 1) coopg(ndim*(i-1)+1) = xpg(i)
        if (ndim .ge. 2) coopg(ndim*(i-1)+2) = ypg(i)
        if (ndim .eq. 3) coopg(ndim*(i-1)+3) = zpg(i)
    end do

end subroutine
