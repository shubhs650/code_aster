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
subroutine elraca(elrefz, &
                  nbfpg_, fapg_, nbpg_, &
                  ndim_, nno_, nnos_, &
                  nodeCoor_, cellVolu_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrfno.h"
#include "MeshTypes_type.h"
!
    character(len=*), intent(in)            :: elrefz
    integer(kind=8), optional, intent(out)  :: ndim_, nno_, nnos_, nbfpg_, nbpg_(MT_NBFAMX)
    real(kind=8), optional, intent(out)     :: nodeCoor_(3*MT_NNOMAX), cellVolu_
    character(len=8), optional, intent(out) :: fapg_(MT_NBFAMX)
!
! --------------------------------------------------------------------------------------------------
!
! Finite elements management
!
! Get list of integration schemes of geometric support
!
! --------------------------------------------------------------------------------------------------
!
! In  elrefa           : name of geometric support
! Out ndim             : topological dimension (0/1/2/3)
! Out nno              : number of nodes
! Out nnos             : number of middle nodes
! Out nbfpg            : number of families of integration schemes
! Out fapg             : name of families for all integration schemes
! Out nbpg             : number of points for all integration schemes
! Out nodeCoor         : coordinates of node of geometric support in parametric space
! Out cellVolu         : volume of geometric support in parametric space
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elrefa
    integer(kind=8) :: i, deca
    integer(kind=8) :: ndim, nno, nnos, nbfpg, nbpg(MT_NBFAMX)
    character(len=8) :: fapg(MT_NBFAMX)
    real(kind=8) :: coorno(3, MT_NNOMAX), nodeCoor(3*MT_NNOMAX), cellVolu
!
! --------------------------------------------------------------------------------------------------
!
    elrefa = elrefz
    ndim = 0
    nno = 0
    nnos = 0
    nbfpg = 0
    nbpg = 0
    nodeCoor = 0.d0
    fapg = ' '
    cellVolu = 0.d0

! - Get main parameters of geometric support in parametric space
    call elrfno(elrefa, nno, nnos, ndim, coorno, cellVolu)

    select case (elrefa)
    case ('HE8')
        nbfpg = 11
        nbpg(1:nbfpg) = [nno, nnos, 1, 8, 27, 16, 64, 125, 216, 343, 512]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG8'
        fapg(5) = 'FPG27'
        fapg(6) = 'FPG8NOS'
        fapg(7) = 'FPG64'
        fapg(8) = 'FPG125'
        fapg(9) = 'FPG216'
        fapg(10) = 'FPG343'
        fapg(11) = 'FPG512'

    case ('HE9')
        nbfpg = 6
        nbpg(1:nbfpg) = [nno, nnos, 1, 5, 7, 8]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'LOB5'
        fapg(5) = 'LOB7'
        fapg(6) = 'FPG8'

    case ('H10')
        nbfpg = 3
        nbpg(1:nbfpg) = [3, 4, 2]
        fapg(1) = 'FIS3'
        fapg(2) = 'FIS4'
        fapg(3) = 'FIS2'

    case ('H20')
        nbfpg = 7
        nbpg(1:nbfpg) = [nno, nnos, 1, 8, 27, 16, 64]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG8'
        fapg(5) = 'FPG27'
        fapg(6) = 'FPG8NOS'
        fapg(7) = 'FPG64'

    case ('H27')
        nbfpg = 10
        nbpg(1:nbfpg) = [nno, nnos, 1, 8, 27, 64, 125, 216, 343, 512]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG8'
        fapg(5) = 'FPG27'
        fapg(6) = 'FPG64'
        fapg(7) = 'FPG125'
        fapg(8) = 'FPG216'
        fapg(9) = 'FPG343'
        fapg(10) = 'FPG512'

    case ('TE4')
        nbfpg = 17
        nbpg(1:nbfpg) = [nno, nnos, 1, 4, 5, 11, 15, 24, 8, &
                         35, 46, 59, 74, 94, 117, 144, 204]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG4'
        fapg(5) = 'FPG5'
        fapg(6) = 'FPG11'
        fapg(7) = 'FPG15'
        fapg(8) = 'FPG24'
        fapg(9) = 'FPG4NOS'
        fapg(10) = 'FPG35'
        fapg(11) = 'FPG46'
        fapg(12) = 'FPG59'
        fapg(13) = 'FPG74'
        fapg(14) = 'FPG94'
        fapg(15) = 'FPG117'
        fapg(16) = 'FPG144'
        fapg(17) = 'FPG204'

    case ('T10')
        nbfpg = 7
        nbpg(1:nbfpg) = [nno, nnos, 1, 4, 5, 15, 8]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG4'
        fapg(5) = 'FPG5'
        fapg(6) = 'FPG15'
        fapg(7) = 'FPG4NOS'

    case ('T15')
        nbfpg = 17
        nbpg(1:nbfpg) = [nno, nnos, 1, 4, 5, 11, 15, 24, 8, &
                         35, 46, 59, 74, 94, 117, 144, 204]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG4'
        fapg(5) = 'FPG5'
        fapg(6) = 'FPG11'
        fapg(7) = 'FPG15'
        fapg(8) = 'FPG24'
        fapg(9) = 'FPG4NOS'
        fapg(10) = 'FPG35'
        fapg(11) = 'FPG46'
        fapg(12) = 'FPG59'
        fapg(13) = 'FPG74'
        fapg(14) = 'FPG94'
        fapg(15) = 'FPG117'
        fapg(16) = 'FPG144'
        fapg(17) = 'FPG204'

    case ('T20')
        nbfpg = 17
        nbpg(1:nbfpg) = [nno, nnos, 1, 4, 5, 11, 15, 24, 8, &
                         35, 46, 59, 74, 94, 117, 144, 204]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG4'
        fapg(5) = 'FPG5'
        fapg(6) = 'FPG11'
        fapg(7) = 'FPG15'
        fapg(8) = 'FPG24'
        fapg(9) = 'FPG4NOS'
        fapg(10) = 'FPG35'
        fapg(11) = 'FPG46'
        fapg(12) = 'FPG59'
        fapg(13) = 'FPG74'
        fapg(14) = 'FPG94'
        fapg(15) = 'FPG117'
        fapg(16) = 'FPG144'
        fapg(17) = 'FPG204'

    case ('PE6')
        nbfpg = 12
        nbpg(1:nbfpg) = [nno, nnos, 1, 6, 6, 8, 21, 12, 29, 52, 95, 168]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG6'
        fapg(5) = 'FPG6B'
        fapg(6) = 'FPG8'
        fapg(7) = 'FPG21'
        fapg(8) = 'FPG6NOS'
        fapg(9) = 'FPG29'
        fapg(10) = 'FPG52'
        fapg(11) = 'FPG95'
        fapg(12) = 'FPG168'

    case ('PE7')
        nbfpg = 5
        nbpg(1:nbfpg) = [nno, nnos, 1, 5, 7]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'LOB5'
        fapg(5) = 'LOB7'

    case ('P15')
        nbfpg = 9
        nbpg(1:nbfpg) = [nno, nnos, 1, 6, 6, 8, 21, 12, 29]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG6'
        fapg(5) = 'FPG6B'
        fapg(6) = 'FPG8'
        fapg(7) = 'FPG21'
        fapg(8) = 'FPG6NOS'
        fapg(9) = 'FPG29'

    case ('P18')
        nbfpg = 9
        nbpg(1:nbfpg) = [nno, nnos, 1, 6, 6, 8, 21, 12, 29]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG6'
        fapg(5) = 'FPG6B'
        fapg(6) = 'FPG8'
        fapg(7) = 'FPG21'
        fapg(8) = 'FPG6NOS'
        fapg(9) = 'FPG29'

    case ('P21')
        nbfpg = 12
        nbpg(1:nbfpg) = [nno, nnos, 1, 6, 6, 8, 21, 12, 29, 52, 95, 168]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG6'
        fapg(5) = 'FPG6B'
        fapg(6) = 'FPG8'
        fapg(7) = 'FPG21'
        fapg(8) = 'FPG6NOS'
        fapg(9) = 'FPG29'
        fapg(10) = 'FPG52'
        fapg(11) = 'FPG95'
        fapg(12) = 'FPG168'

    case ('PY5')
        nbfpg = 14
        nbpg(1:nbfpg) = [nno, nnos, 1, 1, 5, 6, 10, 10, 15, 24, 31, 47, 62, 83]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG1B'
        fapg(5) = 'FPG5'
        fapg(6) = 'FPG6'
        fapg(7) = 'FPG10'
        fapg(8) = 'FPG5NOS'
        fapg(9) = 'FPG15'
        fapg(10) = 'FPG24'
        fapg(11) = 'FPG31'
        fapg(12) = 'FPG47'
        fapg(13) = 'FPG62'
        fapg(14) = 'FPG83'

    case ('P13')
        nbfpg = 11
        nbpg(1:nbfpg) = [nno, nnos, 1, 1, 5, 6, 10, 10, 15, 24, 31]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG1B'
        fapg(5) = 'FPG5'
        fapg(6) = 'FPG6'
        fapg(7) = 'FPG10'
        fapg(8) = 'FPG5NOS'
        fapg(9) = 'FPG15'
        fapg(10) = 'FPG24'
        fapg(11) = 'FPG31'

    case ('P19')
        nbfpg = 14
        nbpg(1:nbfpg) = [nno, nnos, 1, 1, 5, 6, 10, 10, 15, 24, 31, 47, 62, 83]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG1B'
        fapg(5) = 'FPG5'
        fapg(6) = 'FPG6'
        fapg(7) = 'FPG10'
        fapg(8) = 'FPG5NOS'
        fapg(9) = 'FPG15'
        fapg(10) = 'FPG24'
        fapg(11) = 'FPG31'
        fapg(12) = 'FPG47'
        fapg(13) = 'FPG62'
        fapg(14) = 'FPG83'

    case ('TR3')
        nbfpg = 19
        nbpg(1:nbfpg) = [nno, nnos, 1, 3, 4, 6, 7, 12, 3, 6, 13, 16, 6, &
                         19, 25, 28, 33, 37, 42]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG3'
        fapg(5) = 'FPG4'
        fapg(6) = 'FPG6'
        fapg(7) = 'FPG7'
        fapg(8) = 'FPG12'
        fapg(9) = 'COT3'
        fapg(10) = 'FPG3NOS'
        fapg(11) = 'FPG13'
        fapg(12) = 'FPG16'
        fapg(13) = 'SIMP'
        fapg(14) = 'FPG19'
        fapg(15) = 'FPG25'
        fapg(16) = 'FPG28'
        fapg(17) = 'FPG33'
        fapg(18) = 'FPG37'
        fapg(19) = 'FPG42'

    case ('TR6')
        nbfpg = 11
        nbpg(1:nbfpg) = [nno, nnos, 1, 3, 4, 6, 7, 12, 6, 13, 16]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG3'
        fapg(5) = 'FPG4'
        fapg(6) = 'FPG6'
        fapg(7) = 'FPG7'
        fapg(8) = 'FPG12'
        fapg(9) = 'FPG3NOS'
        fapg(10) = 'FPG13'
        fapg(11) = 'FPG16'

    case ('TR7')
        nbfpg = 16
        nbpg(1:nbfpg) = [nno, nnos, 1, 3, 4, 6, 7, 12, 13, 16, &
                         19, 25, 28, 33, 37, 42]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG3'
        fapg(5) = 'FPG4'
        fapg(6) = 'FPG6'
        fapg(7) = 'FPG7'
        fapg(8) = 'FPG12'
        fapg(9) = 'FPG13'
        fapg(10) = 'FPG16'
        fapg(11) = 'FPG19'
        fapg(12) = 'FPG25'
        fapg(13) = 'FPG28'
        fapg(14) = 'FPG33'
        fapg(15) = 'FPG37'
        fapg(16) = 'FPG42'

    case ('TR1')
        nbfpg = 16
        nbpg(1:nbfpg) = [nno, nnos, 1, 3, 4, 6, 7, 12, 13, 16, &
                         19, 25, 28, 33, 37, 42]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG3'
        fapg(5) = 'FPG4'
        fapg(6) = 'FPG6'
        fapg(7) = 'FPG7'
        fapg(8) = 'FPG12'
        fapg(9) = 'FPG13'
        fapg(10) = 'FPG16'
        fapg(11) = 'FPG19'
        fapg(12) = 'FPG25'
        fapg(13) = 'FPG28'
        fapg(14) = 'FPG33'
        fapg(15) = 'FPG37'
        fapg(16) = 'FPG42'

    case ('QU4')
        nbfpg = 12
        nbpg(1:nbfpg) = [nno, nnos, 1, 4, 9, 16, 25, 36, 49, 64, &
                         2, 8]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG4'
        fapg(5) = 'FPG9'
        fapg(6) = 'FPG16'
        fapg(7) = 'FPG25'
        fapg(8) = 'FPG36'
        fapg(9) = 'FPG49'
        fapg(10) = 'FPG64'
        fapg(11) = 'FIS2'
        fapg(12) = 'FPG4NOS'

    case ('QU8')
        nbfpg = 8
        nbpg(1:nbfpg) = [nno, nnos, 1, 4, 9, 9, 8, 16]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG4'
        fapg(5) = 'FPG9'
        fapg(6) = 'FPG9COQ'
        fapg(7) = 'FPG4NOS'
        fapg(8) = 'FPG16'

    case ('QU9')
        nbfpg = 11
        nbpg(1:nbfpg) = [nno, nnos, 1, 4, 9, 9, 16, 25, 36, 49, 64]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG4'
        fapg(5) = 'FPG9'
        fapg(6) = 'FPG9COQ'
        fapg(7) = 'FPG16'
        fapg(8) = 'FPG25'
        fapg(9) = 'FPG36'
        fapg(10) = 'FPG49'
        fapg(11) = 'FPG64'

    case ('Q12')
        nbfpg = 11
        nbpg(1:nbfpg) = [nno, nnos, 1, 4, 9, 9, 16, 25, 36, 49, 64]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG4'
        fapg(5) = 'FPG9'
        fapg(6) = 'FPG9COQ'
        fapg(7) = 'FPG16'
        fapg(8) = 'FPG25'
        fapg(9) = 'FPG36'
        fapg(10) = 'FPG49'
        fapg(11) = 'FPG64'

    case ('SE2')
        nbfpg = 17
        nbpg(1:nbfpg) = [nno, nnos, 1, 2, 3, 4, 5, 6, 7, 8, &
                         3, 5, 4, 5, 10, nnos+2, nnos+3]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG2'
        fapg(5) = 'FPG3'
        fapg(6) = 'FPG4'
        fapg(7) = 'FPG5'
        fapg(8) = 'FPG6'
        fapg(9) = 'FPG7'
        fapg(10) = 'FPG8'
        fapg(11) = 'SIMP'
        fapg(12) = 'SIMP1'
        fapg(13) = 'COTES'
        fapg(14) = 'COTES1'
        fapg(15) = 'COTES2'
        fapg(16) = 'FPG2NOS'
        fapg(17) = 'FPG3NOS'

    case ('SE3')
        nbfpg = 14
        nbpg(1:nbfpg) = [nno, nnos, 1, 2, 3, 4, 3, 4, nnos+2, nnos+3, &
                         5, 6, 7, 8]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG2'
        fapg(5) = 'FPG3'
        fapg(6) = 'FPG4'
        fapg(7) = 'SIMP'
        fapg(8) = 'COTES'
        fapg(9) = 'FPG2NOS'
        fapg(10) = 'FPG3NOS'
        fapg(11) = 'FPG5'
        fapg(12) = 'FPG6'
        fapg(13) = 'FPG7'
        fapg(14) = 'FPG8'

    case ('SE4')
        nbfpg = 10
        nbpg(1:nbfpg) = [nno, nnos, 1, 2, 3, 4, 5, 6, 3, 4]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'
        fapg(4) = 'FPG2'
        fapg(5) = 'FPG3'
        fapg(6) = 'FPG4'
        fapg(7) = 'FPG5'
        fapg(8) = 'FPG6'
        fapg(9) = 'SIMP'
        fapg(10) = 'COTES'

    case ('PO1')
        nbfpg = 3
        nbpg(1:nbfpg) = [1, 1, 1]
        fapg(1) = 'NOEU'
        fapg(2) = 'NOEU_S'
        fapg(3) = 'FPG1'

    case default
        ASSERT(ASTER_FALSE)

    end select
!
    if (present(nodeCoor_)) then
        do i = 1, nno
            deca = ndim*(i-1)
            if (ndim .ge. 1) nodeCoor_(deca+1) = coorno(1, i)
            if (ndim .ge. 2) nodeCoor_(deca+2) = coorno(2, i)
            if (ndim .eq. 3) nodeCoor_(deca+3) = coorno(3, i)
        end do
    end if
    if (present(ndim_)) then
        ndim_ = ndim
    end if
    if (present(nno_)) then
        nno_ = nno
    end if
    if (present(nnos_)) then
        nnos_ = nnos
    end if
    if (present(nbfpg_)) then
        nbfpg_ = nbfpg
    end if
    if (present(nbpg_)) then
        nbpg_ = nbpg
    end if
    if (present(cellVolu_)) then
        cellVolu_ = cellVolu
    end if
    if (present(fapg_)) then
        fapg_ = fapg
    end if
!
end subroutine
