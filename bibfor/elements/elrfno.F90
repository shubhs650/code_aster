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
subroutine elrfno(elrefz, nno, nnos, ndim, nodeCoor, cellVolu)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
!
    character(len=*), intent(in)        :: elrefz
    integer(kind=8), optional, intent(out)      :: nno, ndim, nnos
    real(kind=8), optional, intent(out) :: nodeCoor(3, MT_NNOMAX), cellVolu
!
! --------------------------------------------------------------------------------------------------
!
! Finite elements management
!
! Get main parameters of geometric support in parametric space
!
! --------------------------------------------------------------------------------------------------
!
! In  elrefa           : name of geometric support
! Out ndim             : topological dimension (0/1/2/3)
! Out nno              : number of nodes
! Out nnos             : number of geometric vertex nodes
! Out nodeCoor         : coordinates of node of geometric support in parametric space
! Out cellVolu         : volume of geometric support in parametric space
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nnos_, ndim_, nno_
    real(kind=8), parameter :: untiers = 1.d0/3.d0, deuxtiers = 2.d0/3.d0
    real(kind=8) :: cellVolu_
!
! --------------------------------------------------------------------------------------------------
!
    nnos_ = 0
    ndim_ = 0
    nno_ = 0
    cellVolu_ = 0.d0
!
    if (present(nodeCoor)) then
        nodeCoor = 0.d0
    end if
!
    select case (elrefz)
    case ('HE8')
        nno_ = 8
        nnos_ = 8
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:8) = [-1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0, -1.d0]
            nodeCoor(2, 1:8) = [-1.d0, -1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0]
            nodeCoor(3, 1:8) = [-1.d0, -1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0, +1.d0]
        end if
        cellVolu_ = 8.d0

    case ('HE9')
        nno_ = 9
        nnos_ = 8
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:9) = [-1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0, -1.d0, 0.d0]
            nodeCoor(2, 1:9) = [-1.d0, -1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0, 0.d0]
            nodeCoor(3, 1:9) = [-1.d0, -1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0, +1.d0, 0.d0]
        end if
        cellVolu_ = 8.d0

    case ('H10')
        nno_ = 27
        ! ne sert que pour passer dans inmat5 sans rien faire
        nnos_ = 27
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:8) = [-1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0, -1.d0]
            nodeCoor(2, 1:8) = [-1.d0, -1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0]
            nodeCoor(3, 1:8) = [-1.d0, -1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(1, 9:20) = [0.d0, +1.d0, 0.d0, -1.d0, -1.d0, +1.d0, +1.d0, -1.d0, &
                                 0.d0, +1.d0, 0.d0, -1.d0]
            nodeCoor(2, 9:20) = [-1.d0, 0.d0, +1.d0, 0.d0, -1.d0, -1.d0, +1.d0, +1.d0, &
                                 -1.d0, 0.d0, +1.d0, 0.d0]
            nodeCoor(3, 9:20) = [-1.d0, -1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                                 +1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(1, 21:27) = [0.d0, 0.d0, +1.d0, 0.d0, -1.d0, 0.d0, 0.d0]
            nodeCoor(2, 21:27) = [0.d0, -1.d0, 0.d0, +1.d0, 0.d0, 0.d0, 0.d0]
            nodeCoor(3, 21:27) = [-1.d0, 0.d0, 0.d0, 0.d0, 0.d0, +1.d0, 0.d0]
        end if
        cellVolu_ = 8.d0

    case ('H20')
        nno_ = 20
        nnos_ = 8
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:8) = [-1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0, -1.d0]
            nodeCoor(2, 1:8) = [-1.d0, -1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0]
            nodeCoor(3, 1:8) = [-1.d0, -1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(1, 9:20) = [0.d0, +1.d0, 0.d0, -1.d0, -1.d0, +1.d0, +1.d0, -1.d0, &
                                 0.d0, +1.d0, 0.d0, -1.d0]
            nodeCoor(2, 9:20) = [-1.d0, 0.d0, +1.d0, 0.d0, -1.d0, -1.d0, +1.d0, +1.d0, &
                                 -1.d0, 0.d0, +1.d0, 0.d0]
            nodeCoor(3, 9:20) = [-1.d0, -1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                                 +1.d0, +1.d0, +1.d0, +1.d0]
        end if
        cellVolu_ = 8.d0

    case ('H27')
        nno_ = 27
        nnos_ = 8
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:8) = [-1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0, -1.d0]
            nodeCoor(2, 1:8) = [-1.d0, -1.d0, +1.d0, +1.d0, -1.d0, -1.d0, +1.d0, +1.d0]
            nodeCoor(3, 1:8) = [-1.d0, -1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(1, 9:20) = [0.d0, +1.d0, 0.d0, -1.d0, -1.d0, +1.d0, +1.d0, -1.d0, &
                                 0.d0, +1.d0, 0.d0, -1.d0]
            nodeCoor(2, 9:20) = [-1.d0, 0.d0, +1.d0, 0.d0, -1.d0, -1.d0, +1.d0, +1.d0, &
                                 -1.d0, 0.d0, +1.d0, 0.d0]
            nodeCoor(3, 9:20) = [-1.d0, -1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                                 +1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(1, 21:27) = [0.d0, 0.d0, +1.d0, 0.d0, -1.d0, 0.d0, 0.d0]
            nodeCoor(2, 21:27) = [0.d0, -1.d0, 0.d0, +1.d0, 0.d0, 0.d0, 0.d0]
            nodeCoor(3, 21:27) = [-1.d0, 0.d0, 0.d0, 0.d0, 0.d0, +1.d0, 0.d0]
        end if
        cellVolu_ = 8.d0

    case ('TE4')
        nno_ = 4
        nnos_ = 4
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:4) = [0.d0, 0.d0, 0.d0, +1.d0]
            nodeCoor(2, 1:4) = [+1.d0, 0.d0, 0.d0, 0.d0]
            nodeCoor(3, 1:4) = [0.d0, +1.d0, 0.d0, 0.d0]
        end if
        cellVolu_ = 1.d0/6.d0

    case ('T10')
        nno_ = 10
        nnos_ = 4
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:4) = [0.d0, 0.d0, 0.d0, +1.d0]
            nodeCoor(2, 1:4) = [+1.d0, 0.d0, 0.d0, 0.d0]
            nodeCoor(3, 1:4) = [0.d0, +1.d0, 0.d0, 0.d0]
            nodeCoor(1, 5:10) = [0.d0, 0.d0, 0.d0, 0.5d0, 0.5d0, 0.5d0]
            nodeCoor(2, 5:10) = [0.5d0, 0.d0, 0.5d0, 0.5d0, 0.d0, 0.d0]
            nodeCoor(3, 5:10) = [0.5d0, 0.5d0, 0.d0, 0.d0, 0.5d0, 0.d0]
        end if
        cellVolu_ = 1.d0/6.d0

    case ('T15')
        nno_ = 15
        nnos_ = 4
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:4) = [0.d0, 0.d0, 0.d0, +1.d0]
            nodeCoor(2, 1:4) = [+1.d0, 0.d0, 0.d0, 0.d0]
            nodeCoor(3, 1:4) = [0.d0, +1.d0, 0.d0, 0.d0]
            nodeCoor(1, 5:10) = [0.d0, 0.d0, 0.d0, 0.5d0, 0.5d0, 0.5d0]
            nodeCoor(2, 5:10) = [0.5d0, 0.d0, 0.5d0, 0.5d0, 0.d0, 0.d0]
            nodeCoor(3, 5:10) = [0.5d0, 0.5d0, 0.d0, 0.d0, 0.5d0, 0.d0]
            nodeCoor(1, 11:14) = [0.d0, untiers, untiers, untiers]
            nodeCoor(2, 11:14) = [untiers, untiers, untiers, 0.d0]
            nodeCoor(3, 11:14) = [untiers, untiers, 0.d0, untiers]
            nodeCoor(1:3, 15) = [0.25d0, 0.25d0, 0.25d0]
        end if
        cellVolu_ = 1.d0/6.d0

    case ('T20')
        nno_ = 20
        nnos_ = 4
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1:3, 1) = [0.d0, 1.d0, 0.d0]
            nodeCoor(1:3, 2) = [0.d0, 0.d0, 1.d0]
            nodeCoor(1:3, 3) = [0.d0, 0.d0, 0.d0]
            nodeCoor(1:3, 4) = [1.d0, 0.d0, 0.d0]
            nodeCoor(1:3, 5) = [0.d0, deuxtiers, untiers]
            nodeCoor(1:3, 6) = [0.d0, untiers, deuxtiers]
            nodeCoor(1:3, 7) = [0.d0, 0.d0, deuxtiers]
            nodeCoor(1:3, 8) = [0.d0, 0.d0, untiers]
            nodeCoor(1:3, 9) = [0.d0, deuxtiers, 0.d0]
            nodeCoor(1:3, 10) = [0.d0, untiers, 0.d0]
            nodeCoor(1:3, 11) = [untiers, deuxtiers, 0.d0]
            nodeCoor(1:3, 12) = [deuxtiers, untiers, 0.d0]
            nodeCoor(1:3, 13) = [untiers, 0.d0, deuxtiers]
            nodeCoor(1:3, 14) = [deuxtiers, 0.d0, untiers]
            nodeCoor(1:3, 15) = [untiers, 0.d0, 0.d0]
            nodeCoor(1:3, 16) = [deuxtiers, 0.d0, 0.d0]
            nodeCoor(1:3, 17) = [0.d0, untiers, untiers]
            nodeCoor(1:3, 18) = [untiers, untiers, untiers]
            nodeCoor(1:3, 19) = [untiers, untiers, 0.d0]
            nodeCoor(1:3, 20) = [untiers, 0.d0, untiers]
        end if
        cellVolu_ = 1.d0/6.d0

    case ('PE6')
        nno_ = 6
        nnos_ = 6
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:6) = [-1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(2, 1:6) = [+1.d0, 0.d0, 0.d0, +1.d0, 0.d0, 0.d0]
            nodeCoor(3, 1:6) = [0.d0, +1.d0, 0.d0, 0.d0, +1.d0, 0.d0]
        end if
        cellVolu_ = 1.d0

    case ('PE7')
        nno_ = 7
        nnos_ = 6
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:7) = [0.d0, 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, untiers]
            nodeCoor(2, 1:7) = [0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 1.d0, untiers]
            nodeCoor(3, 1:7) = [-1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0, 0.d0]
        end if
        cellVolu_ = 1.d0

    case ('P15')
        nno_ = 15
        nnos_ = 6
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:6) = [-1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(2, 1:6) = [+1.d0, 0.d0, 0.d0, +1.d0, 0.d0, 0.d0]
            nodeCoor(3, 1:6) = [0.d0, +1.d0, 0.d0, 0.d0, +1.d0, 0.d0]
            nodeCoor(1, 7:15) = [-1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(2, 7:15) = [0.5d0, 0.d0, 0.5d0, +1.d0, 0.d0, 0.d0, 0.5d0, 0.d0, 0.5d0]
            nodeCoor(3, 7:15) = [0.5d0, 0.5d0, 0.d0, 0.d0, +1.d0, 0.d0, 0.5d0, 0.5d0, 0.d0]
        end if
        cellVolu_ = 1.d0

    case ('P18')
        nno_ = 18
        nnos_ = 6
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:6) = [-1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(2, 1:6) = [+1.d0, 0.d0, 0.d0, +1.d0, 0.d0, 0.d0]
            nodeCoor(3, 1:6) = [0.d0, +1.d0, 0.d0, 0.d0, +1.d0, 0.d0]
            nodeCoor(1, 7:15) = [-1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(2, 7:15) = [0.5d0, 0.d0, 0.5d0, +1.d0, 0.d0, 0.d0, 0.5d0, 0.d0, 0.5d0]
            nodeCoor(3, 7:15) = [0.5d0, 0.5d0, 0.d0, 0.d0, +1.d0, 0.d0, 0.5d0, 0.5d0, 0.d0]
            nodeCoor(1, 16:18) = [0.d0, 0.d0, 0.d0]
            nodeCoor(2, 16:18) = [0.5d0, 0.d0, 0.5d0]
            nodeCoor(3, 16:18) = [0.5d0, 0.5d0, 0.d0]
        end if
        cellVolu_ = 1.d0

    case ('P21')
        nno_ = 21
        nnos_ = 6
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:6) = [-1.d0, -1.d0, -1.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(2, 1:6) = [+1.d0, 0.d0, 0.d0, +1.d0, 0.d0, 0.d0]
            nodeCoor(3, 1:6) = [0.d0, +1.d0, 0.d0, 0.d0, +1.d0, 0.d0]
            nodeCoor(1, 7:15) = [-1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, +1.d0, +1.d0, +1.d0]
            nodeCoor(2, 7:15) = [0.5d0, 0.d0, 0.5d0, +1.d0, 0.d0, 0.d0, 0.5d0, 0.d0, 0.5d0]
            nodeCoor(3, 7:15) = [0.5d0, 0.5d0, 0.d0, 0.d0, +1.d0, 0.d0, 0.5d0, 0.5d0, 0.d0]
            nodeCoor(1, 16:20) = [0.d0, 0.d0, 0.d0, -1.d0, 1.d0]
            nodeCoor(2, 16:20) = [0.5d0, 0.d0, 0.5d0, untiers, untiers]
            nodeCoor(3, 16:20) = [0.5d0, 0.5d0, 0.d0, untiers, untiers]
            nodeCoor(1:3, 21) = [0.d0, untiers, untiers]
        end if
        cellVolu_ = 1.d0

    case ('PY5')
        nno_ = 5
        nnos_ = 5
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:5) = [+1.d0, 0.d0, -1.d0, 0.d0, 0.d0]
            nodeCoor(2, 1:5) = [0.d0, +1.d0, 0.d0, -1.d0, 0.d0]
            nodeCoor(3, 1:5) = [0.d0, 0.d0, 0.d0, 0.d0, +1.d0]
        end if
        cellVolu_ = 2.d0/3.d0

    case ('P13')
        nno_ = 13
        nnos_ = 5
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:5) = [+1.d0, 0.d0, -1.d0, 0.d0, 0.d0]
            nodeCoor(2, 1:5) = [0.d0, +1.d0, 0.d0, -1.d0, 0.d0]
            nodeCoor(3, 1:5) = [0.d0, 0.d0, 0.d0, 0.d0, +1.d0]
            nodeCoor(1, 6:13) = [0.5d0, -0.5d0, -0.5d0, 0.5d0, 0.5d0, 0.d0, -0.5d0, 0.d0]
            nodeCoor(2, 6:13) = [0.5d0, 0.5d0, -0.5d0, -0.5d0, 0.d0, 0.5d0, 0.d0, -0.5d0]
            nodeCoor(3, 6:13) = [0.d0, 0.d0, 0.d0, 0.d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0]
        end if
        cellVolu_ = 2.d0/3.d0

    case ('P19')
        nno_ = 19
        nnos_ = 5
        ndim_ = 3
        if (present(nodeCoor)) then
            nodeCoor(1, 1:5) = [+1.d0, 0.d0, -1.d0, 0.d0, 0.d0]
            nodeCoor(2, 1:5) = [0.d0, +1.d0, 0.d0, -1.d0, 0.d0]
            nodeCoor(3, 1:5) = [0.d0, 0.d0, 0.d0, 0.d0, +1.d0]
            nodeCoor(1, 6:13) = [0.5d0, -0.5d0, -0.5d0, 0.5d0, 0.5d0, 0.d0, -0.5d0, 0.d0]
            nodeCoor(2, 6:13) = [0.5d0, 0.5d0, -0.5d0, -0.5d0, 0.d0, 0.5d0, 0.d0, -0.5d0]
            nodeCoor(3, 6:13) = [0.d0, 0.d0, 0.d0, 0.d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0]
            nodeCoor(1, 14:18) = [0.d0, untiers, -untiers, -untiers, untiers]
            nodeCoor(2, 14:18) = [0.d0, untiers, untiers, -untiers, -untiers]
            nodeCoor(3, 14:18) = [0.d0, untiers, untiers, untiers, untiers]
            nodeCoor(1:3, 19) = [0.d0, 0.d0, 0.2d0]
        end if
        cellVolu_ = 2.d0/3.d0

    case ('TR3')
        nno_ = 3
        nnos_ = 3
        ndim_ = 2
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [0.d0, +1.d0, 0.d0]
            nodeCoor(2, 1:nno_) = [0.d0, 0.d0, +1.d0]
        end if
        cellVolu_ = 1.d0/2.d0

    case ('TR6')
        nno_ = 6
        nnos_ = 3
        ndim_ = 2
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [0.d0, +1.d0, 0.d0, 0.5d0, 0.5d0, 0.d0]
            nodeCoor(2, 1:nno_) = [0.d0, 0.d0, +1.d0, 0.d0, 0.5d0, 0.5d0]
        end if
        cellVolu_ = 1.d0/2.d0

    case ('TR7')
        nno_ = 7
        nnos_ = 3
        ndim_ = 2
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [0.d0, +1.d0, 0.d0, 0.5d0, 0.5d0, 0.d0, untiers]
            nodeCoor(2, 1:nno_) = [0.d0, 0.d0, +1.d0, 0.d0, 0.5d0, 0.5d0, untiers]
        end if
        cellVolu_ = 1.d0/2.d0

    case ('TR1')
        nno_ = 10
        nnos_ = 3
        ndim_ = 2
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [0.d0, +1.d0, 0.d0, untiers, 2.d0/3.d0, 2.d0/3.d0, &
                                   untiers, 0.d0, 0.d0, untiers]
            nodeCoor(2, 1:nno_) = [0.d0, 0.d0, +1.d0, 0.d0, 0.0d0, untiers, &
                                   2.d0/3.d0, 2.d0/3.d0, untiers, untiers]
        end if
        cellVolu_ = 1.d0/2.d0

    case ('QU4')
        nno_ = 4
        nnos_ = 4
        ndim_ = 2
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [-1.d0, +1.d0, +1.d0, -1.d0]
            nodeCoor(2, 1:nno_) = [-1.d0, -1.d0, +1.d0, +1.d0]
        end if
        cellVolu_ = 4.d0

    case ('QU8')
        nno_ = 8
        nnos_ = 4
        ndim_ = 2
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [-1.d0, +1.d0, +1.d0, -1.d0, 0.d0, +1.d0, 0.d0, -1.d0]
            nodeCoor(2, 1:nno_) = [-1.d0, -1.d0, +1.d0, +1.d0, -1.d0, 0.d0, +1.d0, 0.d0]
        end if
        cellVolu_ = 4.d0

    case ('QU9')
        nno_ = 9
        nnos_ = 4
        ndim_ = 2
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [-1.d0, +1.d0, +1.d0, -1.d0, 0.d0, +1.d0, 0.d0, -1.d0, 0.d0]
            nodeCoor(2, 1:nno_) = [-1.d0, -1.d0, +1.d0, +1.d0, -1.d0, 0.d0, +1.d0, 0.d0, 0.d0]
        end if
        cellVolu_ = 4.d0

    case ('Q12')
        nno_ = 12
        nnos_ = 4
        ndim_ = 2
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [-1.d0, +1.d0, +1.d0, -1.d0, -untiers, untiers, +1.d0, +1.d0, &
                                   untiers, -untiers, -1.d0, -1.d0]
            nodeCoor(2, 1:nno_) = [-1.d0, -1.d0, +1.d0, +1.d0, -1.d0, -1.d0, -untiers, untiers, &
                                   +1.d0, +1.d0, untiers, -untiers]
        end if
        cellVolu_ = 4.d0

    case ('SE2')
        nno_ = 2
        nnos_ = 2
        ndim_ = 1
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [-1.d0, +1.d0]
        end if
        cellVolu_ = 2.d0

    case ('SE3')
        nno_ = 3
        nnos_ = 2
        ndim_ = 1
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [-1.d0, +1.d0, 0.d0]
        end if
        cellVolu_ = 2.d0

    case ('SE4')
        nno_ = 4
        nnos_ = 2
        ndim_ = 1
        if (present(nodeCoor)) then
            nodeCoor(1, 1:nno_) = [-1.d0, +1.d0, -1.d0/3.d0, 1.d0/3.d0]
        end if
        cellVolu_ = 2.d0

    case ('PO1')
        nno_ = 1
        nnos_ = 1
        ndim_ = 0
        if (present(nodeCoor)) then
            nodeCoor(1, 1) = 0.d0
        end if
        cellVolu_ = 1.d0

    case default
        ASSERT(ASTER_FALSE)

    end select
!
    if (present(nno)) then
        nno = nno_
    end if
    if (present(nnos)) then
        nnos = nnos_
    end if
    if (present(ndim)) then
        ndim = ndim_
    end if
    if (present(cellVolu)) then
        cellVolu = cellVolu_
    end if
!
end subroutine
