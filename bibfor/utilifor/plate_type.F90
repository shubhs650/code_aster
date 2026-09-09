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
! ==================================================================================================
!
! Types for the plates/shells elements
!
! ==================================================================================================
!
module plate_type
! ==================================================================================================

! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/plate_type.h"
! ==================================================================================================
! Global variables
! ==================================================================================================
! ==================================================================================================
! Type: parameters
! ==================================================================================================
    type plateCara_Para
! ----- Parameters have been read from user
        aster_logical :: lRead = ASTER_FALSE
        integer(kind=8) :: type = PLATE_UNKW
        integer(kind=8) :: geom = PLATE_GEOM_UNKW
        integer(kind=8) :: nbLayer = 0
        real(kind=8) :: thick = 0.d0
        real(kind=8) :: shearCoef = 0.d0
        real(kind=8) :: metric = 0.d0
        real(kind=8) :: coefRigiDRZ = 0.d0
        real(kind=8) :: offset = 0.d0
        real(kind=8) :: inerRota = 0.d0
        real(kind=8) :: section = 0.d0
        real(kind=8) :: tension = 0.d0
    end type plateCara_Para
! ==================================================================================================
! Type: orientation
! ==================================================================================================
    type plateOrie_Para
! ----- Local orientation has been read from user
        aster_logical :: lRead = ASTER_FALSE
! ----- Orientation is updated
        aster_logical :: lUpdate = ASTER_FALSE
! ----- Orientation
        real(kind=8) :: alpha = 0.d0, beta = 0.d0
! ----- For COQUE_3D: normal and tangents
        real(kind=8) :: vectNorm(9, 3) = 0.d0
        real(kind=8) :: vectTang(9, 2, 3) = 0.d0
        ! From Local (I) to global (U)
        real(kind=8) :: matevn(2, 2, 10) = 0.d0
        real(kind=8) :: matevg(2, 2, 10) = 0.d0
! ----- For plates
        real(kind=8) :: t2iu(4) = 0.d0, t2ui(4) = 0.d0
        real(kind=8) :: c = 0.d0, s = 0.d0
        real(kind=8) :: t1ve(9) = 0.d0
! ----- For grids
        real(kind=8) :: gridDir11(3) = 0.d0
        real(kind=8) :: gridNorm(3) = 0.d0
    end type plateOrie_Para
!===================================================================================================
    public :: plateCara_Para, plateOrie_Para
contains
!===================================================================================================
end module
