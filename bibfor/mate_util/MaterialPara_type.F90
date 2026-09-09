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
! Types for the type of material parameters in behaviour
!
! ==================================================================================================
!
module MaterialPara_type
! ==================================================================================================

! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/MaterialPara_type.h"
! ==================================================================================================
! Global variables
! ==================================================================================================
! ==================================================================================================
! Type: parameters of local coordinate system
! ==================================================================================================
    type LCS_Para
! ----- Type
        integer(kind=8) :: lcsType = MATER_LCS_UNDEF
! ----- Angles
        real(kind=8) :: lcsAngle(3) = 0.d0
! ----- Angles at Gauss points
        real(kind=8) :: lcsAnglePg(3*27) = 0.d0
    end type LCS_Para
! ==================================================================================================
! Type: parameters of integration scheme
! ==================================================================================================
    type Scheme_Para
! ----- Name of quadrature scheme
        character(len=8) :: fami = " "
! ----- Index of quadrature point
        integer(kind=8) :: kpg = 0
! ----- Index of "sub"-point (plates, pipes, beams, etc.)
        integer(kind=8) :: ksp = 0
    end type Scheme_Para
! ==================================================================================================
! Type: Parameters for material parametres (main object)
! ==================================================================================================
    type Material_Para
! ----- Parameters of integrate scheme
        type(Scheme_Para) :: schemePara

! ----- Parameters of local coordinate system
        type(LCS_Para) :: lcsPara

! ----- Adress for material parameters
        integer(kind=8) :: jvMaterCode = 0

! ----- Material name
        character(len=8) :: matname = " "

! ----- Type for elasticity
        integer(kind=8) :: elasID = ELAS_UNDEF
        character(len=16) :: elasKeyword = " "

! ----- Flag for metallurgical case
        aster_logical :: lElasIsMeta = ASTER_FALSE
        aster_logical :: lMetaLemaAni = ASTER_FALSE

    end type Material_Para
!===================================================================================================
    public :: Material_Para, Scheme_Para, LCS_Para
contains
!===================================================================================================
end module
