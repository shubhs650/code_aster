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
! Module for types in MODI_MAILLAGE command
!
! ==================================================================================================
module mesh_operators_type
! ==================================================================================================
    use Behaviour_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
! ==================================================================================================
! Global variables
! ==================================================================================================
! - Non
! ==================================================================================================
! Define types
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
! For MODI_MAILLAGE/ ORIE_NORM_COQUE operator
! --------------------------------------------------------------------------------------------------
    type MESH_OPER_ORIE_SHELL
        aster_logical :: orieByVect = ASTER_FALSE
        real(kind=8) :: orieVect(3) = 0.d0
        integer(kind=8) :: nbGroupCell = 0
        character(len=24), pointer :: listOfGroupOfCell(:) => null()
        integer(kind=8) :: nodeNume = 0
    end type MESH_OPER_ORIE_SHELL
! --------------------------------------------------------------------------------------------------
! For MODI_MAILLAGE/ ORIE_INTERF_POU operator
! --------------------------------------------------------------------------------------------------
    type MESH_OPER_ORIE_INTERF_POU
        real(kind=8) :: orieVect(3) = 0.d0
        integer(kind=8) :: nbGroupCell = 0
        character(len=24), pointer :: listOfGroupOfCell(:) => null()
    end type MESH_OPER_ORIE_INTERF_POU
! --------------------------------------------------------------------------------------------------
! For MODI_MAILLAGE/ORIE_* operators
! --------------------------------------------------------------------------------------------------
    type MESH_OPER_MODI_PARA
        integer(kind=8) :: orieShell = 0
        type(MESH_OPER_ORIE_SHELL), pointer :: meshOperOrieShell(:) => null()
        integer(kind=8) :: orieInterfPou = 0
        type(MESH_OPER_ORIE_INTERF_POU), pointer :: meshOperOrieInterfPou(:) => null()
        integer(kind=8) :: orieSkin = 0
        integer(kind=8) :: orieLine = 0
    end type MESH_OPER_MODI_PARA
!===================================================================================================
!===================================================================================================
    public :: MESH_OPER_ORIE_SHELL, MESH_OPER_ORIE_INTERF_POU, MESH_OPER_MODI_PARA
contains
!===================================================================================================
!===================================================================================================
!
end module mesh_operators_type
