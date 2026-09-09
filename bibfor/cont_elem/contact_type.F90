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
module contact_type
!
!
    implicit none
!
    private
!
#include "asterf_types.h"
#include "contact_module.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Type
!
! Generic type for contact
!
! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
!
    type ContactParameters
    !! Reso Geom parameters
        integer(kind=8)                             :: reso_geom = 0

    !! Contact parameters
        integer(kind=8)                             :: algo_cont = 0
        integer(kind=8)                             :: type_cont = 0
        integer(kind=8)                             :: vari_cont = 0
        integer(kind=8)                             :: jac_type = 0
        integer(kind=8)                             :: inte_type = 0
        real(kind=8)                        :: vari_cont_coef = 0.d0
        real(kind=8), dimension(MT_NNOMAX2D)          :: coef_cont = 0.d0

    !! Friction paramaters
        aster_logical                       :: l_fric = ASTER_FALSE
        integer(kind=8)                             :: algo_fric = 0
        integer(kind=8)                             :: type_fric = 0
        real(kind=8), dimension(MT_NNOMAX2D)          :: coef_fric = 0.d0
        real(kind=8), dimension(MT_NNOMAX2D)          :: threshold = 0.d0
        real(kind=8)                        :: threshold_given = 0.d0

    !! Other
        real(kind=8)                        :: proj_tole = 0.d0
        integer(kind=8)                             :: cont_init = PAIR_CONT_INTE
        real(kind=8)                        :: E = 0.d0
    end type
!
    type ContactGeom
    !! Slave side parameters
        integer(kind=8)                             :: nb_node_slav = 0
        character(len=8)                    :: elem_slav_code = " "
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: coor_slav_init = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: coor_slav_prev = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: coor_slav_curr = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: coor_slav_pair = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: depl_slav_curr = 0.d0
        real(kind=8), dimension(4)          :: lagc_slav_curr = 0.d0
        real(kind=8), dimension(2, 4)        :: lagf_slav_curr = 0.d0
        integer(kind=8)                             :: nb_lagr_c = 0
        integer(kind=8), dimension(MT_NNOMAX2D)               :: indi_lagc = 0

    !! Slave cell volu
        integer(kind=8)                             :: nb_node_volu = 0
        character(len=8)                    :: elem_volu_code = " "
        real(kind=8), dimension(3, MT_NNOMAX3D)        :: coor_volu_init = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX3D)        :: coor_volu_prev = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX3D)        :: coor_volu_curr = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX3D)        :: coor_volu_pair = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX3D)        :: depl_volu_curr = 0.d0
        integer(kind=8), dimension(MT_NNOMAX2D)        :: mapVolu2Surf = 0

    !! Master side paramaters
        integer(kind=8)                             :: nb_node_mast = 0
        character(len=8)                    :: elem_mast_code = " "
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: coor_mast_init = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: coor_mast_prev = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: coor_mast_curr = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: coor_mast_pair = 0.d0
        real(kind=8), dimension(3, MT_NNOMAX2D)        :: depl_mast_curr = 0.d0

    !! Time
        real(kind=8) :: time_prev = 0.d0, time_curr = 0.d0

    !! Other
        integer(kind=8)                             :: elem_dime = 0
        aster_logical                       :: l_axis = ASTER_FALSE
        integer(kind=8)                             :: nb_dofs = 0
    end type

    type Cell_GeomC
        integer(kind=8) :: cellDime = 0
        integer(kind=8) :: nbNode = 0
        character(len=8) :: cellCode = " "
        real(kind=8), dimension(3, MT_NNOMAX2D) :: coorNodeGlob = 0.d0
        real(kind=8), dimension(2, MT_NNOMAX2D) :: coorNodePara = 0.d0
    end type Cell_GeomC

    type Contact_CellGeom
        integer(kind=8) :: cellDime = 0
! ----- Slave side
        type(Cell_GeomC) :: slav
! ----- Master side
        type(Cell_GeomC) :: mast
    end type Contact_CellGeom

! - Type for parameters of projection algorithm
    type Contact_ProjAlgoPara
        aster_logical :: withPrepLine = ASTER_FALSE
        aster_logical :: newtDebug = ASTER_FALSE
        integer(kind=8) :: newtIterMaxi = 0
        real(kind=8) :: newtTole = 0.d0
    end type Contact_ProjAlgoPara

! - Type for parameters of projection
    type Contact_ProjPara
        real(kind=8) :: pointCoor(3) = 0.d0
        integer(kind=8) :: modelDime = 0
        type(Cell_GeomC) :: geomTarget
        type(Cell_GeomC) :: geomTargetLine
        real(kind=8):: projVect(3) = 0.d0
        real(kind=8) :: ksi(2) = 0.d0
        real(kind=8) :: tau1(3) = 0.d0
        real(kind=8) :: tau2(3) = 0.d0
        integer(kind=8):: errorCode = 0
    end type Contact_ProjPara

!
!===================================================================================================
!
!===================================================================================================
!
    public :: ContactParameters, ContactGeom
    public :: Cell_GeomC, Contact_CellGeom
    public :: Contact_ProjAlgoPara, Contact_ProjPara
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
end module
