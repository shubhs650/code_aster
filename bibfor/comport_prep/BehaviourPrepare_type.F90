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
! Types for preparation of behaviour
!
! ==================================================================================================
!
module BehaviourPrepare_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
! ==================================================================================================
! Type: parameters for external behaviours (MFront/UMAT)
! ==================================================================================================
    type BehaviourPrep_Exte
! ----- Type of external integrator
        integer(kind=8) :: solvBehavType = SOLV_BEHAV_ASTER
! ----- Address to UMAT function
        integer(kind=8) :: adrsUMAT = 0
! ----- Name of subroutine for external UMAT law
        character(len=255) :: subrNameUMAT = ' '
! ----- Name of library for external UMAT law
        character(len=255) :: librNameUMAT = ' '
! ----- Number of internal variables for UMAT
        integer(kind=8) :: nbVariUMAT = 0
! ----- Address to MGISBehaviour object as hexadecimal of integer
        character(len=16) :: adrsMGIS = ' '
        integer(kind=8) :: adrsMGISInteger = 0
! ----- Model for MFront law
        integer(kind=8) :: modelMGIS = MGIS_MODEL_UNSET
! ----- Identifier for strains model
        integer(kind=8) :: strainMGIS = MGIS_STRAIN_UNSET
! ----- stress plane hypothesis (for Deborst)
        character(len=16) :: cplaMGIS = " "
    end type BehaviourPrep_Exte

! ==================================================================================================
! Type: behaviour parameters from user
! ==================================================================================================
    type BehaviourPrep_Para
! ----- Keyword RELATION
        character(len=16) :: rela_comp = ' '
! ----- Keyword DEFORMATION
        character(len=16) :: defo_comp = ' '
! ----- Keyword COMP_INCR/COMP_ELAS
        character(len=16) :: type_comp = ' '
! ----- Keyword DEBORST
        character(len=16) :: type_cpla = ' '
! ----- Keyword KIT
        character(len=16) :: kit_comp(4) = ' '
! ----- Keyword COMPOR
        character(len=16) :: mult_comp = ' '
! ----- Keyword POST_ITER
        character(len=16) :: post_iter = ' '
! ----- Type of strain transmitted to the behaviour law : 'OLD', 'MECANIQUE' or 'TOTALE'
        character(len=16) :: defo_ldc = ' '
! ----- Index of law
        integer(kind=8) :: numeLaw = 0
! ----- Total number of internal state variables
        integer(kind=8) :: nbVari = 0
! ----- Number of internal state variables for kit
        integer(kind=8) :: nbVariKit(4) = 0
! ----- Index of law for kit
        integer(kind=8) :: numeLawKit(4) = 0
! ----- Keyword RIGI_GEOM
        character(len=16) :: rigi_geom = ' '
! ----- Keyword REGU_VISC
        character(len=16) :: regu_visc = ' '
! ----- Mechanical part of behaviour
        character(len=16) :: meca_comp = ' '
! ----- Keyword POST_INCR
        character(len=16) :: post_incr = ' '
! ----- Flag for total strain model cases
        aster_logical :: lTotalStrain = ASTER_FALSE
    end type BehaviourPrep_Para
! ==================================================================================================
! Type: behaviour criteria from user
! ==================================================================================================
    type BehaviourPrep_Crit
! ----- Keyword RELATION
        character(len=16) :: rela_comp = ' '
! ----- Mechanical part of behaviour
        character(len=16) :: meca_comp = ' '
! ----- Parameters for external behaviours
        type(BehaviourPrep_Exte) :: prepExte
! ----- Criteria
        integer(kind=8) :: type_matr_t = 0
        real(kind=8) :: parm_theta = 0.d0
        integer(kind=8) :: iter_inte_pas = 0
        real(kind=8) :: vale_pert_rela = 0.d0
        real(kind=8) :: resi_deborst_max = 0.d0
        integer(kind=8) :: iter_deborst_max = 0
        real(kind=8) :: resi_radi_rela = 0.d0
        integer(kind=8) :: ipostiter = 0
        integer(kind=8) :: iveriborne = 0
        aster_logical :: l_matr_unsymm = ASTER_FALSE
        real(kind=8) :: algo_inte_r = 0.d0
        real(kind=8), pointer :: resi_inte => null()
        integer(kind=8), pointer :: iter_inte_maxi => null()
        integer(kind=8) :: jvariext1 = 0
        integer(kind=8) :: jvariext2 = 0
    end type BehaviourPrep_Crit
! ==================================================================================================
! Type: map for criteria of behaviours (CARCRI)
! ==================================================================================================
    type BehaviourPrep_MapCarcri
! ----- Number of factor keywords
        integer(kind=8) :: nb_comp = 0
! ----- Parameters for THM scheme
        real(kind=8) :: parm_alpha_thm = 0.d0
        real(kind=8) :: parm_theta_thm = 0.d0
! ----- List of criteria (by keyword COMPORTEMENT)
        type(BehaviourPrep_Crit), pointer :: prepCrit(:)
    end type BehaviourPrep_MapCarcri
! ==================================================================================================
! Type: map for parameters of behaviours (COMPOR)
! ==================================================================================================
    type BehaviourPrep_MapCompor
! ----- Number of factor keywords
        integer(kind=8) :: nb_comp = 0
! ----- List of parameters
        type(BehaviourPrep_Para), pointer :: prepPara(:) => null()
! ----- List of parameters for external behaviours
        type(BehaviourPrep_Exte), pointer :: prepExte(:) => null()
! ----- Flag for debug
        aster_logical :: lDebug = ASTER_FALSE
    end type BehaviourPrep_MapCompor
!===================================================================================================
    public :: BehaviourPrep_Exte, BehaviourPrep_MapCarcri, BehaviourPrep_MapCompor
    public :: BehaviourPrep_Crit, BehaviourPrep_Para
contains
!===================================================================================================
end module
