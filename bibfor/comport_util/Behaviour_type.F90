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
! Types for the management of integration of behaviour
!
! ==================================================================================================
!
module Behaviour_type
! ==================================================================================================
    use BehaviourStrain_type
    use MaterialPara_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "MeshTypes_type.h"
! ==================================================================================================
! Global variables
! ==================================================================================================
! ==================================================================================================
! Type: parameters of behaviour
! ==================================================================================================
    type Behaviour_Para
! ----- Times
        real(kind=8) :: timePrev = 0.d0
        real(kind=8) :: timeCurr = 0.d0
! ----- Type of strain decomposition
        aster_logical :: lStrainMeca = ASTER_FALSE
        aster_logical :: lStrainAll = ASTER_FALSE
        aster_logical :: lStrainOld = ASTER_FALSE
! ----- Viscuous model
        aster_logical :: lReguVisc = ASTER_FALSE
! ----- Flags for finite transformation
        aster_logical :: lFiniteStrain = ASTER_FALSE
        aster_logical :: lGdefLog = ASTER_FALSE
! ----- Flag for annealing
        aster_logical :: lAnnealing = ASTER_FALSE
! ----- Flag for non-linear options
        aster_logical :: lNonLinear = ASTER_FALSE
! ----- Phase: prediction
        aster_logical :: lPred = ASTER_FALSE
! ----- Part to compute: internal state variables
        aster_logical :: lVari = ASTER_FALSE
! ----- Part to compute: stress
        aster_logical :: lSigm = ASTER_FALSE
! ----- Part to compute: matrix
        aster_logical :: lMatr = ASTER_FALSE
! ----- Activation of IMPLEX method
        aster_logical :: lImplex = ASTER_FALSE
! ----- Offset for index of behaviour
        integer(kind=8) :: lawIndexOffset = 0
! ----- Automatic mechanism for the size of stress and strain tensors (ndimsi)
        aster_logical :: lNdimsiEnabled = ASTER_FALSE
        aster_logical :: lEpsComEnabled = ASTER_FALSE
        aster_logical :: lDeborstEnabled = ASTER_FALSE
        aster_logical :: lCtrlEpsEnabled = ASTER_FALSE
        aster_logical :: lReguViscEnabled = ASTER_FALSE
! ----- Flags for standard FE
        aster_logical :: lAxis = ASTER_FALSE
        aster_logical :: lThreeDim = ASTER_FALSE
        aster_logical :: lPlaneStress = ASTER_FALSE
        aster_logical :: lPlaneStrain = ASTER_FALSE
        aster_logical :: lTwoDim = ASTER_FALSE
        aster_logical :: lOneDim = ASTER_FALSE
! ----- Flags for other modeliazations
        aster_logical :: lGradVari = ASTER_FALSE
        aster_logical :: lTHM = ASTER_FALSE
! ----- Flags for external solves
        aster_logical :: lExteSolver = ASTER_FALSE
        aster_logical :: lMGIS = ASTER_FALSE
        aster_logical :: lUMAT = ASTER_FALSE
! ----- Index of law
        integer(kind=8) :: numlc = 0
! ----- Number of internal state variable
        integer(kind=8) :: nvi = 0
! ----- Default size of strain and stress tensors
        integer(kind=8) :: ndimsi = 0
! ----- Flag to check bounds
        aster_logical :: lChckBounds = ASTER_FALSE
! ----- Cut level
        integer(kind=8) :: cutLevel = 0
    end type Behaviour_Para
! ==================================================================================================
! Type: External state variables - geometric properties
! ==================================================================================================
    type BehaviourESVA_Geom
! ----- Flag for size of element for BETON_DOUBLE_DP
        aster_logical :: lElemSize1 = ASTER_FALSE
! ----- Size of element for BETON_DOUBLE_DP
        real(kind=8)  :: elemSize1 = 0.d0
! ----- Gradient of velocity for *CRISTAL
        real(kind=8)  :: gradVelo(9) = 0.d0
! ----- Coordinates of all Gauss points
        real(kind=8)  :: coorElga(VARC_GEOM_NBMAXI, 3) = 0.d0
    end type BehaviourESVA_Geom
! ==================================================================================================
! Type: External state variables - Other properties
! ==================================================================================================
    type BehaviourESVA_Other
! ----- For *_JOINT_HYME models : kinematic matrix
        real(kind=8) :: rotpg(3*3) = 0.d0
! ----- For CABLE_GAINE elements : tension of the cable
        real(kind=8) :: tenscab = 0.d0
! ----- For CABLE_GAINE elements : curvature of the cable
        real(kind=8) :: curvcab = 0.d0
! ----- For GRAD_VARI models : non-local variables PHI
        real(kind=8) :: nonloc(2) = 0.d0
! ----- For CZM_*_MIX behaviours : Lagrange penalty coefficient
        real(kind=8) :: r = 0.d0
! ----- Current time
        real(kind=8) :: time = 0.d0
! ----- Hygrometry
        aster_logical :: lHygr = ASTER_FALSE
        real(kind=8) :: hygrPrev = 0.d0
        real(kind=8) :: hygrIncr = 0.d0
    end type BehaviourESVA_Other
! ==================================================================================================
! Type: External state variables - External solver (UMAT, MFRONT/MGIS)
! ==================================================================================================
    type BehaviourESVA_Exte
! ----- Number of external state variables used (as scalar) in external solver
        integer(kind=8) :: nbESVAScal = 0
! ----- Value of external state variables used (as scalar) in external solver
        real(kind=8) :: scalESVAPrev(VARC_EXTE_NBMAXI) = 0.d0
! ----- Incremental value of external state variables used (as scalar)  in external solver
        real(kind=8) :: scalESVAIncr(VARC_EXTE_NBMAXI) = 0.d0
! ----- Address for library of MGIS
        character(len=16) :: mgisAddr = " "
    end type BehaviourESVA_Exte
! ==================================================================================================
! Type: Parameters for External State Variables
! ==================================================================================================
    type BehaviourESVA
! ----- Flag when GEOM external state variable is present
        aster_logical :: lGeomInESVA = ASTER_FALSE

! ----- Flag when inelastic strains are present
        aster_logical :: lhasInelasticStrains = ASTER_FALSE

! ----- Integer coded for presence of external state variables
        integer(kind=8) :: tabcod(60) = 0

! ----- Geometric properties
        type(BehaviourESVA_Geom) :: behavESVAGeom

! ----- Other properties
        type(BehaviourESVA_Other) :: behavESVAOther

! ----- External state variables for external solver (UMAT/MFRONT)
        type(BehaviourESVA_Exte) :: behavESVAExte

! ----- Non-mechanical strains (external state variable as temperature)
        real(kind=8) :: epsi_varc(6) = 0.d0
        real(kind=8) :: depsi_varc(6) = 0.d0
    end type BehaviourESVA
! ==================================================================================================
! Type: Parameters for integration (main object)
! ==================================================================================================
    type Behaviour_Integ
! ----- Material parameters
        type(Material_Para) :: materPara

! ----- Behaviour parameters
        type(Behaviour_Para) :: behavPara

! ----- Parameters for external state variables
        type(BehaviourESVA) :: behavESVA

! ----- All external state variables for anelastic strains
        type(All_Varc_Strain) :: allVarcStrain

    end type Behaviour_Integ
!===================================================================================================
    public :: Behaviour_Integ, BehaviourESVA
    public :: Behaviour_Para
    public :: BehaviourESVA_Geom, BehaviourESVA_Other, BehaviourESVA_Exte
contains
!===================================================================================================
end module
