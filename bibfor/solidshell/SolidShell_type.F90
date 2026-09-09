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
module SolidShell_type
!
    use Behaviour_type
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/SolidShell_type.h"

! ==================================================================================================
!
! Global variables - General
!
! ==================================================================================================

! Type of modelization for behaviour
    character(len=8), parameter :: typmod(2) = (/'3D      ', '        '/)

! Identity tensor in Voigt notation
    real(kind=8), parameter :: tensorIden(SSH_SIZE_TENS) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)

! Identity matrix for size 3
    real(kind=8), parameter :: matr3Iden(3, 3) = &
                               reshape((/1.d0, 0.d0, 0.d0, &
                                         0.d0, 1.d0, 0.d0, &
                                         0.d0, 0.d0, 1.d0/), (/3, 3/))

! ==================================================================================================
!
! Global variables for HEXA
!
! ==================================================================================================

! Vectors to define interpolations for HEXA
    real(kind=8), parameter :: hexaVectG1(8) = (/-1.d0, +1.d0, +1.d0, -1.d0, &
                                                 -1.d0, +1.d0, +1.d0, -1.d0/)/8.d0
    real(kind=8), parameter :: hexaVectG2(8) = (/-1.d0, -1.d0, +1.d0, +1.d0, &
                                                 -1.d0, -1.d0, +1.d0, +1.d0/)/8.d0
    real(kind=8), parameter :: hexaVectG3(8) = (/-1.d0, -1.d0, -1.d0, -1.d0, &
                                                 +1.d0, +1.d0, +1.d0, +1.d0/)/8.d0
    real(kind=8), parameter :: hexaVectH1(8) = (/+1.d0, -1.d0, +1.d0, -1.d0, &
                                                 1.d0, -1.d0, +1.d0, -1.d0/)/8.d0
    real(kind=8), parameter :: hexaVectH2(8) = (/+1.d0, 1.d0, -1.d0, -1.d0, &
                                                 -1.d0, -1.d0, +1.d0, +1.d0/)/8.d0
    real(kind=8), parameter :: hexaVectH3(8) = (/+1.d0, -1.d0, -1.d0, +1.d0, &
                                                 -1.d0, +1.d0, +1.d0, -1.d0/)/8.d0
    real(kind=8), parameter :: hexaVectH4(8) = (/-1.d0, +1.d0, -1.d0, 1.d0, &
                                                 +1.d0, -1.d0, +1.d0, -1.d0/)/8.d0
    real(kind=8), parameter :: hexaVectS1(8) = (/+1.d0, +1.d0, +1.d0, +1.d0, &
                                                 +1.d0, +1.d0, +1.d0, +1.d0/)/8.d0

! Center of cell in parametric frame for HEXA
    real(kind=8), parameter :: hexaCovaCenter(3) = (/0.d0, 0.d0, 0.d0/)

! Coordinates of rectangular faces for HEXA
    real(kind=8), parameter :: hexaQuadXIEH(4, 3) = &
                               reshape((/-1.d0, +1.d0, +1.d0, -1.d0, &
                                         +0.d0, +0.d0, +0.d0, +0.d0, &
                                         -1.d0, -1.d0, +1.d0, +1.d0/), (/4, 3/))
    real(kind=8), parameter :: hexaQuadXIAD(4, 3) = &
                               reshape((/-1.d0, +1.d0, +1.d0, -1.d0, &
                                         -1.d0, -1.d0, +1.d0, +1.d0, &
                                         +0.d0, +0.d0, +0.d0, +0.d0/), (/4, 3/))
    real(kind=8), parameter :: hexaQuadXIJM(4, 3) = &
                               reshape((/+0.d0, +0.d0, +0.d0, +0.d0, &
                                         -1.d0, +1.d0, +1.d0, -1.d0, &
                                         -1.d0, -1.d0, +1.d0, +1.d0/), (/4, 3/))

! Parameter for minimum coefficient of stabilization for HEXA
    real(kind=8), parameter :: hexaStabMini = 0.001d0

! D matrix for stabilization for HEXA
    real(kind=8), parameter :: hexaStabDMatr(6, 6) = &
                               reshape((/+4.d0/3.d0, -2.d0/3.d0, -2.d0/3.d0, 0.d0, 0.d0, 0.d0, &
                                         -2.d0/3.d0, +4.d0/3.d0, -2.d0/3.d0, 0.d0, 0.d0, 0.d0, &
                                         -2.d0/3.d0, -2.d0/3.d0, +4.d0/3.d0, 0.d0, 0.d0, 0.d0, &
                                         0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
                                         0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
                                         0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/6, 6/))

! Coefficient for stabilization for HEXA
    real(kind=8), parameter :: hexaStabCoef1 = 8.d0/3.d0
    real(kind=8), parameter :: hexaStabCoef2 = 8.d0/9.d0

! ==================================================================================================
!
! Derivated types - General
!
! ==================================================================================================

! Type to define properties of integration scheme
    type SSH_ELEM_INTE
! - Name of family for integration scheme
        character(len=4) :: inteFami = ' '
! - Number of integration points
        integer(kind=8)          :: nbIntePoint = 0
! - JEVEUX adress to weight of integrations points
        integer(kind=8)          :: jvWeight = 0
! - JEVEUX adress to coordinates of integration points
        integer(kind=8)          :: jvCoor = 0
! - JEVEUX adress to shape functions at integration points
        integer(kind=8)          :: jvShape = 0
! - JEVEUX adress to derivative of shape functions at integration points
        integer(kind=8)          :: jvDShape = 0
    end type SSH_ELEM_INTE

! Type to define properties of material
    type SSH_MATE_PARA
! ----- Parameters of material
        type(Material_Para) :: materPara
! ----- Elasticity matrix at center of element
        real(kind=8) :: elemHookeMatrix(SSH_SIZE_TENS, SSH_SIZE_TENS) = 0.d0
    end type SSH_MATE_PARA

! Type to define general properties of finite element
    type SSH_ELEM_PROP
! - Properties of integration scheme
        type(SSH_ELEM_INTE) :: elemInte
! - Type of cell
        integer(kind=8)             :: cellType = SSH_CELL_UNDEF
! - Number of total DOF for the finite element
        integer(kind=8)             :: nbDof = 0
! - Number of displacement DOF for the finite element
        integer(kind=8)             :: nbDofGeom = 0
! - Number of total nodes for the finite element
        integer(kind=8)             :: nbNode = 0
! - Number of displacement nodes for the finite element
        integer(kind=8)             :: nbNodeGeom = 0
    end type SSH_ELEM_PROP

! Type to define geometric properties of cell
    type SSH_CELL_GEOM
! - Center of cell in covariant (parametric) frame
        real(kind=8)        :: cellCenterCova(3) = 0.d0
! - JEVEUX adress to initial coordinates of cell
        integer(kind=8)             :: jvGeom = 0
! - Initial coordinates of finite element
        real(kind=8)        :: geomInit(SSH_NBDOFG_MAX) = 0.d0
        real(kind=8)        :: geomInitX(SSH_NBNODEG_MAX) = 0.d0
        real(kind=8)        :: geomInitY(SSH_NBNODEG_MAX) = 0.d0
        real(kind=8)        :: geomInitZ(SSH_NBNODEG_MAX) = 0.d0
! - Jacobian matrix at center of cell on initial configuration
! - dX/dXi | dX/dEta | dX/dZeta
! - dY/dXi | dY/dEta | dY/dZeta
! - dZ/dXi | dZ/dEta | dZ/dZeta
        real(kind=8)        :: Jac0(3, 3) = 0.d0
! - Inverse of jacobian matrix at center of cell on initial configuration
! - dXi/dX | dEta/dX | dZeta/dX
! - dXi/dY | dEta/dY | dZeta/dY
! - dXi/dZ | dEta/dZ | dZeta/dZ
        real(kind=8)        :: JacInv0(3, 3) = 0.d0
! - Determinant of jacobian matrix at center of cell on initial configuration
        real(kind=8)        :: detJac0 = 0.d0
! - Thickness value
        real(kind=8)        :: h0 = 0.d0
    end type SSH_CELL_GEOM

! ==================================================================================================
!
! Derivated types - Geometry for HEXA
!
! ==================================================================================================

! Type to define HEXA
    type SSH_GEOM_HEXA
! - Base properties of cell
        type(SSH_CELL_GEOM) :: cellGeom
! - Reference configuration
        real(kind=8) :: geomCurr(SSH_NBDOFG_HEXA) = 0.d0
! - T matrix relating the covariant and cartesian frames
        real(kind=8) :: T0(6, 6) = 0.d0
        real(kind=8) :: TXI(6, 6) = 0.d0
        real(kind=8) :: TETA(6, 6) = 0.d0
        real(kind=8) :: TZETA(6, 6) = 0.d0
    end type SSH_GEOM_HEXA

! ==================================================================================================
!
! Derivated types - Kinematic for HEXA
!
! ==================================================================================================

! Green-Lagrange strains
    type SSH_EPSG_HEXA
! - Decomposition of Green-Lagrange strains
        real(kind=8) :: ECova0(6) = 0.d0
        real(kind=8) :: ECovaXI(6) = 0.d0
        real(kind=8) :: ECovaETA(6) = 0.d0
        real(kind=8) :: ECovaZETA(6) = 0.d0
        real(kind=8) :: ECovaETAZETA(6) = 0.d0
        real(kind=8) :: ECovaXIZETA(6) = 0.d0
        real(kind=8) :: ECovaZETAZETA(6) = 0.d0
! - Green-Lagrange strains
        real(kind=8) :: vale(SSH_SIZE_TENS) = 0.d0
    end type SSH_EPSG_HEXA

! Logarithmic strains
    type SSH_EPSL_HEXA
! - Eigen decomposition of C = 2E+I strain mesure
        real(kind=8) :: eigenVale(3) = 0.d0
        real(kind=8) :: eigenVect(3, 3) = 0.d0
        real(kind=8) :: logl(3) = 0.d0
! - Logarithmic strains
        real(kind=8) :: vale(SSH_SIZE_TENS) = 0.d0
    end type SSH_EPSL_HEXA

! Kinematic for HEXA cell
    type SSH_KINE_HEXA
! - Flag if large hypothesis for kinematic
        aster_logical :: lLarge = ASTER_FALSE
! - Parts of gradient matrix (covariant/parametric frame)
        real(kind=8) :: BCova0(6, 24) = 0.d0
        real(kind=8) :: BCovaZETA(6, 24) = 0.d0
        real(kind=8) :: BCovaZETAZETA(6, 24) = 0.d0
        real(kind=8) :: BCovaXI(6, 24) = 0.d0
        real(kind=8) :: BCovaETA(6, 24) = 0.d0
        real(kind=8) :: BCovaETAZETA(6, 24) = 0.d0
        real(kind=8) :: BCovaXIZETA(6, 24) = 0.d0
! - Parts of gradient matrix (cartesian frame)
        real(kind=8) :: BCart0(6, 24) = 0.d0
        real(kind=8) :: BCartZ(6, 24) = 0.d0
        real(kind=8) :: BCartZZ(6, 24) = 0.d0
        real(kind=8) :: BCartX(6, 24) = 0.d0
        real(kind=8) :: BCartY(6, 24) = 0.d0
        real(kind=8) :: BCartXZ(6, 24) = 0.d0
        real(kind=8) :: BCartYZ(6, 24) = 0.d0
! - Green-Lagrange strain at beginning of time step
        type(SSH_EPSG_HEXA) :: epsgPrev
! - Green-Lagrange strain at end of time step
        type(SSH_EPSG_HEXA) :: epsgCurr
! - Logarithmic strains at beginning of time step
        type(SSH_EPSL_HEXA) :: epslPrev
! - Decomposition of Green-Lagrange strain at end of time step
        type(SSH_EPSL_HEXA) :: epslCurr
! - B EAS matrix
        real(kind=8) :: BCovaEAS(6) = 0.d0
        real(kind=8) :: BCartEAS(6) = 0.d0
! - B gradient matrix
        real(kind=8) :: B(6, 25) = 0.d0
    end type SSH_KINE_HEXA

! ==================================================================================================
!
! Derivated types - Stabilization for HEXA
!
! ==================================================================================================

    type SSH_STAB_HEXA
        real(kind=8) :: SXI(24, 24) = 0.d0
        real(kind=8) :: SETA(24, 24) = 0.d0
        real(kind=8) :: SETAZETA(24, 24) = 0.d0
        real(kind=8) :: SXIZETA(24, 24) = 0.d0
        real(kind=8) :: sigmStabXI(SSH_SIZE_TENS) = 0.d0
        real(kind=8) :: sigmStabETA(SSH_SIZE_TENS) = 0.d0
        real(kind=8) :: sigmStabETAZETA(SSH_SIZE_TENS) = 0.d0
        real(kind=8) :: sigmStabXIZETA(SSH_SIZE_TENS) = 0.d0
        real(kind=8) :: matrStabMate(24, 24) = 0.d0
        real(kind=8) :: matrStabGeom(24, 24) = 0.d0
        real(kind=8) :: forcStab(24) = 0.d0
    end type SSH_STAB_HEXA

! ==================================================================================================
!
! Derivated types - Non-linear
!
! ==================================================================================================

    type SSH_BEHA_PARA
! - JEVEUX adress to behaviour parameters
        character(len=16), pointer :: compor(:) => null()
        real(kind=8), pointer :: carcri(:) => null()
! - Flag when compute tangent matrix
        aster_logical         :: lMatr = ASTER_FALSE
        aster_logical         :: lMatrSyme = ASTER_TRUE
! - Flag when compute internal forces vector
        aster_logical         :: lVect = ASTER_FALSE
! - Flag when compute internal state variables
        aster_logical         :: lVari = ASTER_FALSE
! - Flag when computes stresses and returns code error
        aster_logical         :: lSigm = ASTER_FALSE
! - Flag when large strain model
        aster_logical         :: lLarge = ASTER_FALSE
! - Type of non-linear relation
        character(len=16)     :: relaComp = ' '
! - Type of strain model
        character(len=16)     :: defoComp = ' '
! - Type of integration scheme (COMP_ELAS/COMP_INCR)
        character(len=16)     :: typeComp = ' '
! - Main behaviour datastructure
        type(Behaviour_Integ) :: BEHInteg
    end type SSH_BEHA_PARA

!
end module SolidShell_type
