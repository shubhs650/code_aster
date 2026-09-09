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
! Coupling module : Parameters <-> integer definitions
! -------------------------------------------------------------------------
!
#include "asterfort/HHO_size_module.h"
#include "FE_module.h"
#include "MeshTypes_type.h"
!
! - Offset for pairing - see CouplingPairing.cxx
! - Add +1 since offset C -> fortran
#define    OFFSET_COEF_PENA  1
#define    OFFSET_NB_PTS_INTER  2
#define    OFFSET_COOR_PTS_X  3
#define    OFFSET_COOR_PTS_Y  11
#define    OFFSET_NB_NODES_NITSCHE  19
#define    OFFSET_NODE_IDX_FACE  20
#define    OFFSET_SIZE  27
!
!
! - Static size - FE methods - General
!
! --- maximum number of basis function
#define MSIZE_CPL_PENA 3*MT_NNOMAX2D+MSIZE_FACE_VEC

! FEM-FEM with QU12/QU12 2*3*9=54
#define MSIZE_CPL_PENA_FEM 2*3*MT_NNOMAX2D

! FEM-FEM with H32/QU12
#define MSIZE_CPL_NITS 3*MT_NNOMAX+MSIZE_FACE_VEC

! FEM-FEM with QU12/QU12 6*9 + 3*9=81
#define MSIZE_CPL_LAGR_FEM 6*MT_NNOMAX2D+3*MT_NNOMAX2D
#define MSIZE_LAGR_FEM 3*MT_NNOMAX2D