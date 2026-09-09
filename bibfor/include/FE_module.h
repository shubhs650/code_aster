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
! FE Size module : Parameters <-> integer definitions
! -------------------------------------------------------------------------
!
#include "MeshTypes_type.h"
!
! - Static size - FE methods - General
!
! --- maximum number of basis function
#define MAX_BS_CG MT_NNOMAX3D
#define MAX_BV_CG 3*MT_NNOMAX3D
!
! --- EF Lagrange
#define EF_LAGRANGE 0
!
! --- maximum number of quadrature points
#define MAX_QP 64
! --- maximum number of quadrature points on a face QUAD = 16
#define MAX_QP_FACE 16
! --- maximum number of quadrature points on a cell HEXA = 64
#define MAX_QP_CELL 64
