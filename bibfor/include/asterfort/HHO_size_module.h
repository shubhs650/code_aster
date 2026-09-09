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
! HHO Size module : Parameters <-> integer definitions
! -------------------------------------------------------------------------
!
! - Static size - HHO methods - General
!
! --- maximum degree of a face
#define MAX_DEGREE_FACE 4
! --- maximum number of a cell
#define MAX_DEGREE_CELL 5

! --- maximum number of face for a cell (6 for a hexahedron)
#define MAX_FACE 6
!
! --- Maximum number of componants for scalar cell function (order=5 and dim=3)
#define MSIZE_CELL_SCAL 56
! --- Maximum number of componants for scalar face function (order=4 and dim=3)
#define MSIZE_FACE_SCAL 15
! --- Maximum number of total dofs for a HHO function (dim=3,  and order=4 for a face)
! --- Max size for a hexahedron:   6 (nb faces) * 15 (face dofs)
#define MSIZE_FDOFS_SCAL 90
! --- Maximum number of total dofs for a HHO function (dim=3, order=4 for a cell
! --- and order=4 for a face)
! --- Max size for a hexahedron:  35 (cell dofs) + 6 (nb faces) * 15 (face dofs)
#define MSIZE_TDOFS_SCAL 125
!
! --- vector function
! --- Maximum number of componants for vector cell function (order=5 and dim=3)
#define MSIZE_CELL_VEC 168
! --- Maximum number of componants for vector face function (order=4 and dim=3)
#define MSIZE_FACE_VEC 45
! --- Maximum number of total dofs for vector faces (dim=3,  and order=4 for a face)
! --- Max size for a hexahedron:   6 (nb faces) * 45 (face dofs)
#define MSIZE_FDOFS_VEC 270
! --- Maximum number of total dofs for a HHO function (dim=3, order=4 for a cell
! --- and order=4 for a face)
! --- Max size for a hexahedron:  105 (cell dofs) + 6 (nb faces) * 45 (face dofs)
#define MSIZE_TDOFS_VEC 375

! --- Maximum number of total dofs for a HHO function (dim=3, order=2 for a cell
! --- and order=2 for a face)
! --- DEPL : Max size for a hexahedron:  30 (cell dofs) + 6 (nb faces) * 18 (face dofs)
! --- VARI : Max size for a hexahedron:  10 (cell dofs) + 6 (nb faces) * 6 (face dofs)
! --- LAGR : Max size for a hexahedron:  10 (cell dofs)
#define MSIZE_TDOFS_MIX 194

!
! --- matrix function
! --- Maximum number of componants for matrix cell function (order=5 and dim=3)
#define MSIZE_CELL_MAT 315
!
! --- maximum number of quadrature points
#define MSIZE_QP 370
! --- maximum number of quadrature points on a face QUAD = 16
#define MSIZE_QP_FACE 60
! --- maximum number of quadrature points on a cell HEXA = 64
#define MSIZE_QP_CELL 370
