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
! Contact module: Parameters <-> integer definitions
! -------------------------------------------------------------------------
!
#define MAX_LAGA_DOFS 66

#define MAX_PENA_DOFS 66

! value for 2D problem
#define MAX_NITS_DOFS 50

#define TOLE_BORNE 0.d0

#define PROJ_TOLE 1.d-8

#define THRES_STICK 10.d300

! See ContactEnum.h
! Reso Geom 
#define RESO_GEOM_NEWT 0
#define RESO_GEOM_FIXP 1

! Contact
#define CONT_ALGO_LAGR 0
#define CONT_ALGO_NITS 1
#define CONT_ALGO_PENA 2

#define CONT_VARI_NONE 0
#define CONT_VARI_RAPI 1
#define CONT_VARI_ROBU 2
#define CONT_VARI_SYME 3
#define CONT_VARI_CLAS 4

#define CONT_TYPE_UNIL 0
#define CONT_TYPE_BILA 1

! Friction
#define FRIC_ALGO_LAGR 0
#define FRIC_ALGO_NITS 1
#define FRIC_ALGO_PENA 2
#define FRIC_ALGO_LAGR_STD 3

#define FRIC_TYPE_NONE 0
#define FRIC_TYPE_TRES 1
#define FRIC_TYPE_COUL 2
#define FRIC_TYPE_STIC 3

! Pairing
#define PAIR_CONT_INTE 0
#define FRIC_ALGO_FALSE 1
#define FRIC_ALGO_TRUE 2
