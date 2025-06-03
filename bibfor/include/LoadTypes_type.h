! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

! Maximum number of maps for a load
#define LOAD_MAP_NBMAX     2
! Maximum number of components in a map for a load
#define LOAD_MAP_NBCMPMAX  10

! Maximum number type in one load (for list of loads)
#define LOAD_NBIDEN_MAXI 99

! Maximum number of input fields to compute Neumann loads (mechanic)
#define LOAD_NEUM_NBMAXIN 48

! Maximum number of Neumann loads type (mechanic)
#define LOAD_NEUM_NBTYPE 21

! Identified types for Neumann loads (mechanic)
#define LOAD_NEUM_UNKNOWN -1
#define LOAD_NEUM_FORC_NODA 1
#define LOAD_NEUM_INTE_3D  2
#define LOAD_NEUM_FORC_FACE 3
#define LOAD_NEUM_FORC_EDGE 4
#define LOAD_NEUM_INTE_2D 5
#define LOAD_NEUM_FORC_WIRE 6
#define LOAD_NEUM_FORC_BEAM 7
#define LOAD_NEUM_GRAVITY 8
#define LOAD_NEUM_ROTATION 9
#define LOAD_NEUM_PRESSURE 10
#define LOAD_NEUM_FORC_ELEC 11
#define LOAD_NEUM_FORC_SHEL_3D 12
#define LOAD_NEUM_FORC_SHEL_2D 13
#define LOAD_NEUM_PRE_EPSI 14
#define LOAD_NEUM_FLUX 15
#define LOAD_NEUM_VECT_ASSE 16
#define LOAD_NEUM_PRE_SIGM 17
#define LOAD_NEUM_EFFE_FOND 18
#define LOAD_NEUM_THM_ECHA 19
#define LOAD_NEUM_THM_ECHAH 20

#define LOAD_NEUM_PWAVE 100
#define LOAD_NEUM_VITE_FACE 101
#define LOAD_NEUM_GROUND 102
#define LOAD_NEUM_WAVE 103

! Maximum number of input fields to compute Neumann loads (thermic)
#define LOAD_NEUT_NBMAXIN 16

! Maximum number of Neumann loads type (thermic)
#define LOAD_NEUT_NBTYPE 11

! Identified types for Neumann loads (thermic)
#define LOAD_NEUT_UNKNOWN -1
#define LOAD_NEUT_ECHANGE 1
#define LOAD_NEUT_FLUX_XYZ 2
#define LOAD_NEUT_FLUX_NORM 3
#define LOAD_NEUT_SOURCE 4
#define LOAD_NEUT_ECH_PAROI 5
#define LOAD_NEUT_PRE_GRAD 6
#define LOAD_NEUT_FLUX_NL 7
#define LOAD_NEUT_SOUR_NL 8
#define LOAD_NEUT_RAYO 9
#define LOAD_NEUT_EVOL_CHAR 10
#define LOAD_NEUT_SOURCE_CALC 11

! Maximum number of Neumann loads type (acoustic)
#define LOAD_NEUA_NBTYPE 1

! Identified types for Neumann loads (acoustic)
#define LOAD_NEUA_UNKNOWN -1
#define LOAD_NEUA_VITE_FACE 1
