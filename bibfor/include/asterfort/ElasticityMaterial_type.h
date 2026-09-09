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
! Types for elasticity material
!
! ==================================================================================================
!
!
! --------------------------------------------------------------------------------------------------
! Type of elasticity
! --------------------------------------------------------------------------------------------------
#define ELAS_UNDEF     0
#define ELAS_ISOT      1
#define ELAS_ORTH      2
#define ELAS_ISTR      3
#define ELAS_VISC_ISOT 4
#define ELAS_VISC_ORTH 5
#define ELAS_VISC_ISTR 6
#define ELAS_SHELL     7
#define ELAS_GLRC      8
#define ELAS_DHRC      9
#define ELAS_MEMBRANE  10
#define ELAS_COMPOSITE 11
#define ELAS_FLUID     12

! --------------------------------------------------------------------------------------------------
! External state variables linked to elasticity
! --------------------------------------------------------------------------------------------------
#define ELAS_VARC_HYDR 1
#define ELAS_VARC_SECH 1
#define ELAS_VARC_PTOT 1
#define ELAS_VARC_EPSA 1
#define ELAS_VARC_THER 1
