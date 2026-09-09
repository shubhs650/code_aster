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
! Indicator select phase of algorithm
#define PRED_EULER          1
#define CORR_NEWTON         2
#define INTE_FORCE          3
#define ACCEL_INIT          4
#define POST_BUCKLING       5

! Set to 1 to activate DEBUG
#define NONLINEAR_DEBUG     0

! Indicator to combine nodal fields for internal forces
#define INTE_FORCE_NONE     0
#define INTE_FORCE_COMB     1
#define INTE_FORCE_INTE     2
#define INTE_FORCE_FNOD     3

! Number of events for management of algorithm
#define ZEVEN               37
#define EVENT_IS_INACTIVE   0
#define EVENT_IS_ACTIVE     1

#define NB_LOOP             5

#define LOOP_RESI           1
#define LOOP_NEWT           2
#define LOOP_FIXE           3
#define LOOP_INST           4
#define LOOP_CALC           5

#define LOOP_STATE_CONTINUE 0
#define LOOP_STATE_CONVERGE 1
#define LOOP_STATE_EVENT    2
#define LOOP_STATE_ERROR    3
#define LOOP_STATE_STOP     4
#define LOOP_STATE_CTCD     5
#define LOOP_STATE_ARCH     6
