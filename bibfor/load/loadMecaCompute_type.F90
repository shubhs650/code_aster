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
module loadMecaCompute_type
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
! ==================================================================================================
!
! Global variables - General
!
! Definition of parameters for all AFFE_CHAR_MECA loads
!
! ==================================================================================================

! - Keyword to define load in AFFE_CHAR_MECA
    character(len=24), parameter :: mecaLoadKeyword(LOAD_NEUM_NBTYPE) = (/ &
                                    'FORCE_NODALE            ', 'FORCE_INTERNE#3D        ', &
                                    'FORCE_FACE              ', 'FORCE_ARETE             ', &
                                    'FORCE_INTERNE#2D        ', 'FORCE_CONTOUR           ', &
                                    'FORCE_POUTRE            ', 'PESANTEUR               ', &
                                    'ROTATION                ', 'PRES_REP                ', &
                                    'FORCE_COQUE#3D          ', &
                                    'FORCE_COQUE#2D          ', 'PRE_EPSI                ', &
                                    'FLUX_THM_REP            ', 'VECT_ASSE               ', &
                                    'PRE_SIGM                ', 'EFFE_FOND               ', &
                                    'ECHA_THM                ', 'ECHA_HR                 ', &
                                    'VITE_FACE               '/)

! - Name of field to define load in AFFE_CHAR_MECA
    character(len=6), parameter :: mecaLoadField(LOAD_NEUM_NBTYPE) = (/ &
                                   '.FORNO', '.F3D3D', &
                                   '.F2D3D', '.F1D3D', &
                                   '.F2D2D', '.F1D2D', &
                                   '.F1D1D', '.PESAN', &
                                   '.ROTAT', '.PRESS', &
                                   '.FCO3D', &
                                   '.FCO2D', '.EPSIN', &
                                   '.FLUX ', '.VEASS', &
                                   '.SIINT', '.EFOND', &
                                   '.ETHM ', '.ETHMH', &
                                   '.VFACE'/)

! - Flag to authorize that load can been "undead' (depending on unknowns)
    aster_logical, parameter :: mecaLoadSuiv(LOAD_NEUM_NBTYPE) = (/ &
                                                                 .false._1, .false._1, &
                                                                 .true._1, .false._1, &
                                                                 .false._1, .false._1, &
                                                                 .true._1, .true._1, &
                                                                 .true._1, .true._1, &
                                                                 .true._1, &
                                                                 .false._1, .false._1, &
                                                                 .false._1, .false._1, &
                                                                 .false._1, .true._1, &
                                                                 .true._1, .true._1, &
                                                                 .false._1/)

! - Flag to authorize that load can been used in continuation methods (PILOTAGE)
    aster_logical, parameter :: mecaLoadPilo(LOAD_NEUM_NBTYPE) = (/ &
                                                                 .true._1, .true._1, &
                                                                 .true._1, .true._1, &
                                                                 .true._1, .true._1, &
                                                                 .true._1, .true._1, &
                                                                 .false._1, .true._1, &
                                                                 .true._1, &
                                                                 .true._1, .false._1, &
                                                                 .true._1, .true._1, &
                                                                 .false._1, .true._1, &
                                                                 .false._1, .false._1, &
                                                                 .false._1/)

! - Name of input parameter field (real coefficient)
    character(len=8), parameter :: mecaLoadParaR(LOAD_NEUM_NBTYPE) = (/ &
                                   'PFORNOR ', 'PFR3D3D ', &
                                   'PFR2D3D ', 'PFR1D3D ', &
                                   'PFR2D2D ', 'PFR1D2D ', &
                                   'PFR1D1D ', 'PPESANR ', &
                                   'PROTATR ', 'PPRESSR ', &
                                   'PFRCO3D ', &
                                   'PFRCO2D ', 'PEPSINR ', &
                                   'PFLUXR  ', 'NoInput ', &
                                   'PSIEFR  ', 'PEFOND  ', &
                                   'PECHTHM ', 'HECHTHM ', &
                                   'PVITEFR '/)

! - Name of input parameter field (function coefficient)
    character(len=8), parameter :: mecaLoadParaF(LOAD_NEUM_NBTYPE) = (/ &
                                   'PFORNOF ', 'PFF3D3D ', &
                                   'PFF2D3D ', 'PFF1D3D ', &
                                   'PFF2D2D ', 'PFF1D2D ', &
                                   'PFF1D1D ', 'PPESANR ', &
                                   'PROTATR ', 'PPRESSF ', &
                                   'PFFCO3D ', &
                                   'PFFCO2D ', 'PEPSINF ', &
                                   'PFLUXF  ', 'NoInput ', &
                                   'NoInput ', 'PEFOND  ', &
                                   'PCHTHMF ', 'HCHTHMF ', &
                                   'PVITEFF '/)

! - Name of option for dead load (real coefficient)
    character(len=16), parameter :: mecaLoadVectR(LOAD_NEUM_NBTYPE) = (/ &
                                    'CHAR_MECA_FORC_R', 'CHAR_MECA_FR3D3D', &
                                    'CHAR_MECA_FR2D3D', 'CHAR_MECA_FR1D3D', &
                                    'CHAR_MECA_FR2D2D', 'CHAR_MECA_FR1D2D', &
                                    'CHAR_MECA_FR1D1D', 'CHAR_MECA_PESA_R', &
                                    'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRES_R', &
                                    'CHAR_MECA_FRCO3D', &
                                    'CHAR_MECA_FRCO2D', 'CHAR_MECA_EPSI_R', &
                                    'CHAR_MECA_FLUX_R', 'Copy_Load       ', &
                                    'FORC_NODA       ', 'CHAR_MECA_EFON_R', &
                                    'CHAR_ECHA_THM_R ', 'CHAR_ECHA_HR_R  ', &
                                    'CHAR_MECA_VFAC_R'/)

! - Name of option for undead load (real coefficient)
    character(len=16), parameter :: mecaLoadVectSR(LOAD_NEUM_NBTYPE) = (/ &
                                    'NoVector        ', 'NoVector        ', &
                                    'CHAR_MECA_FRSU23', 'NoVector        ', &
                                    'NoVector        ', 'NoVector        ', &
                                    'CHAR_MECA_SR1D1D', 'CHAR_MECA_PESA_R', &
                                    'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRSU_R', &
                                    'CHAR_MECA_SRCO3D', &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        ', 'CHAR_MECA_EFSU_R', &
                                    'CHAR_ECHA_THM_R ', 'CHAR_ECHA_HR_R  ', &
                                    'NoVector        '/)

! - Name of option for dead load (function)
    character(len=16), parameter :: mecaLoadVectF(LOAD_NEUM_NBTYPE) = (/ &
                                    'CHAR_MECA_FORC_F', 'CHAR_MECA_FF3D3D', &
                                    'CHAR_MECA_FF2D3D', 'CHAR_MECA_FF1D3D', &
                                    'CHAR_MECA_FF2D2D', 'CHAR_MECA_FF1D2D', &
                                    'CHAR_MECA_FF1D1D', 'CHAR_MECA_PESA_R', &
                                    'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRES_F', &
                                    'CHAR_MECA_FFCO3D', &
                                    'CHAR_MECA_FFCO2D', 'CHAR_MECA_EPSI_F', &
                                    'CHAR_MECA_FLUX_F', 'Copy_Load       ', &
                                    'FORC_NODA       ', 'CHAR_MECA_EFON_F', &
                                    'CHAR_ECHA_THM_F ', 'CHAR_ECHA_HR_F  ', &
                                    'CHAR_MECA_VFAC_F'/)

! - Name of option for undead load (function)
    character(len=16), parameter :: mecaLoadVectSF(LOAD_NEUM_NBTYPE) = (/ &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        ', 'NoVector        ', &
                                    'CHAR_MECA_SF1D1D', 'NoVector        ', &
                                    'NoVector        ', 'CHAR_MECA_PRSU_F', &
                                    'CHAR_MECA_SFCO3D', &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        ', 'CHAR_MECA_EFSU_F', &
                                    'CHAR_ECHA_THM_F ', 'CHAR_ECHA_HR_F  ', &
                                    'NoVector        '/)

! - Name of option for undead load (real coefficient)
    character(len=16), parameter :: mecaLoadMatrR(LOAD_NEUM_NBTYPE) = (/ &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    'RIGI_MECA_RO    ', 'RIGI_MECA_PRSU_R', &
                                    'RIGI_MECA_SRCO3D', &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", 'RIGI_MECA_EFSU_R', &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        "/)

! - Name of option for undead load (real function)
    character(len=16), parameter :: mecaLoadMatrF(LOAD_NEUM_NBTYPE) = (/ &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", 'RIGI_MECA_PRSU_F', &
                                    'RIGI_MECA_SFCO3D', &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        ", 'RIGI_MECA_EFSU_F', &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    "NoMatrix        "/)

! - Type of matrix for undead load
    character(len=8), parameter :: mecaLoadParaM(LOAD_NEUM_NBTYPE) = (/ &
                                   "NoMatrix", "NoMatrix", &
                                   "NoMatrix", "NoMatrix", &
                                   "NoMatrix", "NoMatrix", &
                                   "NoMatrix", "NoMatrix", &
                                   'PMATUUR ', 'PMATUNS ', &
                                   'PMATUNS ', &
                                   "NoMatrix", "NoMatrix", &
                                   "NoMatrix", "NoMatrix", &
                                   "NoMatrix", 'PMATUNS ', &
                                   "NoMatrix", "NoMatrix", &
                                   "NoMatrix"/)
!
end module loadMecaCompute_type
