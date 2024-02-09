/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#pragma once

#include "astercxx.h"

/* keep consistency with bibfor/include/asterfort/Behaviour_type.h */
// modelling
enum AsterHypothesis {
    HYP_NOT_SET = 0,
    TRIDIMENSIONAL = 1,
    AXISYMMETRICAL = 2,
    PLANESTRESS = 3,
    PLANESTRAIN = 4,
};

// strain model
enum AsterStrainModel { MOD_NOT_SET = 0, SMALL = 1, SIMOMIEHE = 2, GREENLAGRANGE = 3 };

extern "C" {

/* ***** Define fortran interfaces ***** */

/* External State Variables */
void DEFSP( MGIS_GET_NUMBER_OF_ESVS, mgis_get_number_of_esvs, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *nbvar );

void DEFSS( MGIS_GET_ESVS, mgis_get_esvs, const char *hexid, STRING_SIZE l_id, char *names,
            STRING_SIZE l_nam );

void DEFSPPP( MGIS_SET_EXTERNAL_STATE_VARIABLES, mgis_set_external_state_variables,
              const char *hexid, STRING_SIZE l_id, ASTERINTEGER *stid, ASTERDOUBLE *value,
              ASTERINTEGER *nbvar );

/* Internal State Variables */
void DEFSP( MGIS_GET_NUMBER_OF_ISVS, mgis_get_number_of_isvs, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *nbvar );

void DEFSP( MGIS_GET_SIZEOF_ISVS, mgis_get_sizeof_isvs, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *vectsize );

void DEFSS( MGIS_GET_ISVS, mgis_get_isvs, const char *hexid, STRING_SIZE l_id, char *names,
            STRING_SIZE l_nam );

void DEFSP( MGIS_GET_ISVS_SIZES, mgis_get_isvs_sizes, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *sizes );

void DEFSPPP( MGIS_SET_INTERNAL_STATE_VARIABLES, mgis_set_internal_state_variables,
              const char *hexid, STRING_SIZE l_id, ASTERINTEGER *stid, ASTERDOUBLE *value,
              ASTERINTEGER *vectsize );

/* Other State Variables */
void DEFSPPP( MGIS_SET_GRADIENTS, mgis_set_gradients, const char *hexid, STRING_SIZE l_id,
              ASTERINTEGER *stid, ASTERDOUBLE *values, ASTERINTEGER *nbval );

void DEFSPPP( MGIS_SET_THERMODYNAMIC_FORCES, mgis_set_thermodynamic_forces, const char *hexid,
              STRING_SIZE l_id, ASTERINTEGER *stid, ASTERDOUBLE *values, ASTERINTEGER *nbval );

void DEFSP( MGIS_SET_ROTATION_MATRIX, mgis_set_rotation_matrix, const char *hexid, STRING_SIZE l_id,
            ASTERDOUBLE *values );

/* Parameters */
void DEFSSP( MGIS_GET_DOUBLE_MFRONT_PARAMETER, mgis_get_double_mfront_parameter, const char *hexid,
             STRING_SIZE l_id, const char *param_, STRING_SIZE l_par, ASTERDOUBLE *value );

void DEFSSP( MGIS_SET_DOUBLE_PARAMETER, mgis_set_double_parameter, const char *hexid,
             STRING_SIZE l_id, const char *param_, STRING_SIZE l_par, ASTERDOUBLE *value );

void DEFSSP( MGIS_GET_INTEGER_MFRONT_PARAMETER, mgis_get_integer_mfront_parameter,
             const char *hexid, STRING_SIZE l_id, const char *param_, STRING_SIZE l_par,
             ASTERINTEGER *value );

void DEFSSP( MGIS_SET_INTEGER_PARAMETER, mgis_set_integer_parameter, const char *hexid,
             STRING_SIZE l_id, const char *param_, STRING_SIZE l_par, ASTERINTEGER *value );

void DEFSP( MGIS_SET_OUTOFBOUNDS_POLICY, mgis_set_outofbounds_policy, const char *hexid,
            STRING_SIZE l_id, ASTERINTEGER *value );

/* Material Properties */
void DEFSP( MGIS_GET_NUMBER_OF_PROPS, mgis_get_number_of_props, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *nbprop );

void DEFSS( MGIS_GET_PROPS, mgis_get_props, const char *hexid, STRING_SIZE l_id, char *names,
            STRING_SIZE l_nam );

void DEFSPPP( MGIS_SET_MATERIAL_PROPERTIES, mgis_set_material_properties, const char *hexid,
              STRING_SIZE l_id, ASTERINTEGER *stid, ASTERDOUBLE *values, ASTERINTEGER *nbprop );

/* ***** General ***** */
void DEFSPP( MGIS_LOAD_LIBRARY, mgis_load_library, const char *hexid, STRING_SIZE l_id,
             ASTERINTEGER *model_, ASTERINTEGER *strain_ );

void DEFSPPPPPP( MGIS_INTEGRATE, mgis_integrate, const char *hexid, STRING_SIZE l_id,
                 ASTERDOUBLE *stress, ASTERDOUBLE *statev, ASTERDOUBLE *ddsdde, ASTERDOUBLE *dtime,
                 ASTERDOUBLE *pnewdt, ASTERINTEGER *retcode );

void DEFSS( MGIS_DEBUG, mgis_debug, const char *hexid, STRING_SIZE l_id, const char *title,
            STRING_SIZE l_tit );
}
