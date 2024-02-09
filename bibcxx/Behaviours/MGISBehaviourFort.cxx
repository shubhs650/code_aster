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

#include "MGISBehaviourFort.h"

#include "MGISBehaviour.h"

#ifdef ASTER_HAVE_MGIS
#include "MGIS/Raise.hxx"
#endif
#include "Utilities/Tools.h"

/* ***** Define fortran interfaces ***** */

#ifdef ASTER_HAVE_MGIS
inline MGISBehaviour *getPtr( const std::string &hexid ) {
    AS_ASSERT( strip( hexid ) != "" && strip( hexid ) != "VIDE" );
    ASTERINTEGER addr = (ASTERINTEGER)std::stol( hexid, nullptr, 16 );
    return (MGISBehaviour *)( addr );
}

inline MGISBehaviour *getPtr( const char *hexid, STRING_SIZE l_id ) {
    std::string str_id( hexid, l_id );
    return getPtr( str_id );
}
#endif

/* External State VariableS */

// call mgis_get_number_of_esvs(hexid, nbvar)
void DEFSP( MGIS_GET_NUMBER_OF_ESVS, mgis_get_number_of_esvs, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *nbvar ) {
#ifdef ASTER_HAVE_MGIS
    *nbvar = ( ASTERINTEGER )( getPtr( hexid, l_id )->getNumberOfExternalStateVariables() );
#endif
}

// call mgis_get_esvs(hexid, names)
void DEFSS( MGIS_GET_ESVS, mgis_get_esvs, const char *hexid, STRING_SIZE l_id, char *names,
            STRING_SIZE l_nam ) {
#ifdef ASTER_HAVE_MGIS
    // output array must be big enough
    const auto varnames = getPtr( hexid, l_id )->getExternalStateVariablesNames();
    vectorStringToFStrArray( names, l_nam, varnames );
#endif
}

void DEFSPPP( MGIS_SET_EXTERNAL_STATE_VARIABLES, mgis_set_external_state_variables,
              const char *hexid, STRING_SIZE l_id, ASTERINTEGER *stid, ASTERDOUBLE *values,
              ASTERINTEGER *nbvar ) {
#ifdef ASTER_HAVE_MGIS
    if ( *nbvar < 1 ) {
        return;
    }
    AS_ASSERT( (int)*nbvar == getPtr( hexid, l_id )->getNumberOfExternalStateVariables() );
    getPtr( hexid, l_id )->setExternalStateVariables( (StateId)*stid, values );
#endif
}

/* Internal State VariableS */

void DEFSP( MGIS_GET_NUMBER_OF_ISVS, mgis_get_number_of_isvs, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *nbvar ) {
#ifdef ASTER_HAVE_MGIS
    *nbvar = ( ASTERINTEGER )( getPtr( hexid, l_id )->getNumberOfInternalStateVariables() );
#endif
}

// call mgis_get_sizeof_isvs( hexid, vectsize )
void DEFSP( MGIS_GET_SIZEOF_ISVS, mgis_get_sizeof_isvs, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *vectsize ) {
#ifdef ASTER_HAVE_MGIS
    *vectsize = ( ASTERINTEGER )( getPtr( hexid, l_id )->getSizeOfInternalStateVariables() );
#endif
}

// call mgis_get_isvs(hexid, names)
void DEFSS( MGIS_GET_ISVS, mgis_get_isvs, const char *hexid, STRING_SIZE l_id, char *names,
            STRING_SIZE l_nam ) {
#ifdef ASTER_HAVE_MGIS
    // output array must be big enough
    const auto varnames = getPtr( hexid, l_id )->getInternalStateVariablesNames();
    vectorStringToFStrArray( names, l_nam, varnames );
#endif
}

void DEFSP( MGIS_GET_ISVS_SIZES, mgis_get_isvs_sizes, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *sizes ) {
#ifdef ASTER_HAVE_MGIS
    // output array must be big enough
    auto varsizes = getPtr( hexid, l_id )->getInternalStateVariablesSizes();
    for ( auto i = 0; i < varsizes.size(); ++i ) {
        sizes[i] = (ASTERINTEGER)varsizes[i];
    }
#endif
}

void DEFSPPP( MGIS_SET_INTERNAL_STATE_VARIABLES, mgis_set_internal_state_variables,
              const char *hexid, STRING_SIZE l_id, ASTERINTEGER *stid, ASTERDOUBLE *values,
              ASTERINTEGER *nbvar ) {
#ifdef ASTER_HAVE_MGIS
    auto size = getPtr( hexid, l_id )->getSizeOfInternalStateVariables();
    if ( size == 0 && (int)( *nbvar ) == 1 ) {
        return;
    }
    AS_ASSERT( (int)( *nbvar ) == size );
    getPtr( hexid, l_id )->setInternalStateVariables( (StateId)*stid, values );
#endif
}

/* Other State Variables */

void DEFSPPP( MGIS_SET_GRADIENTS, mgis_set_gradients, const char *hexid, STRING_SIZE l_id,
              ASTERINTEGER *stid, ASTERDOUBLE *values, ASTERINTEGER *nbval ) {
#ifdef ASTER_HAVE_MGIS
    int size = getPtr( hexid, l_id )->getSizeOfGradients();
    if ( (int)*nbval != size ) {
        std::cout << "WARNING: size of deformation gradients: " << size
                  << " != nbval: " << (int)*nbval << std::endl;
        AS_ASSERT( (int)*nbval >= size );
    }
    getPtr( hexid, l_id )->setGradients( (StateId)*stid, values, (int)*nbval );
#endif
}

void DEFSPPP( MGIS_SET_THERMODYNAMIC_FORCES, mgis_set_thermodynamic_forces, const char *hexid,
              STRING_SIZE l_id, ASTERINTEGER *stid, ASTERDOUBLE *values, ASTERINTEGER *nbval ) {
#ifdef ASTER_HAVE_MGIS
    int size = getPtr( hexid, l_id )->getSizeOfThermodynamicForces();
    if ( (int)*nbval != size ) {
        std::cout << "WARNING: size of thermodynamic forces: " << size
                  << " != nbval: " << (int)*nbval << std::endl;
        AS_ASSERT( (int)*nbval >= size );
    }
    AS_ASSERT( (int)*nbval >= getPtr( hexid, l_id )->getSizeOfThermodynamicForces() );
    getPtr( hexid, l_id )->setThermodynamicForces( (StateId)*stid, values, (int)*nbval );
#endif
}

void DEFSP( MGIS_SET_ROTATION_MATRIX, mgis_set_rotation_matrix, const char *hexid, STRING_SIZE l_id,
            ASTERDOUBLE *values ) {
#ifdef ASTER_HAVE_MGIS
    getPtr( hexid, l_id )->setRotationMatrix( values );
#endif
}

/* Parameters */

void DEFSSP( MGIS_GET_DOUBLE_MFRONT_PARAMETER, mgis_get_double_mfront_parameter, const char *hexid,
             STRING_SIZE l_id, const char *param_, STRING_SIZE l_par, ASTERDOUBLE *value ) {
#ifdef ASTER_HAVE_MGIS
    mgis::ExceptionHandler h = nullptr;
    mgis::setExceptionHandler( h );
    std::string param = strip( std::string( param_, l_par ) );
    try {
        *value = getPtr( hexid, l_id )->getMFrontParameter< ASTERDOUBLE, double >( param );
    } catch ( ... ) {
        *value = 0;
    }
    mgis::setExceptionHandler( MGISExceptionHandler );
#endif
}

void DEFSSP( MGIS_SET_DOUBLE_PARAMETER, mgis_set_double_parameter, const char *hexid,
             STRING_SIZE l_id, const char *param_, STRING_SIZE l_par, ASTERDOUBLE *value ) {
#ifdef ASTER_HAVE_MGIS
    std::string param = strip( std::string( param_, l_par ) );
    getPtr( hexid, l_id )->setParameter( param, *value );
#endif
}

void DEFSSP( MGIS_GET_INTEGER_MFRONT_PARAMETER, mgis_get_integer_mfront_parameter,
             const char *hexid, STRING_SIZE l_id, const char *param_, STRING_SIZE l_par,
             ASTERINTEGER *value ) {
#ifdef ASTER_HAVE_MGIS
    mgis::ExceptionHandler h = nullptr;
    mgis::setExceptionHandler( h );
    std::string param = strip( std::string( param_, l_par ) );
    try {
        *value = getPtr( hexid, l_id )->getMFrontParameter< ASTERINTEGER, int >( param );
    } catch ( ... ) {
        *value = 0;
    }
    mgis::setExceptionHandler( MGISExceptionHandler );
#endif
}

void DEFSSP( MGIS_SET_INTEGER_PARAMETER, mgis_set_integer_parameter, const char *hexid,
             STRING_SIZE l_id, const char *param_, STRING_SIZE l_par, ASTERINTEGER *value ) {
#ifdef ASTER_HAVE_MGIS
    std::string param = strip( std::string( param_, l_par ) );
    getPtr( hexid, l_id )->setParameter( param, *value );
#endif
}

void DEFSP( MGIS_SET_OUTOFBOUNDS_POLICY, mgis_set_outofbounds_policy, const char *hexid,
            STRING_SIZE l_id, ASTERINTEGER *value ) {
#ifdef ASTER_HAVE_MGIS
    // getPtr(hexid, l_id)->setParameter( "OutOfBoundsPolicy", (int)*value );
    std::cout << "<I> setting OutOfBoundsPolicy is not supported." << std::endl;
#endif
}

/* Material Properties */

void DEFSP( MGIS_GET_NUMBER_OF_PROPS, mgis_get_number_of_props, const char *hexid, STRING_SIZE l_id,
            ASTERINTEGER *nbprop ) {
#ifdef ASTER_HAVE_MGIS
    *nbprop = ( ASTERINTEGER )( getPtr( hexid, l_id )->getNumberOfMaterialProperties() );
#endif
}

void DEFSS( MGIS_GET_PROPS, mgis_get_props, const char *hexid, STRING_SIZE l_id, char *names,
            STRING_SIZE l_nam ) {
#ifdef ASTER_HAVE_MGIS
    // output array must be big enough
    const auto propnames = getPtr( hexid, l_id )->getMaterialPropertiesNames();
    VectorString as_names;
    for ( auto prop : propnames ) {
        as_names.push_back( remove_brackets( prop ) );
    }
    vectorStringToFStrArray( names, l_nam, as_names );
#endif
}

void DEFSPPP( MGIS_SET_MATERIAL_PROPERTIES, mgis_set_material_properties, const char *hexid,
              STRING_SIZE l_id, ASTERINTEGER *stid, ASTERDOUBLE *values, ASTERINTEGER *nbprop ) {
#ifdef ASTER_HAVE_MGIS
    AS_ASSERT( (int)*nbprop == getPtr( hexid, l_id )->getNumberOfMaterialProperties() );
    getPtr( hexid, l_id )->setMaterialProperties( (StateId)*stid, values );
#endif
}

/* General */

// call mgis_load_library(hexid, model, strain)
void DEFSPP( MGIS_LOAD_LIBRARY, mgis_load_library, const char *hexid, STRING_SIZE l_id,
             ASTERINTEGER *model_, ASTERINTEGER *strain_ ) {
#ifdef ASTER_HAVE_MGIS
    namespace MGB = mgis::behaviour;

    int model = (int)( *model_ );
    int strain = (int)( *strain_ );
#ifdef ASTER_DEBUG_CXX
    std::string str_hexid( hexid, l_id );
    std::cout << "MGISDBG: request MGISBehaviour from id: " << strip( str_hexid )
              << ", model: " << model << ", strain: " << strain << std::endl
              << std::flush;
#endif

    MGB::Hypothesis hyp;
    if ( model == AsterHypothesis::TRIDIMENSIONAL ) {
        hyp = MGB::Hypothesis::TRIDIMENSIONAL;
    } else if ( model == AsterHypothesis::AXISYMMETRICAL ) {
        hyp = MGB::Hypothesis::AXISYMMETRICAL;
    } else if ( model == AsterHypothesis::PLANESTRESS ) {
        hyp = MGB::Hypothesis::PLANESTRESS;
    } else if ( model == AsterHypothesis::PLANESTRAIN ) {
        hyp = MGB::Hypothesis::PLANESTRAIN;
    } else {
        std::cout << "Unexpected value: model : " << model << std::endl;
        AS_ABORT( "mgis_load_library: invalid value for 'model'" );
    }

    if ( strain == AsterStrainModel::SMALL ) {
        getPtr( hexid, l_id )->load( hyp );
    } else {
        MGB::FiniteStrainBehaviourOptions::StressMeasure stress_measure;
        MGB::FiniteStrainBehaviourOptions::TangentOperator tangent_op;
        switch ( strain ) {
        case AsterStrainModel::SMALL:
            AS_ASSERT( false );
        case AsterStrainModel::SIMOMIEHE:
            stress_measure = MGB::FiniteStrainBehaviourOptions::CAUCHY;
            tangent_op = MGB::FiniteStrainBehaviourOptions::DTAU_DDF;
            break;
        case AsterStrainModel::GREENLAGRANGE:
            stress_measure = MGB::FiniteStrainBehaviourOptions::PK2;
            tangent_op = MGB::FiniteStrainBehaviourOptions::DS_DEGL;
            break;
        default:
            std::cout << "Unexpected value: string : " << strain << std::endl;
            AS_ABORT( "mgis_load_library: invalid value for 'strain'" );
        }
        getPtr( hexid, l_id )->load( hyp, stress_measure, tangent_op );
    }
#endif
}

// call mgis_integrate(data...)
void DEFSPPPPPP( MGIS_INTEGRATE, mgis_integrate, const char *hexid, STRING_SIZE l_id,
                 ASTERDOUBLE *stress, ASTERDOUBLE *statev, ASTERDOUBLE *ddsdde, ASTERDOUBLE *dtime,
                 ASTERDOUBLE *pnewdt, ASTERINTEGER *retcode ) {
#ifdef ASTER_HAVE_MGIS
    auto behaviour = getPtr( hexid, l_id );
    // fill BehaviourData object with inputs
    behaviour->setInitialState( *dtime, ddsdde );

    *retcode = (ASTERINTEGER)behaviour->integrate();

    // decode outputs from BehaviourData
    behaviour->setOutputs( stress, statev, ddsdde, pnewdt );
#endif
}

void DEFSS( MGIS_DEBUG, mgis_debug, const char *hexid, STRING_SIZE l_id, const char *title,
            STRING_SIZE l_tit ) {
#ifdef ASTER_HAVE_MGIS
#ifdef ASTER_DEBUG_CXX
    std::cout << std::string( title, l_tit ) << std::endl << std::flush;
    getPtr( hexid, l_id )->repr();
#endif
#endif
}
