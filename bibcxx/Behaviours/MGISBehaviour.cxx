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

#include "astercxx.h"

#include "MGISBehaviour.h"

#ifdef ASTER_HAVE_MGIS

#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/Behaviour/State.hxx"
#include "MGIS/LibrariesManager.hxx"
#include "MGIS/Raise.hxx"
#include "Messages/Messages.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

#include <algorithm>

namespace MGB = mgis::behaviour;

MGB::State *MGISBehaviour::getState( const StateId stid ) const {
    MGB::State *pstat;
    if ( stid == StateId::s0 ) {
        pstat = &( _data->s0 );
    } else {
        pstat = &( _data->s1 );
    }
    return pstat;
}

void _show_strings( std::string title, VectorString &vectstr ) {
    std::cout << title << ": " << std::endl;
    auto i = 0;
    for ( auto param : vectstr ) {
        if ( i > 0 ) {
            std::cout << ", ";
        }
        std::cout << param;
        i++;
    }
    std::cout << std::endl;
}

void MGISBehaviour::repr() const {
    print_markdown( std::cout, *_behav, *_data, 0 );

    _show_strings( "Parameters (real)", _behav->params );
    _show_strings( "Parameters (integer)", _behav->iparams );
    _show_strings( "Parameters (unsigned short)", _behav->usparams );
    std::cout << "nopts: " << _behav->nopts << std::endl;
}

VectorString MGISBehaviour::internal_state() const {
    // add arguments for `load()`?
    VectorString ret = { _libpath, _bname };
    return ret;
}

void MGISBehaviour::_load_library(
    const MGB::Hypothesis hyp, const bool finiteStrain,
    const MGB::FiniteStrainBehaviourOptions::StressMeasure stress_measure,
    const MGB::FiniteStrainBehaviourOptions::TangentOperator tangent_op ) {
    if ( _behav ) {
        return;
    }

    mgis::setExceptionHandler( MGISExceptionHandler );
#ifdef ASTER_DEBUG_CXX
    std::cout << "MGISDBG: loading library " << _libpath << ", for behaviour: " << _bname << "..."
              << std::endl
              << std::flush;
#endif
    if ( !finiteStrain ) {
        _behav = std::make_shared< MGB::Behaviour >( MGB::load( _libpath, _bname, hyp ) );
    } else {
        auto opts = MGB::FiniteStrainBehaviourOptions {};
        opts.stress_measure = stress_measure;
        opts.tangent_operator = tangent_op;

        _behav = std::make_shared< MGB::Behaviour >( MGB::load( opts, _libpath, _bname, hyp ) );
    }
    std::cout << "<I> Library loaded from " << _libpath << ", behaviour '" << _bname
              << "', based on TFEL version " << _behav->tfel_version << std::endl
              << std::flush;
    AS_ASSERT( _behav );
    _data = std::make_shared< MGB::BehaviourData >( MGB::BehaviourData( *_behav ) );
    AS_ASSERT( _data );
#ifdef ASTER_DEBUG_CXX
    std::cout << "Loaded behaviour:" << std::endl << std::flush;
    repr();
#endif
}

void MGISBehaviour::load( const MGB::Hypothesis hyp ) { _load_library( hyp ); }

void MGISBehaviour::load( const MGB::Hypothesis hyp,
                          const MGB::FiniteStrainBehaviourOptions::StressMeasure stress_measure,
                          const MGB::FiniteStrainBehaviourOptions::TangentOperator tangent_op ) {
    _load_library( hyp, true, stress_measure, tangent_op );
}

/* External State VariableS */

int MGISBehaviour::getNumberOfExternalStateVariables() const {
    int size = _behav->esvs.size();
    return size;
}

VectorString MGISBehaviour::getExternalStateVariablesNames() const {
    VectorString names;
    for ( auto &var : _behav->esvs ) {
        auto name = var.name;
        names.push_back( name );
    }
    return names;
}

void MGISBehaviour::setExternalStateVariables( const StateId stid, const ASTERDOUBLE *values ) {
    MGB::State *state = getState( stid );
    for ( auto i = 0; i < getNumberOfExternalStateVariables(); ++i ) {
        state->external_state_variables[i] = (double)values[i];
    }
}

/* Internal State VariableS */

int MGISBehaviour::getNumberOfInternalStateVariables() const { return _behav->isvs.size(); }

int MGISBehaviour::getSizeOfInternalStateVariables() const {
    return getArraySize( _behav->isvs, _behav->hypothesis );
}

VectorString MGISBehaviour::getInternalStateVariablesNames() const {
    VectorString names;
    for ( auto &var : _behav->isvs ) {
        names.push_back( var.name );
    }
    return names;
}

VectorLong MGISBehaviour::getInternalStateVariablesTypes() const {
    VectorLong types;
    for ( auto &var : _behav->isvs ) {
        types.push_back( var.type );
    }
    return types;
}

VectorLong MGISBehaviour::getInternalStateVariablesSizes() const {
    VectorLong sizes;
    for ( auto &var : _behav->isvs ) {
        sizes.push_back( MGB::getVariableSize( var, _behav->hypothesis ) );
    }
    return sizes;
}

void MGISBehaviour::setInternalStateVariables( const StateId stid, const ASTERDOUBLE *values ) {
    MGB::State *state = getState( stid );
    for ( auto i = 0; i < getSizeOfInternalStateVariables(); ++i ) {
        state->internal_state_variables[i] = (double)values[i];
    }
}

/* Other State Variables */

int MGISBehaviour::getSizeOfGradients() const {
    return getArraySize( _behav->gradients, _behav->hypothesis );
}

void MGISBehaviour::setGradients( const StateId stid, const ASTERDOUBLE *values,
                                  const int insize ) {
    MGB::State *state = getState( stid );
    convertTensorToMgis( values, insize, state->gradients, getSizeOfGradients() );
}

int MGISBehaviour::getSizeOfThermodynamicForces() const {
    return getArraySize( _behav->thermodynamic_forces, _behav->hypothesis );
}

void MGISBehaviour::setThermodynamicForces( const StateId stid, const ASTERDOUBLE *values,
                                            const int insize ) {
    MGB::State *state = getState( stid );
    for ( auto i = 0; i < getSizeOfThermodynamicForces(); ++i ) {
        state->thermodynamic_forces[i] = (double)values[i];
    }
}

void MGISBehaviour::setRotationMatrix( const ASTERDOUBLE *drot ) {
    _rotmatrix.clear();
    for ( auto i = 0; i < 9; ++i ) {
        _rotmatrix.push_back( (double)drot[i] );
    }
    _rotmatrixT = _rotmatrix;
    std::swap( _rotmatrixT[1], _rotmatrixT[3] );
    std::swap( _rotmatrixT[2], _rotmatrixT[6] );
    std::swap( _rotmatrixT[5], _rotmatrixT[7] );
}

/* Parameters */

void MGISBehaviour::setParameter( const std::string param, const ASTERDOUBLE value ) {
    VectorString vect = _behav->params;
    if ( std::find( vect.begin(), vect.end(), param ) != vect.end() ) {
        MGB::setParameter( *_behav, param, (double)value );
#ifdef ASTER_DEBUG_CXX
        std::cout << "MGISDBG: set '" << param << "' (double): " << (double)value << std::endl
                  << std::flush;
#endif
    } else {
        std::cout << "Parameter '" << param << "' (real) not found, ignored." << std::endl;
    }
}

template ASTERDOUBLE
MGISBehaviour::getMFrontParameter< ASTERDOUBLE, double >( const std::string param );
template ASTERINTEGER
MGISBehaviour::getMFrontParameter< ASTERINTEGER, int >( const std::string param );

template < typename T, typename U >
T MGISBehaviour::getMFrontParameter( const std::string param ) {
    return (T)MGB::getParameterDefaultValue< U >( *_behav, param );
}

void MGISBehaviour::setParameter( const std::string param, const ASTERINTEGER value ) {
    VectorString ivect = _behav->iparams;
    VectorString usvect = _behav->usparams;
    if ( std::find( ivect.begin(), ivect.end(), param ) != ivect.end() ) {
        MGB::setParameter( *_behav, param, (int)value );
#ifdef ASTER_DEBUG_CXX
        std::cout << "MGISDBG: set '" << param << "' (int): " << (int)value << std::endl
                  << std::flush;
#endif
    } else if ( std::find( usvect.begin(), usvect.end(), param ) != usvect.end() ) {
        MGB::setParameter( *_behav, param, (unsigned short)value );
#ifdef ASTER_DEBUG_CXX
        std::cout << "MGISDBG: set '" << param << "' (ushort): " << (unsigned short)value
                  << std::endl
                  << std::flush;
#endif
    } else {
        std::cout << "Parameter '" << param << "' (int) not found, ignored." << std::endl;
    }
}

/* Material Properties */

int MGISBehaviour::getNumberOfMaterialProperties() const { return _behav->mps.size(); }

VectorString MGISBehaviour::getMaterialPropertiesNames() const {
    VectorString names;
    for ( auto &var : _behav->mps ) {
        names.push_back( var.name );
    }
    return names;
}

void MGISBehaviour::setMaterialProperties( const StateId stid, const ASTERDOUBLE *values ) {
    MGB::State *state = getState( stid );
    for ( auto i = 0; i < getNumberOfMaterialProperties(); ++i ) {
        state->material_properties[i] = (double)values[i];
    }
}

/* General */
void MGISBehaviour::setInitialState( ASTERDOUBLE dt, ASTERDOUBLE *K ) {
    _data->dt = (double)dt;
    _data->K[0] = (double)K[0];
}

int MGISBehaviour::integrate() {
    /*
     * integrate returns:
     * - -1: integration failed
     * -  0: integration succeeded but results are unreliable
     * -  1: integration succeeded and results are reliable
     */
    if ( _rotmatrix.size() != 0 ) {
        // global frame to material frame
        MGB::rotateGradients( _data->s0.gradients, *_behav, _rotmatrix );
        MGB::rotateGradients( _data->s1.gradients, *_behav, _rotmatrix );
        // rotateThermodynamicForces transposes the rotation matrix
        MGB::rotateThermodynamicForces( _data->s0.thermodynamic_forces, *_behav, _rotmatrixT );
    }

    auto view = make_view( *_data );
    int iret = MGB::integrate( view, *_behav );

    if ( _rotmatrix.size() != 0 ) {
        // material frame to global frame
        MGB::rotateThermodynamicForces( _data->s1.thermodynamic_forces, *_behav, _rotmatrix );
        MGB::rotateTangentOperatorBlocks( _data->K, *_behav, _rotmatrix );
    }
    return iret;
}

void MGISBehaviour::setOutputs( ASTERDOUBLE *stress, ASTERDOUBLE *internal_variables,
                                ASTERDOUBLE *K, ASTERDOUBLE *rdt ) {
    convertTensorFromMgis( _data->s1.thermodynamic_forces, stress );
    for ( auto i = 0; i < getSizeOfInternalStateVariables(); ++i ) {
        internal_variables[i] = ( ASTERDOUBLE )( _data->s1.internal_state_variables[i] );
    }
    convertMatrixFromMgis( _data->K, K );
    *rdt = (ASTERDOUBLE)_data->rdt;
}

void MGISExceptionHandler() {
    try {
        throw;
    } catch ( std::exception &e ) {
        std::cout << e.what() << '\n';
    } catch ( ... ) {
        std::cout << "unknown exception thrown";
    }
    UTMESS( "F", "DVP_98" );
}

#endif // ASTER_HAVE_MGIS

void convertTensorFromMgis( const VectorReal &vect, ASTERDOUBLE *array ) {
    auto size = vect.size();
    // const double sqr2 = sqrt( 2.0 );
    for ( auto i = 0; i < size; ++i ) {
        array[i] = (ASTERDOUBLE)vect[i];
    }

    // switch ( size ) {
    // case 4: {
    //     array[3] = (ASTERDOUBLE)( vect[3] / sqr2 );
    //     break;
    // }
    // case 6: {
    //     for ( auto i : { 3, 4, 5 } ) {
    //         array[i] = (ASTERDOUBLE)( vect[i] / sqr2 );
    //     }
    //     break;
    // }
    // default:
    //     throw std::runtime_error( "wrong size: " + std::to_string( size ) );
    // }
}

void convertMatrixFromMgis( const VectorReal &vect, ASTERDOUBLE *array ) {
    auto size = vect.size();
    const double sqr2 = sqrt( 2.0 );
    if ( size != 81 ) {
        for ( auto i = 0; i < size; ++i ) {
            array[i] = (ASTERDOUBLE)vect[i];
        }
    }

    switch ( size ) {
    case 9:
        break;
    case 16:
        for ( auto i : { 3, 7, 11, 12, 13, 14 } ) {
            array[i] = ( ASTERDOUBLE )( vect[i] / sqr2 );
        }
        array[15] = ( ASTERDOUBLE )( vect[15] / 2.0 );
        break;
    case 25:
        break;
    case 36:
        for ( auto i : { 3, 4, 5, 9, 10, 11, 15, 16, 17, 18, 19, 20, 24, 25, 26, 30, 31, 32 } ) {
            array[i] = ( ASTERDOUBLE )( vect[i] / sqr2 );
        }
        for ( auto i : { 21, 22, 23, 27, 28, 29, 33, 34, 35 } ) {
            array[i] = ( ASTERDOUBLE )( vect[i] / 2.0 );
        }
        break;
    case 81: // dest is 6 x 3 x 3 (see lcmelas.F90)
        // xx, yy, zz, xy
        for ( auto i = 0; i < 36; ++i ) {
            array[i] = (ASTERDOUBLE)vect[i];
        }
        // skip yx
        // xz
        for ( auto i = 36; i < 45; ++i ) {
            array[i] = (ASTERDOUBLE)vect[i + 9];
        }
        // skip zx
        // yz
        for ( auto i = 45; i < 54; ++i ) {
            array[i] = (ASTERDOUBLE)vect[i + 18];
        }
        // skip zy
        break;
    default:
        throw std::runtime_error( "convertMatrixFromMgis: wrong size: " + std::to_string( size ) );
    }
}

void convertTensorToMgis( const ASTERDOUBLE *src, const int insize, VectorReal &dest,
                          const int destsize ) {
    switch ( destsize ) {
    case 5:
        dest[0] = src[0]; // 0, 0
        dest[1] = src[4]; // 1, 1
        dest[2] = src[8]; // 2, 2
        dest[3] = src[1]; // 0, 1
        dest[4] = src[3]; // 1, 0
        break;
    case 3:
    case 4:
    case 6:
        AS_ASSERT( destsize == insize );
        for ( auto i = 0; i < insize; ++i ) {
            dest[i] = (double)src[i];
        }
        break;
    case 9:
        dest[0] = src[0]; // 0, 0
        dest[1] = src[4]; // 1, 1
        dest[2] = src[8]; // 2, 2
        dest[3] = src[3]; // 0, 1
        dest[4] = src[1]; // 1, 0
        dest[5] = src[6]; // 0, 2
        dest[6] = src[2]; // 2, 0
        dest[7] = src[7]; // 1, 2
        dest[8] = src[5]; // 2, 1
        break;
    default:
        throw std::runtime_error(
            "convertTensorToMgis: unexpected sizes: size(src):" + std::to_string( insize ) +
            ", size(dest):" + std::to_string( destsize ) );
    }
}
