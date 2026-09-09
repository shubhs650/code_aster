/**
 * @file Exception.h
 * @brief Definition of code_aster exceptions
 * @author Mathieu Courtois
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Supervis/Exceptions.h"

#include "aster_pybind.h"

#include "Utilities/Tools.h"

#include <exception>
#include <string>

#include <assert.h>

PyObject *AbstractErrorCpp::py_attrs() const {
    int idx = 0;
    PyObject *py_err = PyTuple_New( 4 );

    PyObject *py_id = PyUnicode_FromString( _idmess.c_str() );

    PyObject *py_valk = PyTuple_New( _valk.size() );
    idx = 0;
    for ( auto it = _valk.begin(); it != _valk.end(); ++it, ++idx )
        PyTuple_SetItem( py_valk, idx, PyUnicode_FromString( it->c_str() ) );

    PyObject *py_vali = PyTuple_New( _vali.size() );
    idx = 0;
    for ( auto it = _vali.begin(); it != _vali.end(); ++it, ++idx )
        PyTuple_SetItem( py_vali, idx, PyLong_FromLong( *it ) );

    PyObject *py_valr = PyTuple_New( _valr.size() );
    idx = 0;
    for ( auto it = _valr.begin(); it != _valr.end(); ++it, ++idx )
        PyTuple_SetItem( py_valr, idx, PyFloat_FromDouble( *it ) );

    PyTuple_SetItem( py_err, 0, py_id );
    PyTuple_SetItem( py_err, 1, py_valk );
    PyTuple_SetItem( py_err, 2, py_vali );
    PyTuple_SetItem( py_err, 3, py_valr );

    return py_err;
}

void createExceptions( py::module_ &mod ) {
    static py::exception< AsterErrorCpp > excAsterError( mod, "AsterError" );
    PyObject *pyError = excAsterError.ptr();
    static py::exception< ConvergenceErrorCpp > excConvergenceError( mod, "ConvergenceError",
                                                                     pyError );
    static py::exception< IntegrationErrorCpp > excIntegrationError( mod, "IntegrationError",
                                                                     pyError );
    static py::exception< SolverErrorCpp > excSolverError( mod, "SolverError", pyError );
    static py::exception< ContactErrorCpp > excContactError( mod, "ContactError", pyError );
    static py::exception< TimeLimitErrorCpp > excTimeLimitError( mod, "TimeLimitError", pyError );

    py::register_exception_translator( []( std::exception_ptr ptr ) {
        try {
            if ( ptr )
                std::rethrow_exception( ptr );
        } catch ( const ConvergenceErrorCpp &exc ) {
            PyObject *args = exc.py_attrs();
            PyErr_SetObject( excConvergenceError.ptr(), args );
            Py_DECREF( args );
        } catch ( const IntegrationErrorCpp &exc ) {
            PyObject *args = exc.py_attrs();
            PyErr_SetObject( excIntegrationError.ptr(), args );
            Py_DECREF( args );
        } catch ( const SolverErrorCpp &exc ) {
            PyObject *args = exc.py_attrs();
            PyErr_SetObject( excSolverError.ptr(), args );
            Py_DECREF( args );
        } catch ( const ContactErrorCpp &exc ) {
            PyObject *args = exc.py_attrs();
            PyErr_SetObject( excContactError.ptr(), args );
            Py_DECREF( args );
        } catch ( const TimeLimitErrorCpp &exc ) {
            PyObject *args = exc.py_attrs();
            PyErr_SetObject( excTimeLimitError.ptr(), args );
            Py_DECREF( args );
        } catch ( const AsterErrorCpp &exc ) {
            PyObject *args = exc.py_attrs();
            PyErr_SetObject( excAsterError.ptr(), args );
            Py_DECREF( args );
        }
    } );
}

void raiseAsterError( const std::string idmess, const VectorString &valk, const VectorLong &vali,
                      const VectorReal &valr ) {
#ifdef ASTER_DEBUG_CXX
    std::cout << "Raising C++ AsterError with id '" << idmess << "'..." << std::endl;
#endif
    throw AsterErrorCpp( idmess, valk, vali, valr );
}

extern "C" void DEFPSPSPPPP( UEXCEP, uexcep, _IN ASTERINTEGER *exc_id, _IN char *idmess,
                             _IN STRING_SIZE lidmess, _IN ASTERINTEGER *nbk, _IN char *valk,
                             _IN STRING_SIZE lvk, _IN ASTERINTEGER *nbi, _IN ASTERINTEGER *vali,
                             _IN ASTERINTEGER *nbr, _IN ASTERDOUBLE *valr ) {
    VectorString argk = {};
    VectorLong argi = {};
    VectorReal argr = {};
    for ( int i = 0; i < *nbk; ++i ) {
        argk.push_back( strip( std::string( valk + i * lvk, lvk ) ) );
    }
    for ( int i = 0; i < *nbi; ++i ) {
        argi.push_back( vali[i] );
    }
    for ( int i = 0; i < *nbr; ++i ) {
        argr.push_back( valr[i] );
    }

    // The identifier of each Python exception is defined in 'LibAster.cxx'
    std::string idm( strip( std::string( idmess, lidmess ) ) );

    switch ( *exc_id ) {
    case ASTER_CONVERGENCE_ERROR:
        throw ErrorCpp< ASTER_CONVERGENCE_ERROR >( idm, argk, argi, argr );

    case ASTER_INTEGRATION_ERROR:
        throw ErrorCpp< ASTER_INTEGRATION_ERROR >( idm, argk, argi, argr );

    case ASTER_SOLVER_ERROR:
        throw ErrorCpp< ASTER_SOLVER_ERROR >( idm, argk, argi, argr );

    case ASTER_CONTACT_ERROR:
        throw ErrorCpp< ASTER_CONTACT_ERROR >( idm, argk, argi, argr );

    case ASTER_TIMELIMIT_ERROR:
        throw ErrorCpp< ASTER_TIMELIMIT_ERROR >( idm, argk, argi, argr );

    default:
        throw AsterErrorCpp( idm, argk, argi, argr );
    }
}

// TODO for aster_exceptions, to be removed in the future
extern "C" void _raiseException() { raiseAsterError(); }
