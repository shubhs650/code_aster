#ifndef EXCEPTION_H_
#define EXCEPTION_H_

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

#ifdef __cplusplus

#include "astercxx.h"

#include "aster_pybind.h"

#include <exception>
#include <string>

class AbstractErrorCpp : public std::exception {
  private:
    std::string _idmess;
    VectorString _valk;
    VectorLong _vali;
    VectorReal _valr;

  public:
    AbstractErrorCpp( const std::string idmess, const VectorString &valk = {},
                      const VectorLong &vali = {}, const VectorReal &valr = {} )
        : _idmess( idmess ), _valk( valk ), _vali( vali ), _valr( valr ) {}

    const char *what() const noexcept { return _idmess.c_str(); }

    virtual ~AbstractErrorCpp() noexcept {}

    /* Build arguments for the Python exception */
    PyObject *py_attrs() const;
};

// Subclasses
template < int Id >
class ErrorCpp : public AbstractErrorCpp {
  private:
    int _id = Id;

  public:
    ErrorCpp( const std::string idmess, const VectorString &valk = {}, const VectorLong &vali = {},
              const VectorReal &valr = {} )
        : AbstractErrorCpp( idmess, valk, vali, valr ) {}

    int id() const { return _id; }
};

using AsterErrorCpp = ErrorCpp< ASTER_ERROR >;
using ConvergenceErrorCpp = ErrorCpp< ASTER_CONVERGENCE_ERROR >;
using IntegrationErrorCpp = ErrorCpp< ASTER_INTEGRATION_ERROR >;
using SolverErrorCpp = ErrorCpp< ASTER_SOLVER_ERROR >;
using ContactErrorCpp = ErrorCpp< ASTER_CONTACT_ERROR >;
using TimeLimitErrorCpp = ErrorCpp< ASTER_TIMELIMIT_ERROR >;

void createExceptions( py::module_ &mod );

void raiseAsterError( const std::string idmess = "VIDE_1", const VectorString &valk = {},
                      const VectorLong &vali = {}, const VectorReal &valr = {} );

extern "C" void DEFPSPSPPPP( UEXCEP, uexcep, _IN ASTERINTEGER *exc_id, _IN char *idmess,
                             _IN STRING_SIZE lidmess, _IN ASTERINTEGER *nbk, _IN char *valk,
                             _IN STRING_SIZE lvk, _IN ASTERINTEGER *nbi, _IN ASTERINTEGER *vali,
                             _IN ASTERINTEGER *nbr, _IN ASTERDOUBLE *valr );

extern "C" {
#endif // __cplusplus
void _raiseException();

#ifdef __cplusplus
}
#endif

#endif
