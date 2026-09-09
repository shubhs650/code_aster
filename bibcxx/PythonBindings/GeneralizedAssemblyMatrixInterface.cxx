/**
 * @file GeneralizedAssemblyMatrixInterface.cxx
 * @brief Interface python de GeneralizedAssemblyMatrix
 * @author Nicolas Sellenet
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

#include "PythonBindings/GeneralizedAssemblyMatrixInterface.h"

#include "aster_pybind.h"

#include "Numbering/GeneralizedDOFNumbering.h"
#include "PythonBindings/VariantModalBasisInterface.h"

void exportGeneralizedAssemblyMatrixToPython( py::module_ &mod ) {

    py::class_< GenericGeneralizedAssemblyMatrix, GenericGeneralizedAssemblyMatrixPtr,
                DataStructure >( mod, "GeneralizedAssemblyMatrix" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "getGeneralizedDOFNumbering",
              &GenericGeneralizedAssemblyMatrix::getGeneralizedDOFNumbering )
        .def( "getModalBasis", &getModalBasis< GenericGeneralizedAssemblyMatrixPtr > )
        .def( "setGeneralizedDOFNumbering",
              &GenericGeneralizedAssemblyMatrix::setGeneralizedDOFNumbering )
        .def( "setModalBasis", py::overload_cast< const GeneralizedModeResultPtr & >(
                                   &GenericGeneralizedAssemblyMatrix::setModalBasis ) )
        .def( "setModalBasis", py::overload_cast< const ModeResultPtr & >(
                                   &GenericGeneralizedAssemblyMatrix::setModalBasis ) )
        .def( "exists", &GenericGeneralizedAssemblyMatrix::exists, R"(
Return True if the matrix exists

Returns:
    bool: True if the matrix exists else False.
        )" )
        .def( "isDiagonal", &GenericGeneralizedAssemblyMatrix::isDiagonal, R"(
Return True if the matrix is diagonal

Returns:
    bool: True if the matrix is diagonal else False.
        )" )
        .def( "isDense", &GenericGeneralizedAssemblyMatrix::isDense, R"(
Return True if the matrix is dense

Returns:
    bool: True if the matrix is dense else False.
        )" )
        .def( "size", &GenericGeneralizedAssemblyMatrix::size, R"(
Return the size of the matrix

Returns:
    int: size of the matrix.
        )" );

    py::class_< GeneralizedAssemblyMatrixReal, GeneralizedAssemblyMatrixRealPtr,
                GenericGeneralizedAssemblyMatrix >( mod, "GeneralizedAssemblyMatrixReal" )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyMatrixReal > ) )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyMatrixReal, std::string > ) )
        .def( "allocate", &GeneralizedAssemblyMatrixReal::allocate, R"(
Allocate the matrix
        )",
              py::arg( "isSymmetric" ) = true )
        .def( "isSymmetric", &GeneralizedAssemblyMatrixReal::isSymmetric, R"(
Return True if the matrix is symmetric

Returns:
    bool: True if the matrix is symmetric else False.
        )" )
        .def( "getUpperValues", &GeneralizedAssemblyMatrixReal::getUpperValues, R"(
Return the upper part of the matrix.

Returns:
    list[float]: upper part of the matrix.
        )" )
        .def( "getLowerValues", &GeneralizedAssemblyMatrixReal::getLowerValues, R"(
Return the lower part of the matrix.

Returns:
    list[float]: lower part of the matrix.
        )" )
        .def( "setUpperValues", &GeneralizedAssemblyMatrixReal::setUpperValues, R"(
Set the upper part of the matrix.

Arguments:
    values [list[float]]: set upper part of the matrix.
        )",
              py::arg( "values" ) )
        .def( "setLowerValues", &GeneralizedAssemblyMatrixReal::setLowerValues, R"(
Set the lower part of the matrix.

Arguments:
    values [list[float]]: set lower part of the matrix.
        )",
              py::arg( "values" ) );

    py::class_< GeneralizedAssemblyMatrixComplex, GeneralizedAssemblyMatrixComplexPtr,
                GenericGeneralizedAssemblyMatrix >( mod, "GeneralizedAssemblyMatrixComplex" )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyMatrixComplex > ) )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyMatrixComplex, std::string > ) )
        .def( "isSymmetric", &GeneralizedAssemblyMatrixComplex::isSymmetric, R"(
Return True if the matrix is symmetric

Returns:
    bool: True if the matrix is symmetric else False.
        )" )
        .def( "getUpperValues", &GeneralizedAssemblyMatrixComplex::getUpperValues, R"(
Return the upper part of the matrix.

Returns:
    list[complex]: upper part of the matrix.
        )" )
        .def( "getLowerValues", &GeneralizedAssemblyMatrixComplex::getLowerValues, R"(
Return the lower part of the matrix.

Returns:
    list[complex]: lower part of the matrix.
        )" )
        .def( "setUpperValues", &GeneralizedAssemblyMatrixComplex::setUpperValues, R"(
Set the upper part of the matrix.

Arguments:
    values [list[complex]]: set upper part of the matrix.
        )",
              py::arg( "values" ) )
        .def( "setLowerValues", &GeneralizedAssemblyMatrixComplex::setLowerValues, R"(
Set the lower part of the matrix.

Arguments:
    values [list[complex]]: set lower part of the matrix.
        )",
              py::arg( "values" ) );
};
