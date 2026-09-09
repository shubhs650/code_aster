/**
 * @file ConstantFieldOnCellsInterface.cxx
 * @brief Interface python de ConstantFieldOnCells
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

#include "PythonBindings/ConstantFieldOnCellsInterface.h"

#include "aster_pybind.h"

#include "DataFields/FieldConverter.h"
#include "PythonBindings/DataFieldInterface.h"
#include "PythonBindings/DataStructureInterface.h"

void exportConstantFieldOnCellsToPython( py::module_ &mod ) {

    py::class_< ConstantFieldValues< double >, std::shared_ptr< ConstantFieldValues< double > > >(
        mod, "ConstantFieldValuesReal" )
        // fake initFactoryPtr
        // fake initFactoryPtr
        .def( "getValues", &ConstantFieldValues< double >::getValues, R"(
            Return the field values

            Returns:
                list[float]: List of values
        )" );

    py::class_< ConstantFieldOnCellsReal, ConstantFieldOnCellsRealPtr, DataField >(
        mod, "ConstantFieldOnCellsReal" )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsReal, BaseMeshPtr > ) )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsReal, const py::tuple & > ) )
        .def( define_pickling< ConstantFieldOnCellsReal >() )
        .def( "getMesh", &ConstantFieldOnCellsReal::getMesh )
        .def( "size", &ConstantFieldOnCellsReal::size, R"(
            Return the size of field

            Returns:
                int: size of field
        )" )
        .def( "getValues",
              py::overload_cast< const int & >( &ConstantFieldOnCellsReal::getValues, py::const_ ),
              R"(
            Return the field values

            Returns:
                list[float]: List of values
        )" )
        .def( "setValueOnCells", &ConstantFieldOnCellsReal::setValueOnCells )
        .def(
            "toSimpleFieldOnCells",
            []( const ConstantFieldOnCellsReal &f, const SimpleFieldOnCellsReal &sfm ) {
                return toSimpleFieldOnCells( f, sfm );
            },
            R"(
Convert to SimpleFieldOnCells

Returns:
    SimpleFieldOnCellsReal: field converted
        )" );

    py::class_< ConstantFieldOnCellsChar16, ConstantFieldOnCellsChar16Ptr, DataField >(
        mod, "ConstantFieldOnCellsChar16" )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsChar16, BaseMeshPtr > ) )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsChar16, const py::tuple & > ) )
        .def( define_pickling< ConstantFieldOnCellsChar16 >() )
        .def( "getMesh", &ConstantFieldOnCellsChar16::getMesh );

    py::class_< ConstantFieldOnCellsLong, ConstantFieldOnCellsLongPtr, DataField >(
        mod, "ConstantFieldOnCellsLong" )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsLong, BaseMeshPtr > ) )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsLong, const py::tuple & > ) )
        .def( define_pickling< ConstantFieldOnCellsLong >() )
        .def( "getMesh", &ConstantFieldOnCellsLong::getMesh );
};
