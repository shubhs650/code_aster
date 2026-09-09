/**
 * @file DataFieldInterface.cxx
 * @brief Interface python de DataField
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

#include "PythonBindings/DataFieldInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/DataStructureInterface.h"

void exportDataFieldToPython( py::module_ &mod ) {

    py::class_< DataField, DataFieldPtr, DSWithCppPickling >( mod, "DataField" )
        .def( py::init( &initFactoryPtr< DataField > ) )
        .def( py::init( &initFactoryPtr< DataField > ) )
        .def( py::init( &initFactoryPtr< DataField, std::string > ) )
        .def( py::init( &initFactoryPtr< DataField, std::string, std::string > ) )
        .def( "getFieldType", &DataField::getFieldType, R"(
            Get field type between "ELEM", "ELGA", "ELNO", "NOEU", "CART".

            Returns:
                str: field type
            )" );
};
