/**
 * @file GeneralizedLoadInterface.cxx
 * @brief Interface python de GeneralizedLoad
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

#include "PythonBindings/GeneralizedLoadInterface.h"

#include "aster_pybind.h"

void exportGeneralizedLoadToPython( py::module_ &mod ) {

    py::class_< GeneralizedLoad, GeneralizedLoad::GeneralizedLoadPtr, DataStructure >(
        mod, "GeneralizedLoad" )
        .def( py::init( &initFactoryPtr< GeneralizedLoad > ) )
        .def( py::init( &initFactoryPtr< GeneralizedLoad, std::string > ) )
        .def( "setDOFNumbering", &GeneralizedLoad::setDOFNumbering )
        .def( "getDOFNumbering", &GeneralizedLoad::getDOFNumbering )
        .def( "getAssemblyVector", &GeneralizedLoad::getAssemblyVector )
        .def( "getMPCs", &GeneralizedLoad::getMPCs )
        .def( "setMPCs", &GeneralizedLoad::setMPCs );
};
