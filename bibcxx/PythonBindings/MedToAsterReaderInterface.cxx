/**
 * @file MedToAsterReaderInterface.cxx
 * @brief Interface python de MedToAsterReader
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

#include "PythonBindings/MedToAsterReaderInterface.h"

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

void exportMedToAsterReaderToPython( py::module_ &mod ) {

    py::class_< MedToAsterReader, MedToAsterReader::MedToAsterReaderPtr > c1( mod,
                                                                              "MedToAsterReader" );
    c1.def( py::init( &initFactoryPtr< MedToAsterReader > ) );
    c1.def( "__pickling_disabled__", disable_pickling< MedToAsterReader >() );

#ifdef ASTER_HAVE_MED
    c1.def( "readMeshFromMedFile", &MedToAsterReader::readMeshFromMedFile,
            R"(
Open med file

Arguments:
    Mesh: return mesh to fill
    path (Path|str): path to med file
    mesh_name (str): mesh name (optional)
    verbosity (int): verbosity (optional)
            )",
            py::arg( "mesh" ), py::arg( "path" ), py::arg( "mesh_name" ) = "",
            py::arg( "verbosity" ) = 0 );
#ifdef ASTER_HAVE_MPI
    c1.def( "readIncompleteMeshFromMedFile", &MedToAsterReader::readIncompleteMeshFromMedFile,
            R"(
      Open med file

      Arguments:
          IncompleteMesh: return mesh to fill
          path (Path|str): path to med file
          mesh_name (str): mesh name (optional)
          verbosity (int): verbosity (optional)
                  )",
            py::arg( "mesh" ), py::arg( "path" ), py::arg( "mesh_name" ) = "",
            py::arg( "verbosity" ) = 0 );
    c1.def( "readParallelMeshFromMedFile", &MedToAsterReader::readParallelMeshFromMedFile,
            R"(
      Open med file

      Arguments:
          ParallelMesh: return mesh to fill
          path (Path|str): path to med file
          mesh_name (str): mesh name (optional)
          verbosity (int): verbosity (optional)
                  )",
            py::arg( "mesh" ), py::arg( "path" ), py::arg( "mesh_name" ) = "",
            py::arg( "verbosity" ) = 0 );
#endif
#endif
};
