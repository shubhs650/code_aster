/**
 * @file AsterToMedWriterInterface.cxx
 * @brief Interface python de AsterToMedWriter
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

#include "PythonBindings/AsterToMedWriterInterface.h"

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

void exportAsterToMedWriterToPython( py::module_ &mod ) {

    py::class_< AsterToMedWriter, AsterToMedWriter::AsterToMedWriterPtr > c1( mod,
                                                                              "AsterToMedWriter" );
    c1.def( py::init( &initFactoryPtr< AsterToMedWriter > ) );
    c1.def( "__pickling_disabled__", disable_pickling< AsterToMedWriter >() );

#ifdef ASTER_HAVE_MED
    c1.def( "printMesh",
            py::overload_cast< const MeshPtr &, const std::filesystem::path &, bool,
                               const std::string &, const std::string & >(
                &AsterToMedWriter::printMesh ),
            R"(
Print mesh to med file

Arguments:
    Mesh: mesh to print
    path (Path|str): path to med file
    parallelPrint (bool): false by default. If true print in one parallel file (optional)
    mesh_name (str): mesh name (optional)
    mode (str): 'a' or 'w' ('w' to overwrite file) (optional)
            )",
            py::arg( "mesh" ), py::arg( "path" ), py::arg( "parallelPrint" ) = false,
            py::arg( "mesh_name" ) = "", py::arg( "mode" ) = 'a' );
#ifdef ASTER_HAVE_MPI
    c1.def( "printMesh",
            py::overload_cast< const ParallelMeshPtr &, const std::filesystem::path &, bool,
                               const std::string &, const std::string & >(
                &AsterToMedWriter::printMesh ),
            R"(
Print mesh to med file

Arguments:
    Mesh: mesh to print
    path (Path|str): path to med file
    parallelPrint (bool): false by default. If true print in one parallel file (optional)
    mesh_name (str): mesh name (optional)
    mode (str): 'a' or 'w' ('w' to overwrite file) (optional)
            )",
            py::arg( "mesh" ), py::arg( "path" ), py::arg( "parallelPrint" ) = false,
            py::arg( "mesh_name" ) = "", py::arg( "mode" ) = 'a' );

    c1.def( "printMesh",
            py::overload_cast< const ConnectionMeshPtr &, const std::filesystem::path &, bool,
                               const std::string &, const std::string & >(
                &AsterToMedWriter::printMesh ),
            R"(
Print mesh to med file

Arguments:
    Mesh: mesh to print
    path (Path|str): path to med file
    parallelPrint (bool): false by default. If true print in one parallel file (optional)
    mesh_name (str): mesh name (optional)
    mode (str): 'a' or 'w' ('w' to overwrite file) (optional)
            )",
            py::arg( "mesh" ), py::arg( "path" ), py::arg( "parallelPrint" ) = false,
            py::arg( "mesh_name" ) = "", py::arg( "mode" ) = 'a' );
#endif
    c1.def( "printResult", &AsterToMedWriter::printResult,
            R"(
Print result to med file

Arguments:
    result: result to print
    path (Path|str): path to med file
    parallelPrint (bool): false by default. If true print in one parallel file (optional)
    mode (str): 'a' or 'w' ('w' to overwrite file) (optional)
            )",
            py::arg( "result" ), py::arg( "path" ), py::arg( "parallelPrint" ) = false,
            py::arg( "mode" ) = 'a' );
#endif
};
