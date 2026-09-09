/**
 * @file MeshInterface.cxx
 * @brief Interface python de Mesh
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

#include "PythonBindings/MeshInterface.h"

#include "aster_pybind.h"

#include <Meshes/BaseMesh.h>
#include <Meshes/Mesh.h>
#include <Meshes/MeshEntities.h>

void exportMeshToPython( py::module_ &mod ) {

    py::class_< Mesh, Mesh::MeshPtr, BaseMesh >( mod, "Mesh" )
        .def( py::init( &initFactoryPtr< Mesh > ) )
        .def( py::init( &initFactoryPtr< Mesh, std::string > ) )
        .def( "getGroupsOfCells", &Mesh::getGroupsOfCells, R"(
Return the list of the existing groups of cells.

Returns:
    list[str]: List of groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfCells", &Mesh::hasGroupOfCells, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.
    local (bool): not used (for compatibilty with ParallelMesh)

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "setGroupOfCells", &Mesh::setGroupOfCells, R"(
Set new group of cells in the mesh

Arguments:
    group_name (str): Name of the new group.
    cell_ids (list[int]) : cell ids which are in the group
        )",
              py::arg( "group_name" ), py::arg( "cell_ids" ) )
        .def( "getCells", py::overload_cast< const std::string >( &Mesh::getCells, py::const_ ),
              R"(
Return the list of the indexes of the cells that belong to a group of cells.

Arguments:
    group_name (str): Name of the local group.

Returns:
    list[int]: Indexes of the cells of the local group.
        )",
              py::arg( "group_name" ) )
        .def( "getCells", py::overload_cast< const VectorString & >( &Mesh::getCells, py::const_ ),
              R"(
Return the list of the indexes of the cells that belong to the groups of cells.

Arguments:
    groups_name (str): Name of the local groups.

Returns:
    list[int]: Indexes of the cells of the local groups.
        )",
              py::arg( "groups_name" ) = VectorString() )
        .def( "getGroupsOfNodes", &Mesh::getGroupsOfNodes, R"(
Return the list of the existing groups of nodes.

Arguments:
    local (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[str]: List of groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfNodes", &Mesh::hasGroupOfNodes, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.
    local (bool): not used (for compatibilty with ParallelMesh)

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "setGroupOfNodes", &Mesh::setGroupOfNodes, R"(
Set new group of nodes in the mesh

Arguments:
    group_name (str): Name of the new group.
    node_ids (list[int]) : node ids which are in the group
    localNumbering=false (bool): not used (for compatibilty with ParallelMesh)
        )",
              py::arg( "group_name" ), py::arg( "node_ids" ), py::arg( "localNumbering" ) = false )
        .def( "_getNodes",
              py::overload_cast< const VectorString &, const bool, const ASTERINTEGER >(
                  &Mesh::getNodes, py::const_ ),
              R"(
Return the list of the indexes of the nodes that belong to groups of nodes.
* For internal use only *

Arguments:
    group_name (list[str]): Name of groups (default: "").
    localNumbering (bool): not used (for compatibilty with ParallelMesh)
    same_rank (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: Indexes of the nodes of groups.
        )",
              py::arg( "group_name" ) = VectorString(), py::arg( "localNumbering" ) = false,
              py::arg( "same_rank" ) = PythonBool::None )
        .def( "getInnerNodes", &Mesh::getInnerNodes, R"(
Return the list of the indexes of the nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )" )
        .def( "readAsterFile", &Mesh::readAsterFile, R"(
Read a mesh file from ASTER format.

Arguments:
    filename (Path|str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              py::arg( "filename" ) )
        .def( "readGibiFile", &Mesh::readGibiFile, R"(
Read a mesh file from GIBI format.

Arguments:
    filename (Path|str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              py::arg( "filename" ) )
        .def( "readGmshFile", &Mesh::readGmshFile, R"(
Read a mesh file from GMSH format.

Arguments:
    filename (Path|str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              py::arg( "filename" ) )
        .def( "isQuadratic", &Mesh::isQuadratic, R"(
Tells if the mesh contains quadratic cells.

Arguments:
    local (bool): not used (for compatibilty with ParallelMesh)

Returns:
    bool: *True* if the mesh contains quadratic cells, *False* otherwise.
        )",
              py::arg( "local" ) = false )
        .def( "_getNodesFromCells",
              py::overload_cast< const VectorString &, const bool, const ASTERINTEGER >(
                  &Mesh::getNodesFromCells, py::const_ ),
              R"(
Returns the nodes indexes of a group of cells.
For developpers only.

Arguments:
    group_name (str): name of the group of cells.
    localNumbering (bool): not used (for compatibilty with ParallelMesh)
    same_rank (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: indexes of the nodes.
        )",
              py::arg( "group_name" ), py::arg( "localNumbering" ) = false,
              py::arg( "same_rank" ) = PythonBool::None )
        .def( "fix", &Mesh::fix, R"(
Fix potential problems.

Arguments:
    remove_orphan (bool) : remove orphelan nodes.
    positive_measure (bool) : reorder nodes to have a positive measure of cells.
    outward_normal (bool) : reorder nodes to have an outward normal for boundary faces.
    double_nodes (bool) : merge double nodes with almost same coordinates.
    double_cells (bool) : merge double cells with same nodes.
    tole (float) : tolerance for double nodes
    info (int) : verbosity mode (0 or 1 or 2).

Returns:
    Mesh: fixed mesh
        )",
              py::arg( "remove_orphan" ) = true, py::arg( "positive_measure" ) = true,
              py::arg( "outward_normal" ) = true, py::arg( "double_nodes" ) = true,
              py::arg( "double_cells" ) = true, py::arg( "tole" ) = 1e-7, py::arg( "info" ) = 1 )
        .def( "getOctreeMesh", &Mesh::getOctreeMesh, R"(
Get the octree mesh.

Arguments:
    nb_max_pt (int) : maximum number of points for the last level.
    nb_max_level (int) : maximum number of level.

Returns:
    Mesh: octree mesh.
        )",
              py::arg( "nb_max_pt" ) = 1, py::arg( "nb_max_level" ) = 20 )
        .def( "convertToLinear", &Mesh::convertToLinear, R"(
Convert the mesh to a linear one.

Arguments:
    info (int) : verbosity mode (1 or 2). Default 1.

Returns:
    Mesh: the linearized mesh.
        )",
              py::arg( "info" ) = 1 )
        .def( "convertToBiQuadratic", &Mesh::convertToBiQuadratic, R"(
Convert the mesh to a bi-quadratic one.
For cells that have no bi-quadratic version, the quadratic version is used.

Arguments:
    info (int) : verbosity mode (1 or 2). Default 1.

Returns:
    Mesh: the bi-quadratic mesh.
        )",
              py::arg( "info" ) = 1 )
        .def( "convertToQuadratic", &Mesh::convertToQuadratic, R"(
Convert the mesh to a quadratic one.

Arguments:
    info (int) : verbosity mode (1 or 2). Default 1.

Returns:
    Mesh: the quadratic mesh.
        )",
              py::arg( "info" ) = 1 )
        .def( "convertToCubic", &Mesh::convertToCubic, R"(
Convert the mesh to a cubic one.

Arguments:
    info (int) : verbosity mode (1 or 2). Default 1.

Returns:
    Mesh: the cubic mesh.
        )",
              py::arg( "info" ) = 1 )
        .def( "addNodeLabels", &Mesh::addNodeLabels, R"(
      Add node labels.

      Arguments:
          node_labels (list) : Node labels.
              )",
              py::arg( "node_labels" ) )
        .def( "addCellLabels", &Mesh::addCellLabels, R"(
      Add cell labels.

      Arguments:
          cell_labels (list) : Cell labels.
              )",
              py::arg( "cell_labels" ) )
        .def( "printMedFile", &Mesh::printMedFile, R"(
Print the mesh in the MED format

Arguments:
    filename (Path|str): Name of the file
    local (bool=True) : print local values only (relevant for a ParallelMesh only)
    version (list): list of size 3 ([major, minor, release])
    mode (str): 'a' or 'w' ('w' to overwrite file)

Returns:
    Bool: True if of
            )",
              py::arg( "fileName" ), py::arg( "local" ) = true,
              py::arg( "version" ) = std::array< int, 3 >( { 0, 0, 0 } ), py::arg( "mode" ) = 'a' );
};
