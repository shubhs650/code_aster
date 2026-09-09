/**
 * @file ParallelMeshInterface.cxx
 * @brief Interface python de ParallelMesh
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

#include "PythonBindings/ParallelMeshInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelMeshToPython( py::module_ &mod ) {

    py::class_< ParallelMesh, ParallelMesh::ParallelMeshPtr, BaseMesh >( mod, "ParallelMesh" )
        .def( py::init( &initFactoryPtr< ParallelMesh > ) )
        .def( py::init( &initFactoryPtr< ParallelMesh, std::string > ) )
        .def( "getGroupsOfCells", &ParallelMesh::getGroupsOfCells, R"(
Return the list of the existing (local or global) groups of cells.

Arguments:
    local (bool): search in local or global groups

Returns:
    list[str]: List of (local or global) groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "isQuadratic", &ParallelMesh::isQuadratic, R"(
Tells if the mesh contains quadratic cells.

Arguments:
    local (bool): if *True* only local cells are checked.

Returns:
    bool: *True* if the mesh contains quadratic cells, *False* otherwise.
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfCells", &ParallelMesh::hasGroupOfCells, R"(
The global group exists in the mesh

Arguments:
    group_name (str): Name of the global group.
    local (bool): search in local or global groups

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "setGroupOfCells", &ParallelMesh::setGroupOfCells, R"(
Set new group of cells in the mesh

Arguments:
    group_name (str): Name of the new group.
    cell_ids (list[int]) : cell ids which are in the group
        )",
              py::arg( "group_name" ), py::arg( "cell_ids" ) )
        .def( "setLastGhostsLayer", &ParallelMesh::setLastGhostsLayer, R"(
Set ids in local numbering of ghost nodes on the last layer

Arguments:
    list[int]: List of ghost nodes ids.
        )",
              py::arg( "node_ids" ) )
        .def( "getCells",
              py::overload_cast< const std::string >( &ParallelMesh::getCells, py::const_ ),
              R"(
Return the list of the indexes of the cells that belong to a group of cells.

Arguments:
    group_name (str): Name of the local group.

Returns:
    list[int]: Indexes of the cells of the local group.
        )",
              py::arg( "group_name" ) )
        .def( "getCells",
              py::overload_cast< const VectorString & >( &ParallelMesh::getCells, py::const_ ),
              R"(
Return the list of the indexes of the cells that belong to the groups of cells.

Arguments:
    groups_name (str): Name of the local groups.

Returns:
    list[int]: Indexes of the cells of the local groups.
        )",
              py::arg( "groups_name" ) = VectorString() )
        .def( "getLastGhostsLayer", &ParallelMesh::getLastGhostsLayer, R"(
Return ids in local numbering of ghost nodes on the last layer

Returns:
    list[int]: List of Nodes ids.
        )" )
        .def( "getGroupsOfNodes", &ParallelMesh::getGroupsOfNodes, R"(
Return the list of the existing (local or global) groups of nodes.

Arguments:
    local (bool): search in local or global groups

Returns:
    list[str]: List of (local or global) groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfNodes", &ParallelMesh::hasGroupOfNodes, R"(
The (local or global) group exists in the mesh

Arguments:
    group_name (str): Name of the (local or global) group.
    local (bool): search local or global groups

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "setGroupOfNodes", &ParallelMesh::setGroupOfNodes, R"(
Set new group of nodes in the mesh

Arguments:
    group_name (str): Name of the new group.
    node_ids (list[int]) : node ids which are in the group
    localNumbering=false (bool): ids are given in the local numbering ?
        )",
              py::arg( "group_name" ), py::arg( "node_ids" ), py::arg( "localNumbering" ) = false )

        .def( "_getNodes",
              py::overload_cast< const VectorString &, const bool, const ASTERINTEGER >(
                  &ParallelMesh::getNodes, py::const_ ),
              R"(
Return the list of the indexes of the nodes that belong to a group of nodes
with (local or global) indexing and a restriction to MPI-rank.
* For internal use only *

Arguments:
    group_name (str): Name of the group (default: "" = all nodes).
    localNumbering (bool) : use local or global numbering (default: True)
    same_rank : - None: keep all nodes (default: None)
                - True: keep the nodes which are owned by the current MPI-rank
                - False: keep the nodes which are not owned by the current MPI-rank

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              py::arg( "group_name" ) = VectorString(), py::arg( "localNumbering" ) = false,
              py::arg( "same_rank" ) = PythonBool::None )
        .def( "getInnerNodes", &ParallelMesh::getInnerNodes, R"(
Return the list of the indexes of the inner nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )" )
        .def( "getOuterNodes", &ParallelMesh::getOuterNodes, R"(
Return the list of the indexes of the outer nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )" )
        .def( "getInnerCells", &ParallelMesh::getInnerCells, R"(
Return the list of the indexes of the inner cells in the mesh

Returns:
    list[int]: Indexes of the cells.
        )" )
        .def( "getOuterCells", &ParallelMesh::getOuterCells, R"(
Return the list of the indexes of the outer cells in the mesh

Returns:
    list[int]: Indexes of the cells.
        )" )
        .def( "getNodesOwner", &ParallelMesh::getNodesOwner, R"(
Return the rank of the processor which owns the nodes

Returns:
    list[int]: MPI-Rank of the owner of the nodes
        )" )
        .def( "getNodesRanks", &ParallelMesh::getNodesRanks, R"(
Return the rank of the sub-domains which have the nodes.
The first subdomain given for a node is its owner.

Returns:
    list[list[int]]: MPI-Rank of of the subdomains
        )" )
        .def( "getCellsRanks", &ParallelMesh::getCellsRanks, R"(
Return the rank of the sub-domains which have the cells.
The first subdomain given for a cell is its owner.

Returns:
    list[list[int]]: MPI-Rank of of the subdomains
        )" )
        .def( "getCellsOwner", &ParallelMesh::getCellsOwner, R"(
Return the rank of the processor which owns the cells

Returns:
    list[int]: MPI-Rank of the owners of the cells
        )" )
        .def( "_updateGlobalGroupOfCells", &ParallelMesh::updateGlobalGroupOfCells, R"(
Share and update global groups of cells between MPI process.

This function has to be used by developer only and not user

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )" )
        .def( "_updateGlobalGroupOfNodes", &ParallelMesh::updateGlobalGroupOfNodes, R"(
Share and update global groups of nodes between MPI process.

This function has to be used by developer only and not user

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )" )
        .def( "_getNodesFromCells",
              py::overload_cast< const VectorString &, const bool, const ASTERINTEGER >(
                  &ParallelMesh::getNodesFromCells, py::const_ ),
              R"(
Returns the nodes indexes of a group of cells.
For developers only.

Arguments:
    group_name (str): Name of the group.
    localNumbering (bool) : use local or global numbering (default: True)
    same_rank : - None: keep all nodes (default: None)
                - True keep the nodes which are owned by the current MPI-rank
                - False: keep the nodes which are not owned by the current MPI-rank

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              py::arg( "group_name" ), py::arg( "localNumbering" ) = false,
              py::arg( "same_rank" ) = PythonBool::None )
        .def(
            "getGlobalToLocalNodeIds",
            []( const ParallelMesh &pm ) -> MapLong { return *pm.getGlobalToLocalNodeIds(); },
            R"(
        Returns global to local IDs mapping for nodes

        Returns:
            dict[int]: global to local IDs mapping.
            )" )
        .def( "getOppositeDomains", &ParallelMesh::getOppositeDomains,
              R"(
Returns the list of opposite domains of local process
        )" )
        .def( "getSendJoint", &ParallelMesh::getSendJoint,
              R"(
Returns ids of nodes in joint (inner nodes) for an opposite process

Arguments:
    rank: Rank of opposite domain
        )",
              py::arg( "rank" ) )
        .def( "getReceiveJoint", &ParallelMesh::getReceiveJoint,
              R"(
Returns ids of nodes in joint (inner nodes) for an opposite process

Arguments:
    rank: Rank of opposite domain
        )",
              py::arg( "rank" ) )
        .def( "_create_joints", &ParallelMesh::create_joints, R"(
Create the joints between domains (*for internal use*).

Arguments:
    domains (list[int]): Names of the remote domains.
    globalNodeIds (list[int]): Global IDs of each node.
    nodesOwner (list[int]): Owner of each node.
    globalCellIds (list[int]): Global IDs of each cell.
    joints (list[list[int]]): Definition of *E* mission and *R* eception joints.
    nbLayer (int]): Number of ghost layers.
        )",
              py::arg( "domains" ), py::arg( "globalNodeIds" ), py::arg( "nodesOwner" ),
              py::arg( "globalCellIds" ), py::arg( "joints" ), py::arg( "nbLayer" ) )
        .def( "_endDefinition", &ParallelMesh::endDefinition, R"(
Terminate the mesh creation (*for internal use*).
        )" )
        .def( "fix", &ParallelMesh::fix, R"(
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
        .def( "convertToLinear", &ParallelMesh::convertToLinear, R"(
Convert the mesh to a linear one.

Arguments:
    info (int) : verbosity mode (1 or 2). Default 1.

Returns:
    ParallelMesh: the linearized mesh.
        )",
              py::arg( "info" ) = 1 )
        .def( "convertToBiQuadratic", &ParallelMesh::convertToBiQuadratic, R"(
Convert the mesh to a bi-quadratic one.
For cells that have no bi-quadratic version, the quadratic version is used.

Arguments:
    info (int) : verbosity mode (1 or 2). Default 1.

Returns:
    ParallelMesh: the bi-quadratic mesh.
        )",
              py::arg( "info" ) = 1 )
        .def( "convertToQuadratic", &ParallelMesh::convertToQuadratic, R"(
Convert the mesh to a quadratic one.

Arguments:
    info (int) : verbosity mode (1 or 2). Default 1.

Returns:
    ParallelMesh: the quadratic mesh.
        )",
              py::arg( "info" ) = 1 )
        .def( "convertToCubic", &ParallelMesh::convertToCubic, R"(
Convert the mesh to a cubic one.

Arguments:
    info (int) : verbosity mode (1 or 2). Default 1.

Returns:
    Mesh: the cubic mesh.
        )",
              py::arg( "info" ) = 1 )
        .def( "getAllMedCellsTypes", &ParallelMesh::getAllMedCellsTypes, R"(
Return all Med types available in mesh (for all processors).

Returns:
    list[int]: List of Med types.)" )
        .def( "printMedFile", &ParallelMesh::printMedFile, R"(
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

#endif /* ASTER_HAVE_MPI */
