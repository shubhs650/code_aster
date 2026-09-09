/**
 * @file DOFNumberingInterface.cxx
 * @brief Interface python de DOFNumbering
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

#include "PythonBindings/EquationNumberingInterface.h"

#include "aster_pybind.h"

void exportEquationNumberingToPython( py::module_ &mod ) {

    py::class_< EquationNumbering, EquationNumberingPtr, DataStructure >( mod, "EquationNumbering" )
        .def( py::init( &initFactoryPtr< EquationNumbering > ) )
        .def( py::init( &initFactoryPtr< EquationNumbering, std::string > ) )
        .def( "getModel", &EquationNumbering::getModel )
        .def( "setModel", &EquationNumbering::setModel )
        .def( "getMesh", &EquationNumbering::getMesh )
        .def( "setMesh", &EquationNumbering::setMesh )
        // ---------------------------------------------------------------------
        .def( "useLagrangeDOF", &EquationNumbering::useLagrangeDOF, R"(
Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ---------------------------------------------------------------------
        .def( "useSingleLagrangeDOF", &EquationNumbering::useSingleLagrangeDOF,
              R"(
Single Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ---------------------------------------------------------------------
        .def( "getNumberOfDOFs", &EquationNumbering::getNumberOfDOFs,
              R"(
Returns the number of DOFs.

Arguments:
    local (bool): not used.

Returns:
    int: number of DOFs.
        )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getPhysicalQuantity", &EquationNumbering::getPhysicalQuantity, R"(
Returns the name of the physical quantity that is numbered.

Returns:
    str: physical quantity name.
        )" )
        // ---------------------------------------------------------------------
        .def( "getNoGhostDOFs", &EquationNumbering::getNoGhostDOFs,
              R"(
Returns the indexes of the DOFs owned locally (aka not ghost).

Returns:
int: indexes of the DOFs owned locally.
    )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "isParallel", &EquationNumbering::isParallel, R"(
The numbering is distributed across MPI processes for High Performance Computing.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ---------------------------------------------------------------------
        .def( "getNodeAndComponentFromDOF",
              py::overload_cast< const bool >( &EquationNumbering::getNodeAndComponentFromDOF,
                                               py::const_ ),
              R"(
            Return the list of node id and name of component for each dofs

            Arguments:
                local (bool) = True: if True use local node index else use global index in HPC

            Returns:
                list[tuple[int, str]] : node id and name of component for each dofs
            )",
              py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getNodeAndComponentFromDOF",
              py::overload_cast< const ASTERINTEGER, const bool >(
                  &EquationNumbering::getNodeAndComponentFromDOF, py::const_ ),
              R"(
            Return the node id and name of component for given DOF

            Arguments:
                dof (int): DOF index
                local (bool) = True: if True use local node index else use global index in HPC

            Returns:
                tuple[int, str] : node id and name of component
            )",
              py::arg( "dof" ), py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getNodeAndComponentIdFromDOF",
              py::overload_cast< const bool >( &EquationNumbering::getNodeAndComponentIdFromDOF,
                                               py::const_ ),
              R"(
            Return the list of node id and component id for each dofs

            Arguments:
                local (bool) = True: if True use local node index else use global index in HPC

            Returns:
                list[tuple[int, int]] : node id and component if for each dofs
            )",
              py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getNodeAndComponentIdFromDOF",
              py::overload_cast< const ASTERINTEGER, const bool >(
                  &EquationNumbering::getNodeAndComponentIdFromDOF, py::const_ ),
              R"(
            Return the node id and component id for given DOF

            Arguments:
                dof (int): DOF index
                local (bool) = True: if True use local node index else use global index in HPC

            Returns:
                tuple[int, int] : node id and component if for each dofs
            )",
              py::arg( "dof" ), py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getDOFFromNodeAndComponentId",
              py::overload_cast< const bool >( &EquationNumbering::getDOFFromNodeAndComponentId,
                                               py::const_ ),
              R"(
            Return the dict of dofs with the pair (node id, name id) as keys

            Arguments:
                local (bool) = True: if True use local DOF index else use global index in HPC

            Returns:
                dict[int, str] : dofs id for each node id and component id
            )",
              py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getDOFFromNodeAndComponent",
              py::overload_cast< const bool >( &EquationNumbering::getDOFFromNodeAndComponent,
                                               py::const_ ),
              R"(
           Return the dict of dofs with the pair (node id, component's name) as keys

            Arguments:
                local (bool) = True: if True use local dof index else use global index in HPC

            Returns:
                dict[int, str] : dofs id for each node id and component's name
            )",
              py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getComponents", &EquationNumbering::getComponents, R"(
            Get list of components

            Returns:
                list[str]: list of components
            )" )
        // ---------------------------------------------------------------------
        .def( "getComponentsIdToName", &EquationNumbering::getComponentsIdToName, R"(
            Get map between id and name of components

            Returns:
                dict[int]: map between id and name
            )" )
        // ---------------------------------------------------------------------
        .def( "getComponentsNameToId", &EquationNumbering::getComponentsNameToId, R"(
            Get map between id and name of components

            Returns:
                dict[int]: map between id and name
            )" )
        // ---------------------------------------------------------------------
        .def( "getDOFsWithDescription",
              py::overload_cast< const VectorString &, const VectorString &, const bool,
                                 const ASTERINTEGER >( &EquationNumbering::getDOFsWithDescription,
                                                       py::const_ ),
              R"(
            Get the dofs associated to the given component restricted to the given group.

            Arguments:
                cmps (list[str]): components to extract.
                groupNames (list[str]): group names to filter.
                local (bool): if True use local dof index else use global index in HPC.

            Returns:
                pair[list[int], list[str]]: list of nodes and list of components.
                list[int]: list of dofs.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupNames" ) = VectorString(),
              py::arg( "local" ) = true, py::arg( "same_rank" ) = PythonBool::None )
        // ---------------------------------------------------------------------
        .def( "getDOFsWithDescription",
              py::overload_cast< const VectorString &, const VectorLong &, const bool,
                                 const ASTERINTEGER >( &EquationNumbering::getDOFsWithDescription,
                                                       py::const_ ),
              R"(
            Get the dofs associated to the given component restricted to the given nodes.

            Arguments:
                cmps (list[str]): components to extract.
                nodes (list[int]): list of nodes to filter.
                local (bool): if True use local dof index else use global index in HPC.

            Returns:
                pair[list[int], list[str]]: list of nodes and list of components.
                list[int]: list of dofs.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "nodes" ) = VectorLong(),
              py::arg( "local" ) = true, py::arg( "same_rank" ) = PythonBool::None )
        .def( "getDOFs", &EquationNumbering::getDOFs,
              R"(
            Return list of DOFs

            Arguments:
                sameRank = False: Use only owned nodes / False: Use all nodes
                list_cmp = []: Use all cmp / keep only cmp given
                list_grpno = []: Use all nodes / keep only nodes given

            Returns:
                list[int]: list of dofs.
            )",
              py::arg( "sameRank" ) = false, py::arg( "list_cmp" ) = VectorString(),
              py::arg( "list_grpno" ) = VectorString() );
};
