/**
 * @file FieldOnNodesInterface.cxx
 * @brief Python interface for FieldOnNodes
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

#include "PythonBindings/FieldOnNodesInterface.h"

#include "aster_pybind.h"

#include "DataFields/FieldConverter.h"
#include "DataFields/FieldOnNodesBuilder.h"
#include "PythonBindings/DataStructureInterface.h"

void exportFieldOnNodesToPython( py::module_ &mod ) {
    /**
     * Object FieldOnNodesReal
     */

    py::class_< FieldOnNodesReal, FieldOnNodesRealPtr, DataField >( mod, "FieldOnNodesReal" )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, const FieldOnNodesReal & > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, BaseDOFNumberingPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, EquationNumberingPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesReal, const py::tuple & > ) )
        .def( define_pickling< FieldOnNodesReal >() )
        .def( py::init( []( const BaseMeshPtr mesh, const std::string &quantity,
                            const VectorString &cmps ) {
                  return FieldOnNodesPtrBuilder< ASTERDOUBLE >( mesh, quantity, cmps );
              } ),
              py::arg( "mesh" ), py::arg( "quantity" ), py::arg( "cmps" ) )
        .def( py::init( []( const BaseMeshPtr mesh, const std::string &quantity,
                            const std::map< std::string, ASTERDOUBLE > &values,
                            const VectorString &groupsOfNodes = {},
                            const VectorString &groupsOfCells = {} ) {
                  return FieldOnNodesPtrBuilder< ASTERDOUBLE >( mesh, quantity, values,
                                                                groupsOfNodes, groupsOfCells );
              } ),
              py::arg( "mesh" ), py::arg( "quantity" ), py::arg( "values" ),
              py::arg( "groupsOfNodes" ) = VectorString(),
              py::arg( "groupsOfCells" ) = VectorString() )
        .def( "copy", &FieldOnNodesReal::copy )
#ifdef ASTER_HAVE_MPI
        .def(
            "transfertToConnectionMesh",
            []( const FieldOnNodesRealPtr f, const ConnectionMeshPtr c ) {
                return transferToConnectionMesh( f, c );
            },
            R"(
Transfer SimpleFieldOnNodes to a ConnectionMesh

Returns:
    FieldOnNodesReal: transfered field
        )" )
        .def(
            "transferFromConnectionToParallelMesh",
            []( const FieldOnNodesRealPtr f, const BaseMeshPtr m ) {
                return transferFromConnectionToParallelMesh( f, m );
            },
            R"(
            Transfer FieldOnNodes from a ConnectionMesh to a ParallelMesh

            Arguments:
                mesh [Mesh]: the target mesh

            Returns:
                FieldOnNodesReal: transfered field
        )" )
#endif /* ASTER_HAVE_MPI */
        .def(
            "toSimpleFieldOnNodes",
            []( const FieldOnNodesReal &f ) { return toSimpleFieldOnNodes( f ); },
            R"(
Convert to SimpleFieldOnNodes

Returns:
    SimpleFieldOnNodesReal: field converted
        )" )
        .def(
            "toFieldOnCells",
            []( const FieldOnNodesReal &f, const FiniteElementDescriptorPtr fed,
                const std::string loc ) { return toFieldOnCells( f, fed, loc ); },
            R"(
            Converts to FieldOnCells

            Arguments:
                fed [FiniteElementDescriptor]: finite element descriptor
                loc [str] : name of localization like 'ELGA'.

            Returns:
                FieldOnCellsReal: field converted.
            )",
            py::arg( "fed" ), py::arg( "loc" ) )
        .def( "getPhysicalQuantity", &FieldOnNodesReal::getPhysicalQuantity )
        .def( "getMesh", &FieldOnNodesReal::getMesh )
        .def(
            "__getitem__",
            +[]( const FieldOnNodesReal &v, ASTERINTEGER i ) { return v.operator[]( i ); } )
        .def(
            "__setitem__", +[]( FieldOnNodesReal &v, ASTERINTEGER i,
                                ASTERDOUBLE f ) { return v.operator[]( i ) = f; } )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self / float() )
        .def( py::self /= float() )
        .def( py::self /= py::self )
        .def( -py::self )
        .def( "printMedFile", &FieldOnNodesReal::printMedFile, py::arg( "fileName" ),
              py::arg( "local" ) = true, py::arg( "version" ) = "" )
        .def( "setMesh", &FieldOnNodesReal::setMesh )
        .def( "setDescription", &FieldOnNodesReal::setDescription )
        .def( "build", &FieldOnNodesReal::build, py::arg( "mesh" ) = nullptr )
        .def( "getMesh", &FieldOnNodesReal::getMesh )
        .def( "getDescription", &FieldOnNodesReal::getDescription )
        .def( "getLocalization", &FieldOnNodesReal::getLocalization, R"(
            Get localization = NOEU

            Returns:
                str: "NOEU"
            )" )
        .def( "_restrict", &FieldOnNodesReal::restrict,
              R"(
            Return a new field restricted to the list of components and groups of nodes given

            Arguments:
                cmps[list[str]]: filter on list of components
                If empty, all components are used
                groupsOfNodes[list[str]]: filter on list of groups of nodes (default=" ").
                If empty, the full mesh is used
                same_rank : - None: keep all nodes (default: None)
                            - True: keep the nodes which are owned by the current MPI-rank
                            - False: keep the nodes which are not owned by the current MPI-rank

            Returns:
                FieldOnNodesReal: field restricted.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupsOfNodes" ) = VectorString(),
              py::arg( "same_rank" ) = PythonBool::None )
        .def( "asPhysicalQuantity", &FieldOnNodesReal::asPhysicalQuantity,
              R"(
            Return a new field with a new physical quantity and renamed components.

            Arguments:
                physQuantity [str]: name of the new physical quantity
                map_cmps[dict[str, str]]: dict to rename components
                (only renamed component will be keeped)

            Returns:
                FieldOnNodesReal: field with name physical quantity.
            )",
              py::arg( "physQuantity" ), py::arg( "map_cmps" ) )
        .def(
            "getRealPart", []( const FieldOnNodesReal &f ) { return getRealPart( f ); },
            R"(
Extract the real part of the real field (the field is duplicated)

Returns:
    FieldOnNodesReal: real part
        )" )
        .def(
            "getImaginaryPart", []( const FieldOnNodesReal &f ) { return getImaginaryPart( f ); },
            R"(
Extract the imaginary part of the real field (a 0-filled field is produced)

Returns:
    FieldOnNodesReal: imaginary part
        )" )
        .def( "copyUsingDescription", &FieldOnNodesReal::copyUsingDescription,
              R"(
            Return a new field using the description.
            Be careful, Lagrange DOFs are set to zero. Moreover, components that are
            not present in the field are also set to zero in the output field.

            Arguments:
                desc [EquationNumbering]: description of equations
                warn [bool]: If set to true, raises a warning if values are set to zero

            Returns:
                FieldOnNodesReal: field using new description.
            )",
              py::arg( "desc" ), py::arg( "warn" ) = true )
        .def( "updateValuePointers", &FieldOnNodesReal::updateValuePointers )
        .def( "getComponents", &FieldOnNodesReal::getComponents, R"(
            Get list of components

            Returns:
                list[str]: list of components
            )" )
        .def( "getNumberOfComponents", &FieldOnNodesReal::getNumberOfComponents, R"(
            Get number of components

            Returns:
                int: number of components
            )" )
        .def( "norm", &FieldOnNodesReal::norm, R"(
            Return the euclidean norm of the field

            Arguments:
                normType (str): "NORM_1", "NORM_2", "NORM_INFINITY" (default: "NORM_INFINITY")
                list_cmp (list[str]) : list of components used to compute norm (default: all)

            Returns:
                float: euclidean norm
            )",
              py::arg( "normType" ) = "NORM_INFINITY", py::arg( "list_cmp" ) = VectorString() )
        .def( "dot", &FieldOnNodesReal::dot, R"(
            Return the dot product of two fields

            Arguments:
                other (FieldOnNodes): other field

            Returns:
                float: dot product
            )",
              py::arg( "other" ) )
        .def( "size", &FieldOnNodesReal::size, R"(
            Return the size of the field

            Returns:
                int: number of element in the field
            )" )
        .def( "transform", &FieldOnNodesReal::transform, R"(
            Apply a function to each value of the object.

            Arguments:
                func (*callable*): Callable Python object

            Returns:
                FieldOnNodesReal: New FieldOnNodes object with the transformed values
            )",
              py::arg( "func" ) )
        .def( "scale", &FieldOnNodesReal::scale, R"(
            Scale in-place the field by a diagonal matrix stored as an array

            Arguments:
                vect (float): diagonal matrix stored as an array
            )",
              py::arg( "vect" ) )
        .def( "applyLagrangeScaling", &FieldOnNodesReal::applyLagrangeScaling, R"(
            Multiply in-place the Lagrange multipliers DOFs by the scaling value

            Arguments:
                scaling (float): scaling velue
            )",
              py::arg( "scaling" ) )
#ifdef ASTER_HAVE_PETSC
        .def( "fromPetsc", &FieldOnNodesReal::fromPetsc,
              R"(
            Import a PETSc vector into the field.

            Arguments:
                vec (Vec): The PETSc vector
                scaling (float) : The scaling of the Lagrange DOFs
                local (bool) : Only import the dof that are local to the subdomain
            )",
              py::arg( "vec" ), py::arg( "scaling" ) = 1.0, py::arg( "local" ) = false )
#endif
        .def( "updateGhostValues", &FieldOnNodesReal::updateGhostValues,
              R"(
            Communicates the values of the ghost DOFs on a FieldOnNodes.
  )" )
        .def( "setValues", py::overload_cast< const ASTERDOUBLE & >( &FieldOnNodesReal::setValues ),
              R"(
            Set values of the field

            Arguments:
                value (float): value to set
            )",
              py::arg( "value" ) )
        .def( "setValues", py::overload_cast< const VectorReal & >( &FieldOnNodesReal::setValues ),
              R"(
            Set values of the field

            Arguments:
                values (list[float]): list of values to set
            )",
              py::arg( "values" ) )
        .def( "setValues",
              py::overload_cast< const std::map< std::string, ASTERDOUBLE > &, VectorString >(
                  &FieldOnNodesReal::setValues ),
              R"(
            Set values of the field where components and values are given as a dict.
            If the component is not present in the field then it is discarded
            Example: { "X1" : 0.0, "X3" : 0.0 }

            Arguments:
                value (dict[str, float]): dict of values to set (key: str, value: float)
                groupsOfNodes (list[str]): list of groups. If empty, the full mesh is considered
            )",
              py::arg( "value" ), py::arg( "groupsOfNodes" ) = VectorString() )
        .def( "getValues", py::overload_cast<>( &FieldOnNodesReal::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Returns:
                list[float]: List of values.
            )" )
        .def( "getValues",
              py::overload_cast< const VectorString &, const VectorString & >(
                  &FieldOnNodesReal::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Arguments:
                cmps[list[str]]: filter on list of components
                groupsOfNodes[list[str]]: filter on list of groups of nodes (default=" ").
                If empty, the full mesh is used

            Returns:
                list[double]: List of values.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupsOfNodes" ) = VectorString() )
        .def( "getValues",
              py::overload_cast< const VectorLong & >( &FieldOnNodesReal::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

            Arguments:
                dofs: dofs to extract

            Returns:
                list[double]: List of values.
            )",
              py::arg( "dofs" ) = VectorLong() );
    /**
     * Object FieldOnNodesComplex
     */
    py::class_< FieldOnNodesComplex, FieldOnNodesComplexPtr, DataField >( mod,
                                                                          "FieldOnNodesComplex" )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex > ) )
        .def( py::init< const FieldOnNodesComplex & >() )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex, BaseDOFNumberingPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex, EquationNumberingPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesComplex, const py::tuple & > ) )
        .def( define_pickling< FieldOnNodesComplex >() )
        .def(
            "toSimpleFieldOnNodes",
            []( const FieldOnNodesComplex &f ) { return toSimpleFieldOnNodes( f ); },
            R"(
Convert to SimpleFieldOnNodes

Returns:
    SimpleFieldOnNodesComplex: field converted
        )" )
        .def(
            "getRealPart", []( const FieldOnNodesComplex &f ) { return getRealPart( f ); },
            R"(
Extract the real part of the complex field

Returns:
    FieldOnNodesReal: real part
        )" )
        .def(
            "getImaginaryPart",
            []( const FieldOnNodesComplex &f ) { return getImaginaryPart( f ); },
            R"(
Extract the imaginary part of the complex field

Returns:
    FieldOnNodesReal: imaginary part
        )" )
        .def( "getPhysicalQuantity", &FieldOnNodesComplex::getPhysicalQuantity )
        .def( "getMesh", &FieldOnNodesComplex::getMesh )
        .def(
            "__getitem__",
            +[]( const FieldOnNodesComplex &v, ASTERINTEGER i ) { return v.operator[]( i ); } )
        .def(
            "__setitem__", +[]( FieldOnNodesComplex &v, ASTERINTEGER i,
                                ASTERCOMPLEX f ) { return v.operator[]( i ) = f; } )
        .def( "printMedFile", &FieldOnNodesComplex::printMedFile, py::arg( "fileName" ),
              py::arg( "local" ) = true, py::arg( "version" ) = "" )
        .def( "setMesh", &FieldOnNodesComplex::setMesh )
        .def( "setDescription", &FieldOnNodesComplex::setDescription )
        .def( "build", &FieldOnNodesComplex::build, py::arg( "mesh" ) = nullptr )
        .def( "getMesh", &FieldOnNodesComplex::getMesh )
        .def( "getDescription", &FieldOnNodesComplex::getDescription )
        .def( "getLocalization", &FieldOnNodesComplex::getLocalization, R"(
            Get localization = NOEU

            Returns:
                str: "NOEU"
            )" )
        .def( "_restrict", &FieldOnNodesComplex::restrict,
              R"(
            Return a new field restricted to the list of components and groups of nodes given

            Arguments:
                cmps[list[str]]: filter on list of components
                If empty, all components are used
                groupsOfNodes[list[str]]: filter on list of groups of nodes (default=" ").
                If empty, the full mesh is used
                same_rank : - None: keep all nodes (default: None)
                            - True: keep the nodes which are owned by the current MPI-rank
                            - False: keep the nodes which are not owned by the current MPI-rank

            Returns:
                FieldOnNodesComplex: field restricted.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupsOfNodes" ) = VectorString(),
              py::arg( "same_rank" ) = PythonBool::None )
        .def( "getValues", py::overload_cast<>( &FieldOnNodesComplex::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Returns:
                list[complex]: List of values.
            )" )
        .def( "getValues",
              py::overload_cast< const VectorString &, const VectorString & >(
                  &FieldOnNodesComplex::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Arguments:
                cmps[list[str]]: filter on list of components
                groupsOfNodes[list[str]]: filter on list of groups of nodes (default=" ").
                If empty, the full mesh is used

            Returns:
                list[complex]: List of values.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupsOfNodes" ) = VectorString() )
        .def(
            "getValues",
            py::overload_cast< const VectorLong & >( &FieldOnNodesComplex::getValues, py::const_ ),
            R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

            Arguments:
                dofs: dofs to extract

            Returns:
                list[complex]: List of values.
            )",
            py::arg( "dofs" ) = VectorLong() )
        .def( "getComponents", &FieldOnNodesComplex::getComponents, R"(
            Get list of components

            Returns:
                list[str]: list of components
            )" )
        .def( "getNumberOfComponents", &FieldOnNodesComplex::getNumberOfComponents, R"(
            Get number of components

            Returns:
                int: number of components
            )" )
        .def( "transform", &FieldOnNodesComplex::transform, R"(
            Apply a function to each value of the object.

            Arguments:
                func (*callable*): Callable Python object

            Returns:
                FieldOnNodesComplex: New FieldOnNodes object with the transformed values
            )",
              py::arg( "func" ) )
        .def( "scale", &FieldOnNodesComplex::scale, R"(
            Scale in-place the field by a diagonal matrix stored as an array

            Arguments:
                vect (float): diagonal matrix stored as an array
            )",
              py::arg( "vect" ) )
        .def( "dot", &FieldOnNodesComplex::dot, R"(
            Return the dot product of two complex fields

            Arguments:
                other (FieldOnNodes): other field

            Returns:
                complex: dot product
            )",
              py::arg( "other" ) )
        .def( "norm", &FieldOnNodesComplex::norm, R"(
            Return the euclidean norm of the field

            Arguments:
                normType (str): "NORM_1", "NORM_2", "NORM_INFINITY" (default: "NORM_INFINITY")
                list_cmp (list[str]) : list of components used to compute norm (default: all)

            Returns:
                float: euclidean norm
            )",
              py::arg( "normType" ) = "NORM_INFINITY", py::arg( "list_cmp" ) = VectorString() )
        .def( "setValues",
              py::overload_cast< const ASTERCOMPLEX & >( &FieldOnNodesComplex::setValues ), R"(
            Set values of the field

            Arguments:
                value (complex): value to set
            )",
              py::arg( "value" ) )
        .def( "setValues",
              py::overload_cast< const VectorComplex & >( &FieldOnNodesComplex::setValues ),
              R"(
            Set values of the field

            Arguments:
                values (list[complex]): list of values to set
            )",
              py::arg( "values" ) )
        .def( "updateValuePointers", &FieldOnNodesComplex::updateValuePointers );

    /**
     * Object FieldOnNodesLong
     */
    py::class_< FieldOnNodesLong, FieldOnNodesLongPtr, DataField >( mod, "FieldOnNodesLong" )
        .def( py::init( &initFactoryPtr< FieldOnNodesLong > ) )
        .def( py::init< const FieldOnNodesLong & >() )
        .def( py::init( &initFactoryPtr< FieldOnNodesLong, BaseDOFNumberingPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesLong, const py::tuple & > ) )
        .def( define_pickling< FieldOnNodesLong >() )
        .def( "build", &FieldOnNodesLong::build, py::arg( "mesh" ) = nullptr )
        .def( "setDescription", &FieldOnNodesLong::setDescription )
        .def( "getDescription", &FieldOnNodesLong::getDescription )
        .def( "getMesh", &FieldOnNodesLong::getMesh )
        .def( "setMesh", &FieldOnNodesLong::setMesh );

    /**
     * Object FieldOnNodesChar8
     */
    py::class_< FieldOnNodesChar8, FieldOnNodesChar8Ptr, DataField >( mod, "FieldOnNodesChar8" )
        .def( py::init( &initFactoryPtr< FieldOnNodesChar8 > ) )
        .def( py::init< const FieldOnNodesChar8 & >() )
        .def( py::init( &initFactoryPtr< FieldOnNodesChar8, BaseDOFNumberingPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesChar8, const py::tuple & > ) )
        .def( define_pickling< FieldOnNodesChar8 >() )
        .def( "build", &FieldOnNodesChar8::build, py::arg( "mesh" ) = nullptr )
        .def( "setDescription", &FieldOnNodesChar8::setDescription )
        .def( "getDescription", &FieldOnNodesChar8::getDescription )
        .def( "getMesh", &FieldOnNodesChar8::getMesh )
        .def( "setMesh", &FieldOnNodesChar8::setMesh );
};
