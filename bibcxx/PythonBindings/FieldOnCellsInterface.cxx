/**
 * @file FieldOnCellsInterface.cxx
 * @brief Python interface for FieldOnCells
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

#include "PythonBindings/FieldOnCellsInterface.h"

#include "aster_pybind.h"

#include "DataFields/FieldConverter.h"
#include "DataFields/FieldOnCellsBuilder.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "PythonBindings/DataStructureInterface.h"

void exportFieldOnCellsToPython( py::module_ &mod ) {
    py::class_< FieldOnCellsReal, FieldOnCellsRealPtr, DataField >( mod, "FieldOnCellsReal" )
        .def( py::init( &initFactoryPtr< FieldOnCellsReal > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsReal, ModelPtr, std::string, std::string > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsReal, FiniteElementDescriptorPtr, std::string,
                                         std::string > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsReal, const FieldOnCellsReal & > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsReal, const py::tuple & > ) )
        .def( define_pickling< FieldOnCellsReal >() )
        .def( py::init( []( const ModelPtr model, const std::string &loc,
                            const std::string &quantity, const BehaviourPropertyPtr behaviour,
                            const ElementaryCharacteristicsPtr carael ) {
                  return FieldOnCellsPtrBuilder< ASTERDOUBLE >( model, loc, quantity, behaviour,
                                                                carael );
              } ),
              py::arg( "model" ), py::arg( "loc" ), py::arg( "quantity" ), py::arg( "behaviour" ),
              py::arg( "elem_char" ) )
        .def( py::init( []( const ModelPtr model, const std::string &loc,
                            const std::string &quantity, const BehaviourPropertyPtr behaviour ) {
                  return FieldOnCellsPtrBuilder< ASTERDOUBLE >( model, loc, quantity, behaviour,
                                                                nullptr );
              } ),
              py::arg( "model" ), py::arg( "loc" ), py::arg( "quantity" ), py::arg( "behaviour" ) )
        .def(
            py::init( []( const ModelPtr model, const std::string &loc, const std::string &quantity,
                          const ElementaryCharacteristicsPtr carael ) {
                return FieldOnCellsPtrBuilder< ASTERDOUBLE >( model, loc, quantity, nullptr,
                                                              carael );
            } ),
            py::arg( "model" ), py::arg( "loc" ), py::arg( "quantity" ), py::arg( "elem_char" ) )
        .def( "copy", &FieldOnCellsReal::copy, R"(
            Return a duplicated FieldOnCellsReal as a copy

            Returns:
                FieldOnCellsReal
            )" )
        .def( "exists", &FieldOnCellsReal::exists, R"(
            The field exists or not ?

            Returns:
                Bool
            )" )
        .def(
            "toSimpleFieldOnCells",
            []( const FieldOnCellsReal &f ) { return toSimpleFieldOnCells( f ); },
            R"(
Convert to SimpleFieldOnNodes

Returns:
    SimpleFieldOnNodesReal: field converted
        )" )
        .def(
            "toFieldOnNodes", []( const FieldOnCellsReal &f ) { return toFieldOnNodes( f ); },
            R"(
Convert to FieldOnNodes

Returns:
    FieldOnNodesReal: field converted
        )" )
        .def(
            "toSimpleFieldOnNodes",
            []( const FieldOnCellsReal &f ) { return toSimpleFieldOnNodes( f ); },
            R"(
Convert to SimpleFieldOnNodes

Returns:
    SimpleFieldOnNodesReal: field converted
        )" )
        .def( "asLocalization", &FieldOnCellsReal::asLocalization,
              R"(
            Return a new field interpolated at the given localozation.

            Arguments:
                loc [str]: name of localization (ELEM, ELNO or ELGA)

            Returns:
                FieldOnCellsReal: new field with new localization.
            )",
              py::arg( "loc" ) )
        .def( "getDescription", &FieldOnCellsReal::getDescription, R"(
            Return the descriptor associated with the FieldOnCellsReal object

            Returns:
                FiniteElementDescriptor: FiniteElementDescriptor Object
            )" )
        .def( "getMesh", &FieldOnCellsReal::getMesh, R"(
            Return the Mesh associated with the FieldOnCellsReal object

            Returns:
                BaseMesh: Mesh object
            )" )
        .def( "setDescription", &FieldOnCellsReal::setDescription )
        .def( "getDescription", &FieldOnCellsReal::getDescription )
        .def( "build", &FieldOnCellsReal::build,
              py::arg( "feds" ) = std::vector< FiniteElementDescriptorPtr >() )
        .def(
            "__getitem__", +[]( const FieldOnCellsReal &v, ASTERINTEGER i ) { return v[i]; } )
        .def(
            "__setitem__",
            +[]( FieldOnCellsReal &v, ASTERINTEGER i, float f ) { return v.operator[]( i ) = f; } )
        .def(
            "__len__", +[]( const FieldOnCellsReal &v ) { return v.size(); } )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( -py::self )
        .def( "setValues", py::overload_cast< const ASTERDOUBLE & >( &FieldOnCellsReal::setValues ),
              R"(
            Set values of the field

            Arguments:
                value (float): value to set
            )",
              py::arg( "value" ) )
        .def( "setValues", py::overload_cast< const VectorReal & >( &FieldOnCellsReal::setValues ),
              R"(
            Set values of the field

            Arguments:
                values (list[float]): list of values to set
            )",
              py::arg( "values" ) )
        .def( "getValues", py::overload_cast<>( &FieldOnCellsReal::getValues, py::const_ ), R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Returns:
                list[float]: List of values.
            )" )
        .def( "getValues",
              py::overload_cast< const VectorLong & >( &FieldOnCellsReal::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

            Arguments:
                dofs: dofs to extract

            Returns:
                list[float]: List of values.
            )",
              py::arg( "dofs" ) = VectorLong() )
        .def( "size", &FieldOnCellsReal::size, R"(
            Return the size of the field

            Returns:
                int: number of element in the field
            )" )
        .def( "transform", &FieldOnCellsReal::transform, R"(
            Apply a function to each value of the object.

            Arguments:
                func (*callable*): Callable Python object

            Returns:
                FieldOnCellsReal: New FieldOnCells object with the transformed values
            )",
              py::arg( "func" ) )
        .def( "getPhysicalQuantity", &FieldOnCellsReal::getPhysicalQuantity, R"(
            Get physical quantity

            Returns:
                str: physical quantity
            )" )
        .def( "getLocalization", &FieldOnCellsReal::getLocalization, R"(
            Get localization between ELEM, ELNO and ELGA

            Returns:
                str: localization
            )" )
        .def( "getComponents", &FieldOnCellsReal::getComponents, R"(
            Get list of components

            Returns:
                list[str]: list of components
            )" )
        .def( "getNumberOfComponents", &FieldOnCellsReal::getNumberOfComponents, R"(
            Get number of components

            Returns:
                int: number of components
            )" )
        .def( "printMedFile", &FieldOnCellsReal::printMedFile, R"(
            Print the field in MED format.

            Arguments:
                filename (Path|str): Path to the file to be printed.
                local (bool): Print local values only (relevant for ParallelMesh only,
                    default: *True*)
                version (str): Version of MED file.

            Returns:
                bool: *True* if succeeds, *False* otherwise.
            )",
              py::arg( "filename" ), py::arg( "local" ) = true, py::arg( "version" ) = "" )
        .def( "norm", &FieldOnCellsReal::norm, R"(
            Return the euclidean norm of the field

            Arguments:
                normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"

            Returns:
                float: euclidean norm
            )" )
        .def( "dot", &FieldOnCellsReal::dot, R"(
            Return the dot product of two fields

            Arguments:
                field (FieldOnCells): other field

            Returns:
                float: dot product
            )",
              py::arg( "other" ) )
        .def( "checkInternalStateVariables", &FieldOnCellsReal::checkInternalStateVariables, R"(
            Check consistency of internal states variables with behaviour.
            If you give previous behaviour, check is more precise (name of beahviour for instance)

            Arguments:
                prevBehaviour (ConstantFieldOnCellsChar16): previous behaviour
                currBehaviour (ConstantFieldOnCellsChar16): current behaviour
                newFEDesc (FiniteElementDescriptorPtr): new finite element descriptor


            )",
              py::arg( "prevBehaviour" ), py::arg( "currBehaviour" ), py::arg( "newFEDesc" ) )
        .def( "compareShape", &FieldOnCellsReal::compareShape, R"(
            Compare structure of field with another one and project on new model if require

            Arguments:
                fieldModel (FieldOnCellsRealPtr): field as model
                projectOnLigrel (bool) : project field on new model (from model field)
                paraName (string) : name of parameter to complete the new values in field

            Returns:
                iret (integer) : error code

            )",
              py::arg( "fieldModel" ), py::arg( "projectOnLigrel" ), py::arg( "paraName" ) );

    py::class_< FieldOnCellsComplex, FieldOnCellsComplexPtr, DataField >( mod,
                                                                          "FieldOnCellsComplex" )
        .def( py::init( &initFactoryPtr< FieldOnCellsComplex > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsComplex, FiniteElementDescriptorPtr,
                                         std::string, std::string > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsComplex, const py::tuple & > ) )
        .def( py::init< const FieldOnCellsComplex & >() )
        .def( define_pickling< FieldOnCellsComplex >() )
        .def( "copy", &FieldOnCellsComplex::copy )
        .def( "setDescription", &FieldOnCellsComplex::setDescription )
        .def( "getDescription", &FieldOnCellsComplex::getDescription )
        .def( "getMesh", &FieldOnCellsComplex::getMesh, R"(
            Return the Mesh associated with the FieldOnCellsReal object

            Returns:
                BaseMesh: Mesh object
            )" )
        .def( "build", &FieldOnCellsComplex::build,
              py::arg( "feds" ) = std::vector< FiniteElementDescriptorPtr >() )
        .def( "setValues",
              py::overload_cast< const ASTERCOMPLEX & >( &FieldOnCellsComplex::setValues ),
              R"(
            Set values of the field

            Arguments:
                value (complex): value to set
            )",
              py::arg( "value" ) )
        .def( "setValues",
              py::overload_cast< const VectorComplex & >( &FieldOnCellsComplex::setValues ),
              R"(
            Set values of the field

            Arguments:
                values (list[complex]): list of values to set
            )",
              py::arg( "values" ) )
        .def( "getValues", py::overload_cast<>( &FieldOnCellsComplex::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Returns:
                list[complex]: List of values.
            )" )
        .def(
            "getValues",
            py::overload_cast< const VectorLong & >( &FieldOnCellsComplex::getValues, py::const_ ),
            R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

            Arguments:
                dofs: dofs to extract

            Returns:
                list[complex]: List of values.
            )",
            py::arg( "dofs" ) = VectorLong() )
        .def(
            "__getitem__", +[]( const FieldOnCellsComplex &v, int i ) { return v[i]; } )
        .def(
            "__setitem__", +[]( FieldOnCellsComplex &v, ASTERINTEGER i,
                                ASTERCOMPLEX f ) { return v.operator[]( i ) = f; } )
        .def(
            "__len__", +[]( const FieldOnCellsComplex &v ) { return v.size(); } )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def(
            "toFieldOnNodes", []( const FieldOnCellsComplex &f ) { return toFieldOnNodes( f ); },
            R"(
Convert to FieldOnNodes

Returns:
    FieldOnCellsComplex: field converted
        )" )
        .def( "getPhysicalQuantity", &FieldOnCellsComplex::getPhysicalQuantity, R"(
            Get physical quantity

            Returns:
                str: physical quantity
            )" )
        .def( "getLocalization", &FieldOnCellsComplex::getLocalization, R"(
            Get localization between ELEM, ELNO and ELGA

            Returns:
                str: localization
            )" )
        .def( "size", &FieldOnCellsComplex::size, R"(
            Return the size of the field

            Returns:
                int: number of element in the field
            )" )
        .def( "transform", &FieldOnCellsComplex::transform, R"(
            Apply a function to each value of the object.

            Arguments:
                func (*callable*): Callable Python object

            Returns:
                FieldOnCellsComplex: New FieldOnCells object with the transformed values
            )",
              py::arg( "func" ) )
        .def( "printMedFile", &FieldOnCellsComplex::printMedFile, R"(
            Print the field in MED format.

            Arguments:
                filename (Path|str): Path to the file to be printed.
                local (bool): Print local values only (relevant for ParallelMesh only,
                    default: *True*)
                version (str): Version of MED file.

            Returns:
                bool: *True* if succeeds, *False* otherwise.
            )",
              py::arg( "filename" ), py::arg( "local" ) = true, py::arg( "version" ) = "" );

    py::class_< FieldOnCellsLong, FieldOnCellsLongPtr, DataField >( mod, "FieldOnCellsLong" )
        .def( py::init( &initFactoryPtr< FieldOnCellsLong > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsLong, const py::tuple & > ) )
        .def( py::init< const FieldOnCellsLong & >() )
        .def( define_pickling< FieldOnCellsLong >() )
        .def( "copy", &FieldOnCellsLong::copy )
        .def( "setDescription", &FieldOnCellsLong::setDescription )
        .def( "getDescription", &FieldOnCellsLong::getDescription )
        .def( "getMesh", &FieldOnCellsLong::getMesh, R"(
            Return the Mesh associated with the FieldOnCellsReal object

            Returns:
                BaseMesh: Mesh object
            )" )
        .def( "build", &FieldOnCellsLong::build,
              py::arg( "feds" ) = std::vector< FiniteElementDescriptorPtr >() )
        .def( "setValues",
              py::overload_cast< const ASTERINTEGER & >( &FieldOnCellsLong::setValues ),
              R"(
            Set values of the field

            Arguments:
                value (complex): value to set
            )",
              py::arg( "value" ) )
        .def( "setValues", py::overload_cast< const VectorLong & >( &FieldOnCellsLong::setValues ),
              R"(
            Set values of the field

            Arguments:
                values (list[complex]): list of values to set
            )",
              py::arg( "values" ) )
        .def( "getValues", py::overload_cast<>( &FieldOnCellsLong::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Returns:
                list[int]: List of values.
            )" )
        .def( "getValues",
              py::overload_cast< const VectorLong & >( &FieldOnCellsLong::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...) corresponding to list of dofs

            Arguments:
                dofs: dofs to extract

            Returns:
                list[int]: List of values.
            )",
              py::arg( "dofs" ) = VectorLong() )
        .def(
            "__getitem__", +[]( const FieldOnCellsLong &v, int i ) { return v[i]; } )
        .def(
            "__setitem__", +[]( FieldOnCellsLong &v, ASTERINTEGER i,
                                ASTERINTEGER f ) { return v.operator[]( i ) = f; } )
        .def(
            "__len__", +[]( const FieldOnCellsLong &v ) { return v.size(); } )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self * float() )
        .def( float() * py::self )
        .def( "size", &FieldOnCellsLong::size, R"(
            Return the size of the field

            Returns:
                int: number of element in the field
            )" )
        .def( "printMedFile", &FieldOnCellsLong::printMedFile, R"(
            Print the field in MED format.

            Arguments:
                filename (Path|str): Path to the file to be printed.
                local (bool): Print local values only (relevant for ParallelMesh only,
                    default: *True*)
                version (str): Version of MED file.

            Returns:
                bool: *True* if succeeds, *False* otherwise.
            )",
              py::arg( "filename" ), py::arg( "local" ) = true, py::arg( "version" ) = "" );
    /**
     * Object FieldOnCellsChar8
     */
    py::class_< FieldOnCellsChar8, FieldOnCellsChar8Ptr, DataField >( mod, "FieldOnCellsChar8" )
        .def( py::init( &initFactoryPtr< FieldOnCellsChar8 > ) )
        .def( py::init( &initFactoryPtr< FieldOnCellsChar8, const py::tuple & > ) )
        .def( py::init< const FieldOnCellsChar8 & >() )
        .def( define_pickling< FieldOnCellsChar8 >() )
        .def( "setDescription", &FieldOnCellsChar8::setDescription )
        .def( "getDescription", &FieldOnCellsChar8::getDescription, R"(
            Return the description associated with the FieldOnCellsChar8 object

            Returns:
                FiniteElementDescriptor: FiniteElementDescriptor Object
            )" )
        .def( "getMesh", &FieldOnCellsChar8::getMesh, R"(
            Return the Mesh associated with the FieldOnCellsChar8 object

            Returns:
                BaseMesh: Mesh object
            )" )
        .def( "build", &FieldOnCellsChar8::build,
              py::arg( "feds" ) = std::vector< FiniteElementDescriptorPtr >() );
};
