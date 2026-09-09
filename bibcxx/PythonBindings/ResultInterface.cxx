/**
 * @file ResultInterface.cxx
 * @brief Interface python de Result
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

#include "PythonBindings/ResultInterface.h"

#include "aster_pybind.h"

void exportResultToPython( py::module_ &mod ) {

    py::class_< Result, Result::ResultPtr, DataStructure >( mod, "Result" )
        .def( py::init( &initFactoryPtr< Result, std::string > ) )
        .def( py::init( &initFactoryPtr< Result, std::string, std::string > ) )
        .def( "allocate", &Result::allocate, R"(
Allocate result

Arguments:
    nb_index (int):  number of index to allocate
        )",
              py::arg( "nb_index" ) )
        .def( "exists", &Result::exists, R"(
The result exists or nor

Returns:
    bool: True if the result exists else False
        )" )
        .def( "clear", py::overload_cast<>( &Result::clear ),
              R"(
Clear fields, models, parameters, ... in result
)" )
        .def( "clear", py::overload_cast< const ASTERINTEGER & >( &Result::clear ),
              R"(
Clear fields, models, parameters, ... in result from the given index

Arguments:
    index (int): index from begin cleaning
        )",
              py::arg( "index" ) )
        .def( "setTime", &Result::setTime, R"(
Add time at the specified index

Arguments:
    time (float): time value to save
    index (int):  index where to save time value
        )",
              py::arg( "time" ), py::arg( "index" ) )
        .def( "setParameterValue",
              py::overload_cast< std::string, ASTERDOUBLE, ASTERINTEGER >(
                  &Result::setParameterValue ),

              R"(
Add parameter at the specified index

Arguments:
    para_name (float): parameter name to store
    para_value (float): parameter value to store
    index (int):  index where to save value of parameter
        )",
              py::arg( "para_name" ), py::arg( "para_value" ), py::arg( "index" ) )
        .def( "setParameterValue",
              py::overload_cast< std::string, std::string, ASTERINTEGER >(
                  &Result::setParameterValue ),
              R"(
Add parameter at the specified index

Arguments:
    para_name (float): parameter name to store
    para_value (str): parameter value to store
    index (int):  index where to save value of parameter
        )",
              py::arg( "para_name" ), py::arg( "para_value" ), py::arg( "index" ) )
        .def( "getTime", &Result::getTime, R"(
Get time at the specified index

Arguments:
    index (int):  index where to save time value

Returns
    float: time value
        )",
              py::arg( "index" ) )
        .def( "addEquationNumbering", &Result::addEquationNumbering )
        .def( "setMaterialField",
              py::overload_cast< const MaterialFieldPtr &, bool >( &Result::setMaterialField ), R"(
Set material field on all indexs

Arguments:
    mater (MaterialField): material field to set.
    exists_ok (bool): If *True*, pass silently if a Model is already defined. *False* by default.
        )",
              py::arg( "mater" ), py::arg( "exists_ok" ) = false )
        .def( "setMaterialField",
              py::overload_cast< const MaterialFieldPtr &, ASTERINTEGER >(
                  &Result::setMaterialField ),
              R"(
Set material field on the specified index

Arguments:
    mater (MaterialField): material field to set.
    index (int): index to set
        )",
              py::arg( "mater" ), py::arg( "index" ) )
        .def( "setListOfLoads", &Result::setListOfLoads, R"(
Set list of loads on the specified index

Arguments:
    load (ListOfLoads): list of loads to set.
    index (int): index to set
        )",
              py::arg( "load" ), py::arg( "index" ) )
        .def( "setModel", py::overload_cast< const ModelPtr &, bool >( &Result::setModel ), R"(
Set model on all indexs

Arguments:
    model (Model): Model to be assigned.
    exists_ok (bool): If *True*, pass silently if a Model is already defined. *False* by default.
        )",
              py::arg( "model" ), py::arg( "exists_ok" ) = false )
        .def( "setModel", py::overload_cast< const ModelPtr &, ASTERINTEGER >( &Result::setModel ),
              R"(
Set model on the specified index

Arguments:
    model (Model): model to set
    index (int): index to set
        )",
              py::arg( "model" ), py::arg( "index" ) )
        .def( "setElementaryCharacteristics",
              py::overload_cast< const ElementaryCharacteristicsPtr &, bool >(
                  &Result::setElementaryCharacteristics ),
              R"(
Set elementary characterictics on all indexs

Arguments:
    cara_elem (ElementaryCharacteristics): elementary characterictics to set.
    exists_ok (bool): If *True*, pass silently if a Model is already defined. *False* by default.
        )",
              py::arg( "cara_elem" ), py::arg( "exists_ok" ) = false )
        .def( "setElementaryCharacteristics",
              py::overload_cast< const ElementaryCharacteristicsPtr &, ASTERINTEGER >(
                  &Result::setElementaryCharacteristics ),
              R"(
Set elementary characterictics on the specified index

Arguments:
    cara_elem (ElementaryCharacteristics): elementary characterictics to set.
    index (int): index to set
        )",
              py::arg( "cara_elem" ), py::arg( "index" ) )
        .def( "printListOfFields", &Result::printListOfFields, R"(
Print the names of all fields (real, complex, ...) stored in the result.
        )" )
        .def( "getAllElementaryCharacteristics", &Result::getAllElementaryCharacteristics, R"(
Return the list of all elementary characteristics used in the result

Returns:
    list[ElementaryCharacteristics]: list of ElementaryCharacteristics.
        )" )
        .def(
            "getElementaryCharacteristics",
            py::overload_cast< ASTERINTEGER >( &Result::getElementaryCharacteristics, py::const_ ),
            R"(
Get elementary characterictics at the specfied index

Arguments:
    index (int): index

Returns:
    ElementaryCharacteristics: a pointer to elementary characterictics.
        )",
            py::arg( "index" ) )
        .def( "getElementaryCharacteristics",
              py::overload_cast<>( &Result::getElementaryCharacteristics, py::const_ ) )
        .def(
            "hasElementaryCharacteristics",
            py::overload_cast< ASTERINTEGER >( &Result::hasElementaryCharacteristics, py::const_ ),
            R"(
Test if a elementary characterictics is used at the specfied index

Arguments:
    index (int): index

Returns:
    bool: *True* if at least one elementary characterictics used else *False*.
        )",
            py::arg( "index" ) )
        .def( "hasElementaryCharacteristics",
              py::overload_cast<>( &Result::hasElementaryCharacteristics, py::const_ ) )
        .def( "hasListOfLoads",
              py::overload_cast< const ASTERINTEGER & >( &Result::hasListOfLoads, py::const_ ), R"(
Test if a list of loads is used at the specfied index

Arguments:
    index (int): index

Returns:
    bool: *True* if at least one list of loads is used else *False*.
        )",
              py::arg( "index" ) )
        .def( "hasListOfLoads", py::overload_cast<>( &Result::hasListOfLoads, py::const_ ) )
        .def( "hasModel", &Result::hasModel, R"(
Test if a model is used at the specfied index

Arguments:
    index (int): index

Returns:
    bool: *True* if at a model used else *False*.
        )",
              py::arg( "index" ) )
        .def( "hasMaterialField", &Result::hasMaterialField, R"(
Test if a material field is used at the specfied index

Arguments:
    index (int): index

Returns:
    bool: *True* if at a material field used else *False*.
        )",
              py::arg( "index" ) )
        .def( "getMaterialFields", &Result::getMaterialFields, R"(
Return the list of all material fields used in the result

Returns:
    list[MaterialField]: list of material field.
        )" )
        .def( "getMaterialField",
              py::overload_cast< ASTERINTEGER >( &Result::getMaterialField, py::const_ ),
              R"(
Return the material field for the given index.

Arguments:
    index (int): index

Returns:
    MaterialField: Material field.
              )",
              py::arg( "index" ) )
        .def( "getMaterialField", py::overload_cast<>( &Result::getMaterialField, py::const_ ) )
        .def( "getListOfLoads", &Result::getListOfLoads, R"(
Get list of loads on the specified index

Arguments:
    index (int): index to get

Returns:
    ListOfLoads: a pointer to list of loads.

        )",
              py::arg( "index" ) )
        .def( "getMesh", &Result::getMesh, R"(
Return a pointer to mesh

Returns:
    mesh (Mesh): a pointer to the mesh.
        )" )
        .def( "getModels", &Result::getModels, R"(
Return the list of all models used in the result

Returns:
    list[Model]: list of models.
        )" )
        .def( "getModel", py::overload_cast< ASTERINTEGER >( &Result::getModel, py::const_ ), R"(
Return the model for the given index.

Arguments:
    index (int): index

Returns:
    Model: Model object.
              )",
              py::arg( "index" ) )
        .def( "getModel", py::overload_cast<>( &Result::getModel, py::const_ ) )
        .def( "getNumberOfIndexes", &Result::getNumberOfIndexes, R"(
Get the number of index stored in the result

Returns:
    int: number of index stored.
        )" )
        .def( "getLastIndex", &Result::getLastIndex, R"(
Get the last index stored in the result

Returns:
    int: last index stored.
        )" )
        .def( "getFirstIndex", &Result::getFirstIndex, R"(
Get the first index stored in the result

Returns:
    int: first index stored.
        )" )
        .def( "getLastTime", &Result::getLastTime, R"(
Get the last time value stored in the result

Returns:
    float: last time value.
        )" )
        .def( "getParameters", &Result::getParameters, R"(
Return the parameters of the result as Python dict.

Returns:
    dict{str : list[int,float,str]}: Dict of values for each parameter variable.
        )",
              py::arg( "only_access" ) = false )
        .def( "getAccessParameters", &Result::getParameters, R"(
Return the access parameters of the result as Python dict.

Returns:
    dict{str : list[int,float,str]}: Dict of values for each access variable.
        )",
              py::arg( "only_access" ) = true )
        .def( "createIndexFromParameter",
              py::overload_cast< const std::string &, const std::string & >(
                  &Result::createIndexFromParameter ),
              R"(
Create an index in the result

Arguments:
    para_name (str): parameter name to store
    para_value (str): parameter value to store
        )",
              py::arg( "para_name" ), py::arg( "para_value" ) )
        .def( "getFieldsOnNodesRealNames", &Result::getFieldsOnNodesRealNames, R"(
Return the names of the real fields on nodes as Python list.

Returns:
    list(str): List of names of the real fields on nodes.
        )" )
        .def( "getFieldsOnCellsRealNames", &Result::getFieldsOnCellsRealNames, R"(
Return the names of the real fields on cells as Python list.

Returns:
    list(str): List of names of the real fields on cells.
        )" )
        .def( "getFieldsOnNodesComplexNames", &Result::getFieldsOnNodesComplexNames, R"(
Return the names of the complex fields on nodes as Python list.

Returns:
    list(str): List of names of the complex fields on nodes.
        )" )
        .def( "getFieldsOnCellsComplexNames", &Result::getFieldsOnCellsComplexNames, R"(
Return the names of the complex fields on cells as Python list.

Returns:
    list(str): List of names of the complex fields on cells.
        )" )
        .def( "getFieldsOnCellsLongNames", &Result::getFieldsOnCellsLongNames, R"(
Return the names of the integer fields on cells as Python list.

Returns:
    list(str): List of names of the integer fields on cells.
        )" )
        .def( "getConstantFieldsOnCellsChar16Names", &Result::getConstantFieldsOnCellsChar16Names,
              R"(
Return the names of the contant char16 fields on cells as Python list.

Returns:
    list(str): List of names of the contant fields on cells.
        )" )
        .def( "getConstantFieldsOnCellsRealNames", &Result::getConstantFieldsOnCellsRealNames,
              R"(
Return the names of the contant real fields on cells as Python list.

Returns:
    list(str): List of names of the contant fields on cells.
        )" )
        .def( "getGeneralizedVectorRealNames", &Result::getGeneralizedVectorRealNames, R"(
Return the names of the real generalized vectors as Python list.

Returns:
    list(str): List of names of the real generalized vectors.
        )" )
        .def( "getGeneralizedVectorComplexNames", &Result::getGeneralizedVectorComplexNames, R"(
Return the names of the complex generalized vectors as Python list.

Returns:
    list(str): List of names of the complex generalized vectors.
        )" )
        .def( "getIndexes", &Result::getIndexes, R"(
Return the list of indexes used to store fields

Returns:
    list[int]: List of indexs used to store fields.
        )" )
        .def( "getIndexesForFieldName", &Result::getIndexesForFieldName, R"(
            Returns the list of indexes used to store a specific field,
             indicated by its name.

            Returns:
                list[int]: List of indexs used to store fields.
                    )" )
        .def( "_interpolateFieldOnNodesReal", &Result::interpolateFieldOnNodesReal, R"(
Interpolate a FieldOnNodesReal from result.

* For internal use only - prefer to use interpolateField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    value (float) : acces variable value
    para (str, optional) : acces variable name (ex : 'INST'), (default : 'INST')
    left (str, optional) : extrapolation ('EXCLU', 'CONSTANT', 'LINEAIRE'), (default : 'EXCLU')
    right (str, optional) : extrapolation ('EXCLU', 'CONSTANT', 'LINEAIRE'), (default : 'EXCLU')
    crit (str, optional) : search precision criterion ('RELATIF', 'ABSOLU'), (default : 'RELATIF')
    prec (float, optional) : search precision value, (default : 1.e-6)
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    FieldOnNodesReal: interpolated field
        )",
              py::arg( "name" ), py::arg( "value" ), py::arg( "para" ) = "INST",
              py::arg( "left" ) = "EXCLU", py::arg( "right" ) = "EXCLU",
              py::arg( "crit" ) = "RELATIF", py::arg( "prec" ) = 1.e-6,
              py::arg( "updatePtr" ) = true )
        .def( "_interpolateFieldOnCellsReal", &Result::interpolateFieldOnCellsReal, R"(
Interpolate a FieldOnCellsReal from result.

* For internal use only - prefer to use interpolateField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    value (float) : acces variable value
    para (str, optional) : acces variable name (ex : 'INST'), (default : 'INST')
    left (str, optional) : extrapolation ('EXCLU', 'CONSTANT', 'LINEAIRE'), (default : 'EXCLU')
    right (str, optional) : extrapolation ('EXCLU', 'CONSTANT', 'LINEAIRE'), (default : 'EXCLU')
    crit (str, optional) : search precision criterion ('RELATIF', 'ABSOLU'), (default : 'RELATIF')
    prec (float, optional) : search precision value, (default : 1.e-6)
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    FieldOnCellsReal: interpolated field
        )",
              py::arg( "name" ), py::arg( "value" ), py::arg( "para" ) = "INST",
              py::arg( "left" ) = "EXCLU", py::arg( "right" ) = "EXCLU",
              py::arg( "crit" ) = "RELATIF", py::arg( "prec" ) = 1.e-6,
              py::arg( "updatePtr" ) = true )
        .def( "_getFieldOnNodesReal", &Result::getFieldOnNodesReal, R"(
Get a FieldOnNodesReal from result.
* For internal use only - prefer to use getField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    index (int): index to set the field
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    FieldOnNodesReal: field to get
        )",
              py::arg( "name" ), py::arg( "index" ), py::arg( "updatePtr" ) = true )
        .def( "_getFieldOnCellsReal", &Result::getFieldOnCellsReal, R"(
Get a FieldOnCellsReal from result.
* For internal use only - prefer to use getField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    index (int): index to set the field
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    FieldOnCellsReal: field to get
        )",
              py::arg( "name" ), py::arg( "index" ), py::arg( "updatePtr" ) = true )
        .def( "_getFieldOnNodesComplex", &Result::getFieldOnNodesComplex, R"(
Get a FieldOnNodesComplex from result.
* For internal use only - prefer to use getField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    index (int): index to set the field
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    FieldOnNodesComplex: field to get
        )",
              py::arg( "name" ), py::arg( "index" ), py::arg( "updatePtr" ) = true )
        .def( "_getFieldOnCellsComplex", &Result::getFieldOnCellsComplex, R"(
Get a FieldOnCellsComplex from result.
* For internal use only - prefer to use getField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    index (int): index to set the field
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    FieldOnCellsComplex: field to get
        )",
              py::arg( "name" ), py::arg( "index" ), py::arg( "updatePtr" ) = true )
        .def( "_getFieldOnCellsLong", &Result::getFieldOnCellsLong, R"(
Get a FieldOnCellsLong from result.
* For internal use only - prefer to use getField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    index (int): index to set the field
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    FieldOnCellsLong: field to get
        )",
              py::arg( "name" ), py::arg( "index" ), py::arg( "updatePtr" ) = true )
        .def( "_getConstantFieldOnCellsChar16", &Result::getConstantFieldOnCellsChar16, R"(
Get a ConstantFieldOnCellsChar16 from result.
* For internal use only - prefer to use getField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'COMPORTEMENT', ...)
    index (int): index to set the field
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    ConstantFieldOnCellsChar16: field to get
        )",
              py::arg( "name" ), py::arg( "index" ), py::arg( "updatePtr" ) = true )
        .def( "_getConstantFieldOnCellsReal", &Result::getConstantFieldOnCellsReal, R"(
Get a ConstantFieldOnCellsReal from result.
* For internal use only - prefer to use getField() *

Arguments:
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    index (int): index to set the field
    updatePtr (bool, optional): update pointer on values (default: True)

Returns:
    ConstantFieldOnCellsReal: field to get
        )",
              py::arg( "name" ), py::arg( "index" ), py::arg( "updatePtr" ) = true )
        .def( "printMedFile", &Result::printMedFile,
              R"(
Print the result in a MED file.

Args:
    filename (Path|str): Path to the output file.
    medname (str): Name of the result in the MED file. (default: "")
    local (bool): Print only the local domain if *True*. (default: True)
    fields (list[str]): Name of fields to save. (default: all)
    version (str): Version of MED file.
              )",
              py::arg( "filename" ), py::arg( "medname" ) = "", py::arg( "local" ) = true,
              py::arg( "internalVar" ) = true, py::arg( "fields" ) = VectorString(),
              py::arg( "version" ) = "" )
        .def( "setMesh", &Result::setMesh, R"(
Set the mesh used by the result.

Arguments:
    mesh (BaseMesh): mesh to set

        )",
              py::arg( "mesh" ) )
        .def( "_build", &Result::build, R"(
Build the result from the name of the result. It stores fields which are setted in c++ or
created in fortran

Arguments:
    feds (list[FiniteElementDescriptor]) : list of additional finite element descriptor used to
        build FieldOnCells
    fnds (list[EquationNumbering]) : list of additional field description used to
        build FieldOnNodes

Returns:
    bool: *True* if ok.
        )",
              py::arg( "feds" ) = std::vector< FiniteElementDescriptorPtr >(),
              py::arg( "fnds" ) = std::vector< EquationNumberingPtr >() )
        .def( "printInfo", &Result::printInfo )
        .def( "getFieldsNames", &Result::getFieldsNames, R"(
Return the list of names of stored fields

Returns:
    list[str]: List of names of stored fields.
        )" )
        .def( "resize", &Result::resize, R"(
Resize the object.

Arguments:
    nbIndexes (int): new expected size. Should be greater than the current size,
        otherwise the size is unchanged.
        )",
              py::arg( "nbIndexes" ) )
        .def(
            "_setField",
            py::overload_cast< const FieldOnNodesRealPtr, const std::string &, const ASTERINTEGER >(
                &Result::setField ),
            R"(
Set a real FieldOnNodes to result.

Arguments:
    field (FieldOnNodesReal): field to set
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    index (int): index to set the field
        )",
            py::arg( "field" ), py::arg( "name" ), py::arg( "index" ) )
        .def( "_setField",
              py::overload_cast< const FieldOnNodesComplexPtr, const std::string &,
                                 const ASTERINTEGER >( &Result::setField ),
              R"(
Set a complex FieldOnNodes to result.

Arguments:
    field (FieldOnNodesComplex): field to set
    name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
    index (int): index to set the field
        )",
              py::arg( "field" ), py::arg( "name" ), py::arg( "index" ) )
        .def(
            "_setField",
            py::overload_cast< const FieldOnCellsRealPtr, const std::string &, const ASTERINTEGER >(
                &Result::setField ),
            R"(
Set a real FieldOnCells to result

Arguments:
    field (FieldOnCellsReal): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    index (int): index to set the field
        )",
            py::arg( "field" ), py::arg( "name" ), py::arg( "index" ) )
        .def( "_setField",
              py::overload_cast< const FieldOnCellsComplexPtr, const std::string &,
                                 const ASTERINTEGER >( &Result::setField ),
              R"(
Set a complex FieldOnCells to result

Arguments:
    field (FieldOnCellsComplex): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    index (int): index to set the field
        )",
              py::arg( "field" ), py::arg( "name" ), py::arg( "index" ) )
        .def(
            "_setField",
            py::overload_cast< const FieldOnCellsLongPtr, const std::string &, const ASTERINTEGER >(
                &Result::setField ),
            R"(
Set a long FieldOnCells to result

Arguments:
    field (FieldOnCellsLong): field to set
    name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
    index (int): index to set the field
        )",
            py::arg( "field" ), py::arg( "name" ), py::arg( "index" ) )
        .def( "_setField",
              py::overload_cast< const ConstantFieldOnCellsChar16Ptr, const std::string &,
                                 const ASTERINTEGER >( &Result::setField ),
              R"(
Set a ConstantFieldOnCellsChar16 to result

Arguments:
    field (ConstantFieldOnCellsChar16): field to set
    name (str): symbolic name of the field in the result (ex: 'COMPOR', ...)
    index (int): index to set the field
        )",
              py::arg( "field" ), py::arg( "name" ), py::arg( "index" ) )
        .def( "_setField",
              py::overload_cast< const ConstantFieldOnCellsRealPtr, const std::string &,
                                 const ASTERINTEGER >( &Result::setField ),
              R"(
Set a ConstantFieldOnCellsReal to result

Arguments:
    field (ConstantFieldOnCellsReal): field to set
    name (str): symbolic name of the field in the result (ex: 'COMPOR', ...)
    index (int): index to set the field
        )",
              py::arg( "field" ), py::arg( "name" ), py::arg( "index" ) )
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) )
        .def( "setTable", &ListOfTables::setTable, R"(
Set a Table to the datastructure.

Arguments:
    identifier (str): Table identifier.
    Table: Table to be stored with the given identifier.
        )",
              py::arg( "identifier" ), py::arg( "table" ) )
        .def( "getFieldsNames", &Result::getFieldsNames, R"(
Return the list of names of stored fields

Returns:
    list[str]: List of names of stored fields.
        )" )
        .def( "getFiniteElementDescriptors", &Result::getFiniteElementDescriptors, R"(
Get list of finite element descriptor to build internal FieldOnCells

Returns:
    list[FiniteElementDescriptor]: list of finite element descriptor
        )" )
        .def( "getEquationNumberings", &Result::getEquationNumberings, R"(
Get list of field's description to build internal FieldOnNodes

Returns:
    list[EquationNumbering]: list of field's description
        )" )
        .def( "addFiniteElementDescriptor", &Result::addFiniteElementDescriptor );
};
