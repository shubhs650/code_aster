/**
 * @file Result.cxx
 * @brief Implementation de Result
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

#include "Results/Result.h"

#include "aster_fort_ds.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "Messages/Messages.h"
#include "ParallelUtilities/AsterMPI.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

ASTERINTEGER Result::_getInternalIndex( const ASTERINTEGER &storageIndex ) const {
    _serialNumber->updateValuePointer();
    const auto nbIndexes = getNumberOfIndexes();
    for ( ASTERINTEGER internalIndex = 0; internalIndex < nbIndexes; internalIndex++ ) {
        if ( ( *_serialNumber )[internalIndex] == storageIndex ) {
            return internalIndex;
        }
    }
    return -1;
};

std::string Result::_generateFieldName( const ASTERINTEGER &indexSymbName,
                                        const ASTERINTEGER &internalIndex ) const {

    // Generate name for field symbol
    auto nuch = to_string( indexSymbName, 3 );

    // Generate name for index
    auto chford = to_string( internalIndex, 6 );

    // Generate name of field
    auto fieldName = std::string( strip( getName() ) + "." + nuch + "." + chford );

    return fieldName;
};

template < typename T >
void Result::_setFieldBase(
    const std::string &symbName, const ASTERINTEGER &storageIndex, std::shared_ptr< T > field,
    std::map< std::string, std::map< ASTERINTEGER, std::shared_ptr< T > > > &dict ) {
    if ( !field )
        raiseAsterError( "ValueError: field is empty" );

    // Check mesh
    setMesh( field->getMesh() );

    // Check existence
    if ( !exists() ) {
        UTMESS( "F", "RESULT2_10" );
    }

    // Get index of this symbolic name
    auto indexSymbName = _symbolicNamesOfFields->getIndexFromString( strip( symbName ) );

    if ( indexSymbName == 0 ) {
        UTMESS( "F", "RESULT2_4", { strip( symbName ) } );
    }

    // Get internal index
    auto internalIndex = _getInternalIndex( storageIndex );

    // Add storage index
    if ( internalIndex < 0 ) {
        addStorageIndex( storageIndex );
        internalIndex = _getInternalIndex( storageIndex );
        AS_ASSERT( internalIndex >= 0 );
    }

    // Generate internal name of field
    std::string internalName;
    internalName = _generateFieldName( indexSymbName, internalIndex );

    if ( !_namesOfFields->isBuilt() ) {
        _namesOfFields->build( true );
    }

    // Note field in datastructure
    JeveuxCollectionObjectChar24 storageStructure = ( *_namesOfFields )[indexSymbName];
    storageStructure->updateValuePointer();
    storageStructure[internalIndex] = ljust( internalName, 24, ' ' );

    // Save field in dictionnary
    if ( dict.count( strip( symbName ) ) == 0 ) {
        dict[strip( symbName )] = std::map< ASTERINTEGER, std::shared_ptr< T > >();
    }

    // if field already exist, destroy it before to create new one
    if ( dict[strip( symbName )].count( storageIndex ) > 0 ) {
        dict[strip( symbName )][storageIndex] = nullptr;
    }

    // Create smart pointer
    if ( internalName == field->getName() ) {
        dict[strip( symbName )][storageIndex] = field;
    } else {
        auto result = std::make_shared< T >( internalName, *field );
        dict[strip( symbName )][storageIndex] = result;
    }
};

void Result::_checkMesh( const BaseMeshPtr &mesh ) const {
    if ( !mesh )
        raiseAsterError( "ValueError: Mesh is empty" );

    if ( _mesh ) {
        if ( _mesh != mesh )
            raiseAsterError( "Incompatible meshes" );
    }
}

void Result::setMesh( const BaseMeshPtr &mesh ) {
    _checkMesh( mesh );
    _mesh = mesh;
};

void Result::setElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara,
                                           ASTERINTEGER storageIndex, bool exists_ok ) {
    if ( !cara )
        raiseAsterError( "ValueError: ElementaryCharacteristics is empty" );
    if ( _mapElemCara.count( storageIndex ) != 0 ) {
        if ( exists_ok || _mapElemCara[storageIndex] == cara )
            return;
        raiseAsterError( "ValueError: ElementaryCharacteristics already assigned at index " +
                         std::to_string( storageIndex ) );
    }
    // Check existence
    if ( !exists() ) {
        UTMESS( "F", "RESULT2_10" );
    }
    _mapElemCara[storageIndex] = cara;
    std::string type( "CARAELEM" );
    std::string cel( "E" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &storageIndex, cara->getName(), type, cel );
    setMesh( cara->getMesh() );
};

void Result::setListOfLoads( const ListOfLoadsPtr &load, ASTERINTEGER storageIndex ) {
    if ( !load )
        raiseAsterError( "ValueError: Load is empty" );
    if ( _mapLoads.count( storageIndex ) != 0 ) {
        if ( _mapLoads[storageIndex] == load )
            return;
        raiseAsterError( "ValueError: Load already assigned at index " +
                         std::to_string( storageIndex ) );
    }
    // Check existence
    if ( !exists() ) {
        UTMESS( "F", "RESULT2_10" );
    }
    _mapLoads[storageIndex] = load;
    std::string type( "EXCIT" );
    std::string cel( "E" );
    CALLO_RSADPA_ZK24_WRAP( getName(), &storageIndex, load->getName(), type, cel );
};

void Result::setMaterialField( const MaterialFieldPtr &mater, ASTERINTEGER storageIndex,
                               bool exists_ok ) {
    if ( !mater )
        raiseAsterError( "ValueError: MaterialField is empty" );
    if ( _mapMaterial.count( storageIndex ) != 0 ) {
        if ( exists_ok || _mapMaterial[storageIndex] == mater )
            return;
        raiseAsterError( "ValueError: MaterialField already assigned at index " +
                         std::to_string( storageIndex ) );
    }
    // Check existence
    if ( !exists() ) {
        UTMESS( "F", "RESULT2_10" );
    }
    _mapMaterial[storageIndex] = mater;
    std::string type( "CHAMPMAT" );
    std::string cel( "E" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &storageIndex, mater->getName(), type, cel );
    setMesh( mater->getMesh() );
};

void Result::setModel( const ModelPtr &model, ASTERINTEGER storageIndex, bool exists_ok ) {
    if ( !model )
        raiseAsterError( "ValueError: Model is empty" );
    if ( _mapModel.count( storageIndex ) != 0 ) {
        if ( exists_ok || _mapModel[storageIndex] == model )
            return;
        raiseAsterError( "ValueError: Model already assigned at index " +
                         std::to_string( storageIndex ) );
    }
    // Check existence
    if ( !exists() ) {
        UTMESS( "F", "RESULT2_10" );
    }

    _mapModel[storageIndex] = model;
    std::string type( "MODELE" );
    std::string cel( "E" );
    CALLO_RSADPA_ZK8_WRAP( getName(), &storageIndex, model->getName(), type, cel );
    const auto fed = model->getFiniteElementDescriptor();
    _fieldBuilder.addFiniteElementDescriptor( fed );
    setMesh( model->getMesh() );
};

void Result::setParameterValue( std::string paraName, ASTERDOUBLE paraValue,
                                ASTERINTEGER storageIndex ) {
    CALLO_RSADPA_ZR_WRAP( getName(), &storageIndex, &paraValue, paraName );
};

void Result::setParameterValue( std::string paraName, std::string paraValue,
                                ASTERINTEGER storageIndex ) {

    auto paraIndx = _dictParameters.find( paraName );
    if ( paraIndx == _dictParameters.end() ) {
        raiseAsterError( "Parameter not available" );
    }

    auto paraType = _dictParameters.find( paraName )->second;
    std::string cel( "E" );
    if ( paraType == "Char8" ) {
        CALLO_RSADPA_ZK8_WRAP( getName(), &storageIndex, paraValue, paraName, cel );
    } else if ( paraType == "Char16" ) {
        CALLO_RSADPA_ZK16_WRAP( getName(), &storageIndex, paraValue, paraName, cel );
    } else if ( paraType == "Char24" ) {
        CALLO_RSADPA_ZK24_WRAP( getName(), &storageIndex, paraValue, paraName, cel );
    } else {
        raiseAsterError( "Wrapper not available" );
    }
};

ASTERDOUBLE Result::getTime( ASTERINTEGER storageIndex ) const {

    _rspr->updateValuePointer();
    if ( !_calculationParameter->isBuilt() ) {
        _calculationParameter->build( true );
    }

    int i = 1;
    for ( const auto &item : *_calculationParameter ) {
        item->updateValuePointer();
        auto typevar = strip( ( *item )[3].toString() );

        if ( typevar == "ACCES" ) {
            auto var_name = strip( _accessVariables->getStringFromIndex( i ) );
            if ( var_name == "INST" ) {
                auto nosuff = strip( ( *item )[0].toString() );
                auto ivar = std::stoi( strip( ( *item )[1].toString() ) );
                auto nmax = std::stoi( strip( ( *item )[2].toString() ) );

                AS_ASSERT( nosuff == ".RSPR" )

                auto internalIndex = _getInternalIndex( storageIndex );
                if ( internalIndex < 0 ) {
                    AS_ABORT( "Error: internal index not found" );
                }

                return ( *_rspr )[nmax * internalIndex + ivar - 1];
            }
        }
        ++i;
    }
    UTMESS( "F", "RESULT2_9" );
    return 0.0;
};

void Result::allocate( ASTERINTEGER nbIndexes ) {

    std::string base( JeveuxMemoryTypesNames[Permanent] );
    ASTERINTEGER nbordr = nbIndexes;
    CALLO_RSCRSD( base, getName(), getType(), &nbordr );

    AS_ASSERT( _calculationParameter->build( true ) );
    AS_ASSERT( _namesOfFields->build( true ) );
    _listOfParameters();
};

void Result::_listOfParameters() {

    _calculationParameter->updateValuePointer();

    int i = 1;
    for ( const auto &item : *_calculationParameter ) {
        item->updateValuePointer();
        auto objectSuffix = strip( ( *item )[0].toString() );
        auto paraName = strip( _accessVariables->getStringFromIndex( i ) );
        std::string paraType;
        if ( objectSuffix == ".RSPR" ) {
            paraType = "Real";
        } else if ( objectSuffix == ".RSPI" ) {
            paraType = "Integer";
        } else if ( objectSuffix == ".RSP8" ) {
            paraType = "Char8";
        } else if ( objectSuffix == ".RS16" ) {
            paraType = "Char16";
        } else if ( objectSuffix == ".RS24" ) {
            paraType = "Char24";
        } else {
            AS_ABORT( "Unknown type" );
        }
        _dictParameters.insert( std::make_pair( paraName, paraType ) );
        ++i;
    }
}

void Result::setElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara,
                                           bool exists_ok ) {
    auto allStorageIndexes = getIndexes();
    for ( auto &storageIndex : allStorageIndexes ) {
        setElementaryCharacteristics( cara, storageIndex, exists_ok );
    }
};

void Result::setMaterialField( const MaterialFieldPtr &mater, bool exists_ok ) {
    auto allStorageIndexes = getIndexes();
    for ( auto &storageIndex : allStorageIndexes ) {
        setMaterialField( mater, storageIndex, exists_ok );
    }
};

void Result::setModel( const ModelPtr &model, bool exists_ok ) {
    const auto fed = model->getFiniteElementDescriptor();
    _fieldBuilder.addFiniteElementDescriptor( fed );
    auto allStorageIndexes = getIndexes();
    for ( auto &storageIndex : allStorageIndexes ) {
        setModel( model, storageIndex, exists_ok );
    }
};

std::vector< ElementaryCharacteristicsPtr > Result::getAllElementaryCharacteristics() const {
    return unique( _mapElemCara );
};

ElementaryCharacteristicsPtr Result::getElementaryCharacteristics() const {
    const auto cara = getAllElementaryCharacteristics();
    if ( cara.size() > 1 ) {
        UTMESS( "F", "RESULT2_8" );
    }

    if ( cara.size() == 1 )
        return cara[0];

    return ElementaryCharacteristicsPtr( nullptr );
};

ElementaryCharacteristicsPtr
Result::getElementaryCharacteristics( ASTERINTEGER storageIndex ) const {
    return _mapElemCara.at( storageIndex );
};

bool Result::hasElementaryCharacteristics( ASTERINTEGER storageIndex ) const {
    return _mapElemCara.count( storageIndex ) > 0;
};

bool Result::hasElementaryCharacteristics() const { return !_mapElemCara.empty(); };

bool Result::hasListOfLoads() const { return !_mapLoads.empty(); };

bool Result::hasListOfLoads( const ASTERINTEGER &storageIndex ) const {
    return _mapLoads.count( storageIndex ) > 0;
};

ListOfLoadsPtr Result::getListOfLoads( ASTERINTEGER storageIndex ) const {
    return _mapLoads.at( storageIndex );
};

std::vector< MaterialFieldPtr > Result::getMaterialFields() const {
    return unique( _mapMaterial );
};

MaterialFieldPtr Result::getMaterialField() const {
    const auto mate = getMaterialFields();
    if ( mate.size() > 1 ) {
        UTMESS( "F", "RESULT2_7" );
    };

    if ( mate.size() == 1 )
        return mate[0];

    return MaterialFieldPtr( nullptr );
};

MaterialFieldPtr Result::getMaterialField( ASTERINTEGER storageIndex ) const {
    return _mapMaterial.at( storageIndex );
};

BaseMeshPtr Result::getMesh() const {
    if ( _mesh != nullptr )
        return _mesh;
    const auto model = getModel();
    if ( model != nullptr )
        return model->getMesh();
    return BaseMeshPtr( nullptr );
};

bool Result::hasMultipleModel() const {
    std::string name( "" );
    for ( const auto &curIter : _mapModel ) {
        if ( name == "" ) {
            name = curIter.second->getName();
        }
        if ( name != curIter.second->getName() )
            return true;
    }
    return false;
}

bool Result::hasModel( const ASTERINTEGER &storageIndex ) const {
    return _mapModel.count( storageIndex ) > 0;
}

bool Result::hasMaterialField( const ASTERINTEGER &storageIndex ) const {
    return _mapMaterial.count( storageIndex ) > 0;
}

std::vector< ModelPtr > Result::getModels() const { return unique( _mapModel ); };

ModelPtr Result::getModel() const {
    if ( hasMultipleModel() ) {
        UTMESS( "F", "RESULT2_6" );
    }

    const auto models = getModels();
    AS_ASSERT( models.size() <= 1 );
    auto indexes = getIndexes();
    if ( models.size() == 1 ) {
        ModelPtr toReturn = _mapModel.at( indexes[0] );
        return toReturn;
    }

    return ModelPtr( nullptr );
};

ModelPtr Result::getModel( ASTERINTEGER storageIndex ) const {
    return _mapModel.at( storageIndex );
};

ASTERINTEGER Result::getNumberOfIndexes() const { return _serialNumber->size(); };

VectorLong Result::getIndexes() const { return _serialNumber->toVector(); };

ASTERINTEGER Result::getLastIndex() const {
    _serialNumber->updateValuePointer();
    return ( *_serialNumber )[_serialNumber->size() - 1];
};

ASTERDOUBLE Result::getLastTime() const { return getTime( getLastIndex() ); };

ASTERINTEGER Result::getFirstIndex() const {
    _serialNumber->updateValuePointer();
    return ( *_serialNumber )[0];
};

FieldOnCellsRealPtr Result::getFieldOnCellsReal( const std::string name,
                                                 const ASTERINTEGER storageIndex,
                                                 const bool updatePtr ) const {
    FieldOnCellsRealPtr result = _dictOfMapOfFieldOnCellsReal.at( name ).at( storageIndex );
    if ( updatePtr )
        result->updateValuePointers();
    return result;
};

FieldOnCellsComplexPtr Result::getFieldOnCellsComplex( const std::string name,
                                                       const ASTERINTEGER storageIndex,
                                                       const bool updatePtr ) const {
    FieldOnCellsComplexPtr result = _dictOfMapOfFieldOnCellsComplex.at( name ).at( storageIndex );
    if ( updatePtr )
        result->updateValuePointers();
    return result;
};

FieldOnCellsLongPtr Result::getFieldOnCellsLong( const std::string name,
                                                 const ASTERINTEGER storageIndex,
                                                 const bool updatePtr ) const {
    FieldOnCellsLongPtr result = _dictOfMapOfFieldOnCellsLong.at( name ).at( storageIndex );
    if ( updatePtr )
        result->updateValuePointers();
    return result;
};

ConstantFieldOnCellsChar16Ptr
Result::getConstantFieldOnCellsChar16( const std::string name, const ASTERINTEGER storageIndex,
                                       const bool updatePtr ) const {
    ConstantFieldOnCellsChar16Ptr result =
        _dictOfMapOfConstantFieldOnCellsChar16.at( name ).at( storageIndex );
    if ( updatePtr )
        result->updateValuePointers();
    return result;
};

ConstantFieldOnCellsRealPtr Result::getConstantFieldOnCellsReal( const std::string name,
                                                                 const ASTERINTEGER storageIndex,
                                                                 const bool updatePtr ) const {
    ConstantFieldOnCellsRealPtr result =
        _dictOfMapOfConstantFieldOnCellsReal.at( name ).at( storageIndex );
    if ( updatePtr )
        result->updateValuePointers();
    return result;
};

py::dict Result::getParameters( bool only_access ) const {

    py::dict returnDict;
    std::string var_name, str_val, typevar, nosuff;
    ASTERINTEGER ivar, nmax, index;

    CALL_JEMARQ();
    _serialNumber->updateValuePointer();
    if ( _rspr.exists() )
        _rspr->updateValuePointer();
    if ( _rspi.exists() )
        _rspi->updateValuePointer();
    if ( _rsp8.exists() )
        _rsp8->updateValuePointer();
    if ( _rs16.exists() )
        _rs16->updateValuePointer();
    if ( _rs24.exists() )
        _rs24->updateValuePointer();

    AS_ASSERT( _calculationParameter->build() );

    ASTERINTEGER nbIndexes = getNumberOfIndexes();

    var_name = "NUME_ORDRE";
    py::list listValues;
    for ( ASTERINTEGER internalIndex = 0; internalIndex < nbIndexes; ++internalIndex ) {
        listValues.append( ( *_serialNumber )[internalIndex] );
    }
    returnDict[var_name.c_str()] = listValues;

    auto items = _calculationParameter->getObjects();
    for ( auto &item : items ) {
        item->updateValuePointer();
        typevar = strip( ( *item )[3].toString() );

        if ( typevar == "ACCES" or ( typevar == "PARA" and !only_access ) ) {
            var_name = strip( _accessVariables->getStringFromIndex( item->getIndex() ) );
            nosuff = strip( ( *item )[0].toString() );
            ivar = std::stoi( strip( ( *item )[1].toString() ) );
            nmax = std::stoi( strip( ( *item )[2].toString() ) );

            py::list listV;

            if ( nosuff == ".RSPI" ) {
                for ( ASTERINTEGER j = 0; j < nbIndexes; ++j ) {
                    index = nmax * ( j ) + ivar - 1;
                    listV.append( ( *_rspi )[index] );
                }
            }

            else if ( nosuff == ".RSPR" ) {
                for ( ASTERINTEGER j = 0; j < nbIndexes; ++j ) {
                    index = nmax * ( j ) + ivar - 1;
                    listV.append( ( *_rspr )[index] );
                }
            }

            else {
                for ( ASTERINTEGER j = 0; j < nbIndexes; ++j ) {
                    index = nmax * ( j ) + ivar - 1;
                    if ( nosuff == ".RSP8" ) {
                        str_val = strip( ( ( *_rsp8 )[index] ).toString() );
                    } else if ( nosuff == ".RS16" ) {
                        str_val = strip( ( ( *_rs16 )[index] ).toString() );
                    } else if ( nosuff == ".RS24" ) {
                        str_val = strip( ( ( *_rs24 )[index] ).toString() );
                    } else {
                        AS_ASSERT( false );
                    }

                    if ( str_val.length() == 0 ) {
                        listV.append( py::none() );
                    } else {
                        listV.append( str_val );
                    }
                }
            }
            returnDict[var_name.c_str()] = listV;
        }
    }
    CALL_JEDEMA();
    return returnDict;
}

VectorString Result::getFieldsOnNodesRealNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnNodesReal.size() );

    for ( auto &it : _dictOfMapOfFieldOnNodesReal ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnNodesComplexNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnNodesComplex.size() );

    for ( auto &it : _dictOfMapOfFieldOnNodesComplex ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnCellsRealNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnCellsReal.size() );

    for ( auto &it : _dictOfMapOfFieldOnCellsReal ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnCellsComplexNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnCellsComplex.size() );

    for ( auto &it : _dictOfMapOfFieldOnCellsComplex ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

VectorString Result::getFieldsOnCellsLongNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfFieldOnCellsLong.size() );

    for ( auto &it : _dictOfMapOfFieldOnCellsLong ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

VectorString Result::getConstantFieldsOnCellsChar16Names() const {
    VectorString names;
    names.reserve( _dictOfMapOfConstantFieldOnCellsChar16.size() );

    for ( auto &it : _dictOfMapOfConstantFieldOnCellsChar16 ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

VectorString Result::getConstantFieldsOnCellsRealNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfConstantFieldOnCellsReal.size() );

    for ( auto &it : _dictOfMapOfConstantFieldOnCellsReal ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

VectorString Result::getGeneralizedVectorRealNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfGeneralizedVectorReal.size() );

    for ( auto &it : _dictOfMapOfGeneralizedVectorReal ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

VectorString Result::getGeneralizedVectorComplexNames() const {
    VectorString names;
    names.reserve( _dictOfMapOfGeneralizedVectorComplex.size() );

    for ( auto &it : _dictOfMapOfGeneralizedVectorComplex ) {
        std::string name = it.first;
        names.push_back( strip( name ) );
    }
    return names;
};

FieldOnNodesRealPtr Result::getFieldOnNodesReal( const std::string name,
                                                 const ASTERINTEGER storageIndex,
                                                 const bool updatePtr ) const {
    FieldOnNodesRealPtr result = _dictOfMapOfFieldOnNodesReal.at( name ).at( storageIndex );
    if ( updatePtr )
        result->updateValuePointers();
    return result;
};

FieldOnNodesComplexPtr Result::getFieldOnNodesComplex( const std::string name,
                                                       const ASTERINTEGER storageIndex,
                                                       const bool updatePtr ) const {

    FieldOnNodesComplexPtr result = _dictOfMapOfFieldOnNodesComplex.at( name ).at( storageIndex );
    if ( updatePtr )
        result->updateValuePointers();
    return result;
};

void Result::setField( const FieldOnNodesRealPtr field, const std::string &name,
                       const ASTERINTEGER storageIndex ) {
    _setFieldBase( name, storageIndex, field, _dictOfMapOfFieldOnNodesReal );
    _fieldBuilder.addEquationNumbering( field->getDescription() );
};

void Result::setField( const FieldOnNodesComplexPtr field, const std::string &name,
                       const ASTERINTEGER storageIndex ) {
    _setFieldBase( name, storageIndex, field, _dictOfMapOfFieldOnNodesComplex );
    _fieldBuilder.addEquationNumbering( field->getDescription() );
};

void Result::setField( const FieldOnCellsRealPtr field, const std::string &name,
                       const ASTERINTEGER storageIndex ) {
    _setFieldBase( name, storageIndex, field, _dictOfMapOfFieldOnCellsReal );
    _fieldBuilder.addFiniteElementDescriptor( field->getDescription() );
};

void Result::setField( const FieldOnCellsComplexPtr field, const std::string &name,
                       const ASTERINTEGER storageIndex ) {
    _setFieldBase( name, storageIndex, field, _dictOfMapOfFieldOnCellsComplex );
    _fieldBuilder.addFiniteElementDescriptor( field->getDescription() );
};

void Result::setField( const FieldOnCellsLongPtr field, const std::string &name,
                       const ASTERINTEGER storageIndex ) {
    _setFieldBase( name, storageIndex, field, _dictOfMapOfFieldOnCellsLong );
    _fieldBuilder.addFiniteElementDescriptor( field->getDescription() );
};

void Result::setField( const ConstantFieldOnCellsChar16Ptr field, const std::string &name,
                       const ASTERINTEGER storageIndex ) {
    _setFieldBase( name, storageIndex, field, _dictOfMapOfConstantFieldOnCellsChar16 );
};

void Result::setField( const ConstantFieldOnCellsRealPtr field, const std::string &name,
                       const ASTERINTEGER storageIndex ) {
    _setFieldBase( name, storageIndex, field, _dictOfMapOfConstantFieldOnCellsReal );
};

VectorString Result::getFieldsNames() const {
    VectorString vField;
    for ( auto curIter : _dictOfMapOfFieldOnNodesReal ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfFieldOnNodesComplex ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfFieldOnCellsReal ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfFieldOnCellsComplex ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfFieldOnCellsLong ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfConstantFieldOnCellsChar16 ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfConstantFieldOnCellsReal ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfGeneralizedVectorReal ) {
        vField.push_back( curIter.first );
    }
    for ( auto curIter : _dictOfMapOfGeneralizedVectorComplex ) {
        vField.push_back( curIter.first );
    }
    return vField;
}

VectorLong Result::getIndexesForFieldName( const std::string &name ) const {
    VectorLong keys;

    // Lambda to extract keys
    auto extractKeys = [&keys]( const auto &dict, const std::string &field_name ) -> bool {
        auto it = dict.find( field_name );
        if ( it != dict.end() ) {
            keys.reserve( it->second.size() ); // Optimize memory allocation
            std::transform( it->second.begin(), it->second.end(), std::back_inserter( keys ),
                            []( const auto &pair ) { return pair.first; } );
            return true;
        }
        return false;
    };

    if ( extractKeys( _dictOfMapOfFieldOnNodesComplex, name ) ||
         extractKeys( _dictOfMapOfFieldOnNodesReal, name ) ||
         extractKeys( _dictOfMapOfFieldOnCellsReal, name ) ||
         extractKeys( _dictOfMapOfFieldOnCellsComplex, name ) ||
         extractKeys( _dictOfMapOfFieldOnCellsLong, name ) ||
         extractKeys( _dictOfMapOfConstantFieldOnCellsChar16, name ) ||
         extractKeys( _dictOfMapOfConstantFieldOnCellsReal, name ) ||
         extractKeys( _dictOfMapOfGeneralizedVectorReal, name ) ||
         extractKeys( _dictOfMapOfGeneralizedVectorComplex, name ) ) {
        return keys;
    }

    return keys;
}

void Result::printListOfFields() const {
    auto vField = getFieldsNames();
    std::cout << "Content of DataStructure : ";
    for ( auto field : vField ) {
        std::cout << field << " - ";
    }
    std::cout << std::endl;
}

void Result::printInfo() const {
    ASTERINTEGER umess( 6 );
    CALLO_RSINFO( getName(), &umess );
}

FieldOnNodesRealPtr
Result::interpolateFieldOnNodesReal( const std::string name, const ASTERDOUBLE value,
                                     const std::string para, const std::string left,
                                     const std::string right, const std::string crit,
                                     const ASTERDOUBLE prec, const bool updatePtr ) const {

    VectorString namesok = getFieldsOnNodesRealNames();
    VectorString paraok = { "INST" };
    VectorString extrok = { "EXCLU", "LINEAIRE", "CONSTANT" };
    VectorString critok = { "ABSOLU", "RELATIF" };

    if ( std::count( namesok.begin(), namesok.end(), name ) == 0 ) {
        AS_ABORT( "Invalid input " + name );
    }
    if ( std::count( paraok.begin(), paraok.end(), para ) == 0 ) {
        AS_ABORT( "Invalid input " + para );
    }
    if ( std::count( extrok.begin(), extrok.end(), left ) == 0 ) {
        AS_ABORT( "Invalid input " + left );
    }
    if ( std::count( extrok.begin(), extrok.end(), right ) == 0 ) {
        AS_ABORT( "Invalid input " + right );
    }
    if ( std::count( critok.begin(), critok.end(), crit ) == 0 ) {
        AS_ABORT( "Invalid input " + crit );
    }
    if ( prec < 0.0 ) {
        AS_ABORT( "Invalid input " + std::to_string( prec ) );
    }

    const auto model = this->getModel();
    FieldOnNodesRealPtr result = std::make_shared< FieldOnNodesReal >( model );

    std::string base( "G" );
    ASTERINTEGER iret = 100;
    ASTERINTEGER stop = 0;

    CALLO_RSINCH( this->getName(), name, para, &value, result->getName(), right, left, &stop, base,
                  &prec, crit, &iret );
    if ( iret > 2 ) {
        AS_ABORT( "Interpolation error " + std::to_string( iret ) );
    }
    if ( updatePtr )
        result->updateValuePointers();
    return result;
}

FieldOnCellsRealPtr
Result::interpolateFieldOnCellsReal( const std::string name, const ASTERDOUBLE value,
                                     const std::string para, const std::string left,
                                     const std::string right, const std::string crit,
                                     const ASTERDOUBLE prec, const bool updatePtr ) const {

    VectorString namesok = getFieldsOnCellsRealNames();
    VectorString paraok = { "INST" };
    VectorString extrok = { "EXCLU", "LINEAIRE", "CONSTANT" };
    VectorString critok = { "ABSOLU", "RELATIF" };

    if ( std::count( namesok.begin(), namesok.end(), name ) == 0 ) {
        AS_ABORT( "Invalid input " + name );
    }
    if ( std::count( paraok.begin(), paraok.end(), para ) == 0 ) {
        AS_ABORT( "Invalid input " + para );
    }
    if ( std::count( extrok.begin(), extrok.end(), left ) == 0 ) {
        AS_ABORT( "Invalid input " + left );
    }
    if ( std::count( extrok.begin(), extrok.end(), right ) == 0 ) {
        AS_ABORT( "Invalid input " + right );
    }
    if ( std::count( critok.begin(), critok.end(), crit ) == 0 ) {
        AS_ABORT( "Invalid input " + crit );
    }
    if ( prec < 0.0 ) {
        AS_ABORT( "Invalid input " + std::to_string( prec ) );
    }

    const auto fed = this->getModel()->getFiniteElementDescriptor();
    FieldOnCellsRealPtr result = std::make_shared< FieldOnCellsReal >( fed );

    std::string base( "G" );
    ASTERINTEGER iret = 100;
    ASTERINTEGER stop = 0;

    CALLO_RSINCH( this->getName(), name, para, &value, result->getName(), right, left, &stop, base,
                  &prec, crit, &iret );
    if ( iret > 2 ) {
        AS_ABORT( "Interpolation error " + std::to_string( iret ) );
    }
    if ( updatePtr )
        result->updateValuePointers();
    return result;
};

void Result::printMedFile( const std::filesystem::path &fileName, std::string medName, bool local,
                           bool internalVar, const VectorString &fields,
                           std::string version ) const {
    const auto rank = getMPIRank();
    LogicalUnitFile a;
    ASTERINTEGER retour = -1;
    // In case that the print file (single and absolute path) is unique between processors,
    // it must only be created on proc 0.
    if ( getMesh()->isParallel() || rank == 0 ) {
        if ( rank == 0 )
            a.openFile( fileName, Binary, New );
        if ( getMesh()->isParallel() ) {
#ifdef ASTER_HAVE_MPI
            AsterMPI::barrier();
#endif /* ASTER_HAVE_MPI */
            if ( rank != 0 ) {
                auto mode = local ? New : Old;
                a.openFile( fileName, Binary, mode );
            }
        }
        retour = a.getLogicalUnit();
    }

    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;
    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = retour;
    if ( !version.empty() ) {
        dict.container["VERSION_MED"] = version;
    }

    if ( getMesh()->isParallel() ) {
        dict.container["PROC0"] = "NON";
        if ( !local )
            dict.container["FICHIER_UNIQUE"] = "OUI";
    } else
        dict.container["PROC0"] = "OUI";

    ListSyntaxMapContainer listeResu;
    SyntaxMapContainer dict2;
    dict2.container["RESULTAT"] = getName();
    dict2.container["TOUT_ORDRE"] = "OUI";
    if ( !medName.empty() )
        dict2.container["NOM_RESU_MED"] = medName.substr( 0, 8 );
    if ( !internalVar ) {
        dict2.container["IMPR_NOM_VARI"] = "NON";
    }
    if ( fields.empty() ) {
        dict2.container["TOUT_CHAM"] = "OUI";
    } else {
        dict2.container["NOM_CHAM"] = fields;
    }
    listeResu.push_back( dict2 );
    dict.container["RESU"] = listeResu;

    cmdSt.define( dict );

    ASTERINTEGER op = 39;
    CALL_EXECOP( &op );
};

bool Result::addFiniteElementDescriptor( const FiniteElementDescriptorPtr curFED ) {
    _fieldBuilder.addFiniteElementDescriptor( curFED );
    return true;
}

bool Result::build( const std::vector< FiniteElementDescriptorPtr > feds,
                    const std::vector< EquationNumberingPtr > fnds ) {
    CALL_JEMARQ();
    _serialNumber->updateValuePointer();

    AS_ASSERT( _calculationParameter->build( true ) );
    AS_ASSERT( _namesOfFields->build( true ) );

    auto nbIndexes = getNumberOfIndexes();

    for ( auto &fed : feds ) {
        _fieldBuilder.addFiniteElementDescriptor( fed );
    }

    for ( auto &fnd : fnds ) {
        _fieldBuilder.addEquationNumbering( fnd );
    }

    ASTERINTEGER cmpt = 1;
    for ( const auto &obj : _namesOfFields ) {
        obj->updateValuePointer();
        auto nomSymb = strip( _symbolicNamesOfFields->getStringFromIndex( cmpt ) );
        AS_ASSERT( nbIndexes <= obj->size() );
        for ( ASTERINTEGER internalIndex = 0; internalIndex < nbIndexes; ++internalIndex ) {
            std::string name( strip( ( *obj )[internalIndex].toString() ) );
            if ( name != "" ) {
                const ASTERINTEGER storageIndex = ( *_serialNumber )[internalIndex];
                addFieldFromString( nomSymb, name, storageIndex );
            }
        }
        ++cmpt;
    }

    // add of listofloads
    std::string type = "EXCIT";
    if ( _accessVariables->getIndexFromString( type ) > 0 ) {
        auto indexes = getIndexes();
        std::string cel( "L" );
        for ( auto &storageIndex : indexes ) {
            std::string paraValue( 24, ' ' );
            CALLO_RSADPA_ZK24_WRAP( getName(), &storageIndex, paraValue, type, cel );
            std::string name = paraValue.substr( 0, 19 );
            // only if created by a fortran command
            if ( name.substr( 0, 8 ) != getName().substr( 0, 8 ) )
                continue;
            mapIndexLoads::iterator it;
            for ( it = _mapLoads.begin(); it != _mapLoads.end(); it++ ) {
                if ( name == it->second->getName() )
                    break;
            }
            if ( it == _mapLoads.end() ) {
                _mapLoads[storageIndex] =
                    std::make_shared< ListOfLoads >( name, this->getModel( storageIndex ) );
            } else
                _mapLoads[storageIndex] = it->second;
        }
    }

    CALL_JEDEMA();
    return update_tables();
};

void Result::addFieldFromString( const std::string &nomSymb, const std::string &name,
                                 ASTERINTEGER storageIndex ) {
    CALL_JEMARQ();
    std::string questi( "TYPE_CHAMP" );
    const std::string typeco( "CHAMP" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );

    CALLO_DISMOI( questi, name, typeco, &repi, repk, arret, &ier );
    const std::string resu( strip( repk.toString() ) );

    questi = "TYPE_SCA";
    repk = " ";
    CALLO_DISMOI( questi, name, typeco, &repi, repk, arret, &ier );
    const std::string scalaire( strip( repk.toString() ) );

    if ( resu == "NOEU" ) {
        AS_ASSERT( _mesh != nullptr );
        if ( scalaire == "R" ) {
            if ( _dictOfMapOfFieldOnNodesReal.count( nomSymb ) == 0 ) {
                _dictOfMapOfFieldOnNodesReal[nomSymb] = MapOfFieldOnNodesReal();
            }

            if ( _dictOfMapOfFieldOnNodesReal[nomSymb].count( storageIndex ) == 0 ) {
                FieldOnNodesRealPtr result =
                    _fieldBuilder.buildFieldOnNodes< ASTERDOUBLE >( name, _mesh );
                _dictOfMapOfFieldOnNodesReal[nomSymb][storageIndex] = result;
            }
        } else if ( scalaire == "C" ) {
            if ( _dictOfMapOfFieldOnNodesComplex.count( nomSymb ) == 0 ) {
                _dictOfMapOfFieldOnNodesComplex[nomSymb] = MapOfFieldOnNodesComplex();
            }

            if ( _dictOfMapOfFieldOnNodesComplex[nomSymb].count( storageIndex ) == 0 ) {
                FieldOnNodesComplexPtr result =
                    _fieldBuilder.buildFieldOnNodes< ASTERCOMPLEX >( name, _mesh );
                _dictOfMapOfFieldOnNodesComplex[nomSymb][storageIndex] = result;
            }
        } else {
            AS_ABORT( "Type not supported: " + scalaire );
        }

    } else if ( resu == "ELEM" || resu == "ELNO" || resu == "ELGA" ) {
        if ( scalaire == "R" ) {
            if ( _dictOfMapOfFieldOnCellsReal.count( nomSymb ) == 0 ) {
                _dictOfMapOfFieldOnCellsReal[nomSymb] = MapOfFieldOnCellsReal();
            }

            if ( _dictOfMapOfFieldOnCellsReal[nomSymb].count( storageIndex ) == 0 ) {
                AS_ASSERT( _mesh != nullptr );
                auto result = _fieldBuilder.buildFieldOnCells< ASTERDOUBLE >( name, _mesh );
                _dictOfMapOfFieldOnCellsReal[nomSymb][storageIndex] = result;
            }
        } else if ( scalaire == "C" ) {
            if ( _dictOfMapOfFieldOnCellsComplex.count( nomSymb ) == 0 ) {
                _dictOfMapOfFieldOnCellsComplex[nomSymb] = MapOfFieldOnCellsComplex();
            }

            if ( _dictOfMapOfFieldOnCellsComplex[nomSymb].count( storageIndex ) == 0 ) {
                AS_ASSERT( _mesh != nullptr );
                auto result = _fieldBuilder.buildFieldOnCells< ASTERCOMPLEX >( name, _mesh );
                _dictOfMapOfFieldOnCellsComplex[nomSymb][storageIndex] = result;
            }
        } else if ( scalaire == "I" ) {
            if ( _dictOfMapOfFieldOnCellsLong.count( nomSymb ) == 0 ) {
                _dictOfMapOfFieldOnCellsLong[nomSymb] = MapOfFieldOnCellsLong();
            }

            if ( _dictOfMapOfFieldOnCellsLong[nomSymb].count( storageIndex ) == 0 ) {
                AS_ASSERT( _mesh != nullptr );
                auto result = _fieldBuilder.buildFieldOnCells< ASTERINTEGER >( name, _mesh );
                _dictOfMapOfFieldOnCellsLong[nomSymb][storageIndex] = result;
            }
        } else {
            AS_ABORT( "Type not supported: " + scalaire );
        }
    } else if ( resu == "CART" ) {
        if ( scalaire == "K16" ) {
            if ( _dictOfMapOfConstantFieldOnCellsChar16.count( nomSymb ) == 0 ) {
                _dictOfMapOfConstantFieldOnCellsChar16[nomSymb] = MapOfConstantFieldOnCellsChar16();
            }

            if ( _dictOfMapOfConstantFieldOnCellsChar16[nomSymb].count( storageIndex ) == 0 ) {
                AS_ASSERT( _mesh != nullptr );
                auto result =
                    _fieldBuilder.buildConstantFieldOnCells< JeveuxChar16 >( name, _mesh );
                _dictOfMapOfConstantFieldOnCellsChar16[nomSymb][storageIndex] = result;
            }
        } else if ( scalaire == "R" ) {
            if ( _dictOfMapOfConstantFieldOnCellsReal.count( nomSymb ) == 0 ) {
                _dictOfMapOfConstantFieldOnCellsReal[nomSymb] = MapOfConstantFieldOnCellsReal();
            }

            if ( _dictOfMapOfConstantFieldOnCellsReal[nomSymb].count( storageIndex ) == 0 ) {
                AS_ASSERT( _mesh != nullptr );
                auto result = _fieldBuilder.buildConstantFieldOnCells< ASTERDOUBLE >( name, _mesh );
                _dictOfMapOfConstantFieldOnCellsReal[nomSymb][storageIndex] = result;
            }
        } else {
            AS_ABORT( "Type not supported: " + scalaire );
        }
    } else if ( resu == "VGEN" ) {
        if ( scalaire == "R" ) {
            if ( _dictOfMapOfGeneralizedVectorReal.count( nomSymb ) == 0 ) {
                _dictOfMapOfGeneralizedVectorReal[nomSymb] = MapOfGeneralizedVectorReal();
            }

            if ( _dictOfMapOfGeneralizedVectorReal[nomSymb].count( storageIndex ) == 0 ) {
                GeneralizedAssemblyVectorRealPtr result(
                    new GeneralizedAssemblyVectorReal( name ) );
                _dictOfMapOfGeneralizedVectorReal[nomSymb][storageIndex] = result;
            }
        } else if ( scalaire == "C" ) {
            if ( _dictOfMapOfGeneralizedVectorComplex.count( nomSymb ) == 0 ) {
                _dictOfMapOfGeneralizedVectorComplex[nomSymb] = MapOfGeneralizedVectorComplex();
            }

            if ( _dictOfMapOfGeneralizedVectorComplex[nomSymb].count( storageIndex ) == 0 ) {
                GeneralizedAssemblyVectorComplexPtr result(
                    new GeneralizedAssemblyVectorComplex( name ) );
                _dictOfMapOfGeneralizedVectorComplex[nomSymb][storageIndex] = result;
            }
        } else {
            AS_ABORT( "Type not supported: " + scalaire );
        }
    } else {
        std::cout << "Field not build : " << name << " (" << resu << ")" << std::endl;
    }
    CALL_JEDEMA();
};

bool Result::exists() const { return _symbolicNamesOfFields.exists(); };

void Result::resize( ASTERINTEGER nbIndexes ) {
    if ( !exists() ) {
        allocate( nbIndexes );
    } else {
        CALLO_RSAGSD( getName(), &nbIndexes );
        AS_ASSERT( _calculationParameter->build( true ) );
        AS_ASSERT( _namesOfFields->build( true ) );
    }
}

void Result::clear( const ASTERINTEGER &storageIndex ) {

    auto old_index = getIndexes();

    ASTERINTEGER nume_ordre = storageIndex;
    CALLO_RSRUSD( getName(), &nume_ordre );

    for ( auto &index_2 : old_index ) {
        if ( index_2 >= storageIndex ) {
            _mapModel.erase( index_2 );
            _mapMaterial.erase( index_2 );
            _mapLoads.erase( index_2 );
            _mapElemCara.erase( index_2 );

            for ( auto &[key, fields] : _dictOfMapOfFieldOnNodesReal ) {
                fields.erase( index_2 );
            }
            for ( auto &[key, fields] : _dictOfMapOfFieldOnNodesComplex ) {
                fields.erase( index_2 );
            }
            for ( auto &[key, fields] : _dictOfMapOfFieldOnCellsReal ) {
                fields.erase( index_2 );
            }
            for ( auto &[key, fields] : _dictOfMapOfFieldOnCellsComplex ) {
                fields.erase( index_2 );
            }
            for ( auto &[key, fields] : _dictOfMapOfFieldOnCellsLong ) {
                fields.erase( index_2 );
            }
            for ( auto &[key, fields] : _dictOfMapOfConstantFieldOnCellsReal ) {
                fields.erase( index_2 );
            }
            for ( auto &[key, fields] : _dictOfMapOfConstantFieldOnCellsChar16 ) {
                fields.erase( index_2 );
            }
            for ( auto &[key, fields] : _dictOfMapOfGeneralizedVectorReal ) {
                fields.erase( index_2 );
            }
            for ( auto &[key, fields] : _dictOfMapOfGeneralizedVectorComplex ) {
                fields.erase( index_2 );
            }
        }
    }
};

void Result::clear() {
    auto storageIndex = getFirstIndex();
    CALLO_RSRUSD( getName(), &storageIndex );

    _mapModel.clear();
    _mapMaterial.clear();
    _mapLoads.clear();
    _mapElemCara.clear();

    _dictOfMapOfFieldOnNodesReal.clear();
    _dictOfMapOfFieldOnNodesComplex.clear();
    _dictOfMapOfFieldOnCellsReal.clear();
    _dictOfMapOfFieldOnCellsComplex.clear();
    _dictOfMapOfFieldOnCellsLong.clear();
    _dictOfMapOfConstantFieldOnCellsReal.clear();
    _dictOfMapOfConstantFieldOnCellsChar16.clear();
};

std::vector< FiniteElementDescriptorPtr > Result::getFiniteElementDescriptors() const {
    return _fieldBuilder.getFiniteElementDescriptors();
};

std::vector< EquationNumberingPtr > Result::getEquationNumberings() const {
    return _fieldBuilder.getEquationNumberings();
};

void Result::setTime( const ASTERDOUBLE &time, ASTERINTEGER storageIndex ) {
    setParameterValue( "INST", time, storageIndex );
}

void Result::addStorageIndex( const ASTERINTEGER storageIndex ) {

    CALL_JEMARQ();
    auto nbIndexes = getNumberOfIndexes();
    auto nbIndexMaxi = _serialNumber->capacity();
    if ( nbIndexes >= nbIndexMaxi ) {
        Result::resize( std::max( (ASTERINTEGER)1, 2 * nbIndexes ) );
    }
    _serialNumber->updateValuePointer();
    nbIndexes = getNumberOfIndexes();
    if ( nbIndexes >= 1 ) {
        auto lastIndex = getLastIndex();
        if ( storageIndex <= lastIndex ) {
            UTMESS( "F", "RESULT2_5" );
        }
    }
    _serialNumber->push_back( storageIndex );
    CALL_JEDEMA();
};

ASTERINTEGER Result::createIndexFromParameter( const std::string &paraName,
                                               const std::string &paraValue ) {
    CALL_JEMARQ();

    ASTERINTEGER internalIndex;
    internalIndex = getParameterIndex( paraName, paraValue );
    if ( internalIndex != -1 ) {
        AS_ABORT( "Parameter exists" );
    }

    // Test and resize object if require
    auto nbIndexes = getNumberOfIndexes();
    auto nbIndexMaxi = _serialNumber->capacity();
    if ( nbIndexes >= nbIndexMaxi ) {
        Result::resize( std::max( (ASTERINTEGER)1, 2 * getNumberOfIndexes() ) );
    }

    _serialNumber->updateValuePointer();
    nbIndexes = getNumberOfIndexes();

    // Generate new storage index
    ASTERINTEGER storageIndex;
    if ( nbIndexes >= 1 ) {
        auto lastIndex = getLastIndex();
        storageIndex = lastIndex + 1;
    } else {
        storageIndex = 1;
    }

    // Save new storage index
    addStorageIndex( storageIndex );

    // Add parameter
    setParameterValue( paraName, paraValue, storageIndex );
    CALL_JEDEMA();
    return storageIndex;
};

ASTERINTEGER Result::getParameterIndex( std::string paraName, std::string paraValue ) {

    ASTERINTEGER internalIndex;

    internalIndex = -1;

    auto paraIndx = _dictParameters.find( paraName );
    if ( paraIndx == _dictParameters.end() ) {
        raiseAsterError( "Parameter not available" );
    }

    auto paraType = _dictParameters.find( paraName )->second;
    std::string cel( "L" );
    std::string valueInResult( 24, ' ' );

    auto allStorageIndexes = getIndexes();
    for ( auto storageIndex : allStorageIndexes ) {

        if ( paraType == "Char8" ) {

            CALLO_RSADPA_ZK8_WRAP( getName(), &storageIndex, valueInResult, paraName, cel );
        } else if ( paraType == "Char16" ) {
            CALLO_RSADPA_ZK16_WRAP( getName(), &storageIndex, valueInResult, paraName, cel );
        } else if ( paraType == "Char24" ) {
            CALLO_RSADPA_ZK24_WRAP( getName(), &storageIndex, valueInResult, paraName, cel );
        } else {
            raiseAsterError( "Wrapper not available" );
        }
        if ( strip( paraValue ) == strip( valueInResult ) ) {

            return _getInternalIndex( storageIndex );
        }
    }
    return internalIndex;
}
