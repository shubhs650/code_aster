/**
 * @file Calcul.cxx
 * @brief Implementation of Calcul
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

#include "Discretization/Calcul.h"

#include "aster_fort_utils.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"

/** @brief Constructor */
Calcul::Calcul( const std::string &option )
    : _option( option ),
      _stopCompute( true ),
      _completeField( true ),
      _FEDesc( nullptr ),
      _mesh( nullptr ) {};

/** @brief Compute on partial model */
void Calcul::setFiniteElementDescriptor( const FiniteElementDescriptorPtr FEDesc ) {
    _FEDesc = FEDesc;
    _mesh = _FEDesc->getMesh();
    if ( _mesh->isParallel() ) {
        _completeField = false;
    }
}

/** @brief Compute on all model */
void Calcul::setModel( const ModelPtr &model ) {
    _FEDesc = model->getFiniteElementDescriptor();
    _mesh = model->getMesh();
    if ( _mesh->isParallel() ) {
        _completeField = false;
    }
}

/** @brief Compute on a part of the model */
void Calcul::setGroupsOfCells( const ModelPtr &model, const VectorString &groupOfCells ) {
    _mesh = model->getMesh();
    _FEDesc = std::make_shared< FiniteElementDescriptor >( model->getFiniteElementDescriptor(),
                                                           groupOfCells );
    if ( _mesh->isParallel() ) {
        _completeField = false;
    }
}

/** @brief Add input field */
void Calcul::addInputField( const std::string &parameterName, const DataFieldPtr field ) {
    if ( field ) {
        _inputFields.insert( listFields::value_type( parameterName, field ) );
    }
}

/** @brief Add output field */
void Calcul::addOutputField( const std::string &parameterName, const DataFieldPtr field ) {
    AS_ASSERT( field );
    _outputFields.insert( listFields::value_type( parameterName, field ) );
    _outputFieldsExist.insert( listExists::value_type( parameterName, false ) );
}

/** @brief Add input fields for elementary characteristics */
void Calcul::addElementaryCharacteristicsField( const ElementaryCharacteristicsPtr elemChara ) {
    if ( elemChara ) {
        addInputField( "PCADISK", elemChara->getDiscreteRigidity() );
        addInputField( "PCADISM", elemChara->getDiscreteMass() );
        addInputField( "PCADISA", elemChara->getDiscreteDamping() );
        addInputField( "PCAGEPO", elemChara->getBeamGeometry() );
        addInputField( "PCAGNPO", elemChara->getBeamSection() );
        addInputField( "PCAARPO", elemChara->getFlexibilityCoefficients() );
        addInputField( "PCACABL", elemChara->getCableParameters() );
        addInputField( "PCAGNBA", elemChara->getBarParameters() );
        addInputField( "PCAPOUF", elemChara->getFluidBeamParameters() );
        addInputField( "PNBSP_I", elemChara->getNumberOfSubpoints() );
        addInputField( "PFIBRES", elemChara->getFibers() );
        addInputField( "PCINFDI", elemChara->getDiscreteParameters() );
    }
}

/** @brief Add input fields for orientation */
void Calcul::addOrientationField( const ElementaryCharacteristicsPtr elemChara ) {
    if ( elemChara ) {
        addInputField( "PCAORIE", elemChara->getLocalBasis() );
        addInputField( "PCACOQU", elemChara->getShellParameters() );
        addInputField( "PCAMASS", elemChara->getMaterialBase() );
    }
}

/** @brief Create and add input field for Fourier */
void Calcul::addFourierModeField( const ASTERINTEGER &nh ) {
    auto _FourierField = std::make_shared< ConstantFieldOnCellsLong >( _mesh );
    const std::string physicalName( "HARMON" );
    _FourierField->allocate( physicalName );
    ConstantFieldOnZone a( _mesh );
    ConstantFieldValues< ASTERINTEGER > b( { "NH" }, { nh } );
    _FourierField->setValueOnZone( a, b );
    addInputField( "PHARMON", _FourierField );
}

/** @brief Create and add input field for current time */
void Calcul::addTimeField( const std::string &parameterName, const ASTERDOUBLE time_value ) {
    auto _timeField = std::make_shared< ConstantFieldOnCellsReal >( _mesh );
    const std::string physicalName( "INST_R" );
    _timeField->allocate( physicalName );
    ConstantFieldOnZone a( _mesh );
    ConstantFieldValues< ASTERDOUBLE > b( { "INST" }, { time_value } );
    _timeField->setValueOnZone( a, b );
    addInputField( parameterName, _timeField );
}

/** @brief Create and add input field for current time */
void Calcul::addTimeField( const std::string &parameterName, const ASTERDOUBLE &time_value,
                           const ASTERDOUBLE &time_delta, const ASTERDOUBLE &time_theta ) {
    auto _timeField = std::make_shared< ConstantFieldOnCellsReal >( _mesh );
    const std::string physicalName( "INST_R" );
    _timeField->allocate( physicalName );
    ConstantFieldOnZone a( _mesh );
    ConstantFieldValues< ASTERDOUBLE > b( { "INST", "DELTAT", "THETA", "KHI", "R", "RHO" },
                                          { time_value, time_delta, time_theta, 0.0, 0.0, 0.0 } );
    _timeField->setValueOnZone( a, b );
    addInputField( parameterName, _timeField );
}

/** @brief Create and add input fields for XFEM */
void Calcul::addXFEMField( const XfemModelPtr xfemModel ) {
    if ( xfemModel ) {
        addInputField( "PAINTER", xfemModel->getField( "AINTER" ) );
        addInputField( "PPINTER", xfemModel->getField( "PINTER" ) );
        addInputField( "PPINTTO", xfemModel->getField( "PINTTO" ) );
        addInputField( "PCNSETO", xfemModel->getField( "CNSETO" ) );
        addInputField( "PHEAVTO", xfemModel->getField( "HEAVTO" ) );
        addInputField( "PLONCHA", xfemModel->getField( "LONCHA" ) );
        addInputField( "PBASLOR", xfemModel->getField( "BASLOC" ) );
        addInputField( "PLSN", xfemModel->getField( "LSN" ) );
        addInputField( "PLST", xfemModel->getField( "LST" ) );
        addInputField( "PSTANO", xfemModel->getField( "STANO" ) );
        addInputField( "PPMILTO", xfemModel->getField( "PMILT" ) );
        addInputField( "PFISNO", xfemModel->getField( "FISSNO" ) );
        addInputField( "PHEA_NO", xfemModel->getField( "HEAVNO" ) );
        addInputField( "PHEA_SE", xfemModel->getField( "HEAVSE" ) );
        addInputField( "PHEA_FA", xfemModel->getField( "HEAVFA" ) );
        addInputField( "PCFACE", xfemModel->getField( "CFACE" ) );
        addInputField( "PLONGCO", xfemModel->getField( "LONGCO" ) );
        addInputField( "PBASECO", xfemModel->getField( "BASECO" ) );
    }
}

/** @brief Create and add input fields for HHO */
void Calcul::addXFEMField( const ModelPtr model ) {
    if ( model->existsXfem() ) {
        addXFEMField( model->getXfemModel() );
    }
}

/** @brief Add input fields for non-linear behaviours */
void Calcul::addBehaviourField( const BehaviourPropertyPtr behaviour ) {
    if ( behaviour ) {
        addInputField( "PCOMPOR", behaviour->getBehaviourField() );
        addInputField( "PCARCRI", behaviour->getConvergenceCriteria() );
        addInputField( "PMULCOM", behaviour->getMultipleBehaviourField() );
    }
}

/** @brief Create and add input fields for HHO */
void Calcul::addHHOField( const HHOModelPtr HHOModel ) {
    if ( HHOModel ) {
        addInputField( "PCHHOGT", HHOModel->getGradient() );
        addInputField( "PCHHOST", HHOModel->getStabilization() );
        addInputField( "PCHHOBS", HHOModel->getBasis() );
    }
}

/** @brief Create and add input fields for HHO */
void Calcul::addHHOField( const ModelPtr model ) {
    if ( model && model->existsHHO() ) {
        addHHOField( model->getHHOModel() );
    }
}

/** @brief Compute option */
void Calcul::compute() {

#ifdef ASTER_DEBUG_CXX
    std::cout << "Computing option >> " << _option << std::endl;
#endif

    ASTERINTEGER inputNb = _inputFields.size() + _inputElemTerms.size();
    VectorString inputFields, inputParams;
    inputFields.reserve( inputNb );
    inputParams.reserve( inputNb );
    for ( const auto &[parameterName, field] : _inputFields ) {
#ifdef ASTER_DEBUG_CXX
        std::cout << "Input Field :" << parameterName << ", " << field->getName() << ", "
                  << field->exists() << std::endl;
#endif
        inputParams.push_back( parameterName );
        inputFields.push_back( field->getName() );
    }
    for ( const auto &[parameterName, elemTerm] : _inputElemTerms ) {
#ifdef ASTER_DEBUG_CXX
        std::cout << "Input Term :" << parameterName << ", " << elemTerm->getName() << std::endl;
#endif
        inputParams.push_back( parameterName );
        inputFields.push_back( elemTerm->getName() );
    }

    ASTERINTEGER outputNb = _outputFields.size() + _outputElemTerms.size();
    VectorString outputFields, outputParams;
    outputFields.reserve( outputNb );
    outputParams.reserve( outputNb );
    for ( const auto &[parameterName, field] : _outputFields ) {
        // #ifdef ASTER_DEBUG_CXX
        //         std::cout << "Output Field :" << parameterName << ", " << field->getName() <<
        //         std::endl;
        // #endif
        outputParams.push_back( parameterName );
        outputFields.push_back( field->getName() );
    }
    for ( const auto &[parameterName, elemTerm] : _outputElemTerms ) {
        // #ifdef ASTER_DEBUG_CXX
        //         std::cout << "Output Term :" << parameterName << ", " << elemTerm->getName() <<
        //         std::endl;
        // #endif
        outputParams.push_back( parameterName );
        outputFields.push_back( elemTerm->getName() );
    }

    std::string calculContinue = "C";
    if ( _stopCompute )
        calculContinue = "S";

    std::string calculMPI = "NON";
    if ( _completeField )
        calculMPI = "OUI";

    std::string calculLigrel = _FEDesc->getName();

    const std::string calculBase = JeveuxMemoryTypesNames[Permanent];

    char *inputFields_c = vectorStringAsFStrArray( inputFields, 19 );
    char *inputParams_c = vectorStringAsFStrArray( inputParams, 8 );
    char *outputFields_c = vectorStringAsFStrArray( outputFields, 19 );
    char *outputParams_c = vectorStringAsFStrArray( outputParams, 8 );

    // Compute
    CALLO_CALCUL( calculContinue.c_str(), _option.c_str(), calculLigrel.c_str(), &inputNb,
                  inputFields_c, inputParams_c, &outputNb, outputFields_c, outputParams_c,
                  calculBase.c_str(), calculMPI.c_str() );

    FreeStr( inputFields_c );
    FreeStr( inputParams_c );
    FreeStr( outputFields_c );
    FreeStr( outputParams_c );

    // Create C++ objects for output fields
    postCompute();
}

void Calcul::postCompute() {

    const std::string questi1( "TYPE_CHAMP" );
    const std::string typeco( "CHAMP" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );

    for ( auto &[parameterName, field] : _outputFields ) {
        std::string fieldName = field->getName();
        std::string fieldType = field->getFieldType();
        if ( fieldType == "ELEM" || fieldType == "ELNO" || fieldType == "ELGA" ) {
            std::string fieldScalar = field->getFieldScalar();
            if ( fieldScalar == "R" ) {
                std::static_pointer_cast< FieldOnCellsReal >( field )->setDescription( _FEDesc );
            } else if ( fieldScalar == "C" ) {
                std::static_pointer_cast< FieldOnCellsComplex >( field )->setDescription( _FEDesc );
            } else if ( fieldScalar == "I" ) {
                std::static_pointer_cast< FieldOnCellsLong >( field )->setDescription( _FEDesc );
            } else {
                AS_ABORT( "Field no supported" );
            }
            _outputFieldsExist.at( parameterName ) = true;
        } else {
            _outputFieldsExist.at( parameterName ) = false;
        }
    }

    for ( auto &[parameterName, elemTerm] : _outputElemTerms ) {
        std::string elemTermName = elemTerm->getName();
        CALLO_DISMOI( questi1, elemTermName, typeco, &repi, repk, arret, &ier );

        std::string fieldType( strip( repk.toString() ) );
        if ( fieldType == "RESL" ) {
            _outputElemTermsExist.at( parameterName ) = true;
            std::static_pointer_cast< ElementaryTermReal >( elemTerm )
                ->setFiniteElementDescriptor( _FEDesc );
        } else {
            _outputElemTermsExist.at( parameterName ) = false;
        }
    }
}

/** @brief Is output is elementary term */
bool Calcul::hasOutputElementaryTerm( const std::string &parameterName ) const {
    return _outputElemTermsExist.at( parameterName );
};

/** @brief Get output if is elementary term */
ElementaryTermRealPtr
Calcul::getOutputElementaryTermReal( const std::string &parameterName ) const {
    AS_ASSERT( hasOutputElementaryTerm( parameterName ) );
    return std::static_pointer_cast< ElementaryTermReal >( _outputElemTerms.at( parameterName ) );
};

/** @brief Get output if is elementary term */
ElementaryTermComplexPtr
Calcul::getOutputElementaryTermComplex( const std::string &parameterName ) const {
    AS_ASSERT( hasOutputElementaryTerm( parameterName ) );
    return std::static_pointer_cast< ElementaryTermComplex >(
        _outputElemTerms.at( parameterName ) );
};

/** @brief Detect input fields with complex values */
bool Calcul::hasComplexInputFields( void ) const {
    for ( const auto &[parameterName, field] : _inputFields ) {
        std::string fieldName = field->getName();
        std::string fieldType = field->getFieldType();
#ifdef ASTER_DEBUG_CXX
        std::cout << "Input Field :" << parameterName << ", " << fieldName << ", "
                  << field->exists() << std::endl;
#endif
        if ( fieldType == "ELEM" || fieldType == "ELNO" || fieldType == "ELGA" ) {
            std::string fieldScalar = field->getFieldScalar();
            if ( fieldScalar == "C" ) {
                return true;
            }
        }
    }
    return false;
};
