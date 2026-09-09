#ifndef CALCUL_H_
#define CALCUL_H_

/**
 * @file Calcul.h
 * @brief Header of Calcul class
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

#include "astercxx.h"

#include "Behaviours/BehaviourProperty.h"
#include "DataFields/ElementaryTerm.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "MemoryManager/JeveuxObject.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/HHOModel.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"

/**
 * @class Calcul
 * @brief This class allow computation of elementary terms (CALCUL.F90)
 */
class Calcul {
  private:
    typedef std::map< std::string, DataFieldPtr > listFields;
    typedef std::map< std::string, DataFieldPtr > listElemTerms;
    typedef std::map< std::string, bool > listExists;

    /** @brief Option to compute */
    std::string _option;

    /** @brief Flag to complete elementary terms on all CPU */
    bool _completeField;

    /** @brief Flag for stop if option doesn't exist for one element */
    bool _stopCompute;

    /** @brief Where to compute elementary terms */
    BaseMeshPtr _mesh;
    FiniteElementDescriptorPtr _FEDesc;

    /** @brief Input fields */
    listFields _inputFields;
    listElemTerms _inputElemTerms;

    /** @brief Output fields */
    listFields _outputFields;
    listExists _outputFieldsExist;

    /** @brief Output elementary terms */
    listElemTerms _outputElemTerms;
    listExists _outputElemTermsExist;

  private:
    /** @brief Create C++ objects for output fields */
    void postCompute();

  public:
    /** @typedef CalculPtr */
    typedef std::shared_ptr< Calcul > CalculPtr;

    /** @brief Constructor */
    Calcul( const std::string &option );

    /** @brief Destructor */
    ~Calcul() {};

    /** @brief Set option to compute */
    void setOption( const std::string &option ) { _option = option; };

    /** @brief Set flag for complete field on all CPUs */
    void setCompleteField( const bool flag ) { _completeField = flag; };

    /** @brief Set flag to stop if option doesn't exist for one element */
    void setStopCompute( const bool flag ) { _stopCompute = flag; };

    /** @brief Set computation on finite element descriptor */
    void setFiniteElementDescriptor( const FiniteElementDescriptorPtr FEDesc );

    /** @brief Set computation on model */
    void setModel( const ModelPtr &model );

    /** @brief Set computation on group of Cells */
    void setGroupsOfCells( const ModelPtr &model, const VectorString &groupOfCells );

    /** @brief Add input field */
    void addInputField( const std::string &parameterName, const DataFieldPtr field );

    /** @brief Add output field */
    void addOutputField( const std::string &parameterName, const DataFieldPtr field );

    /** @brief Add input elementary term */
    template < typename ValueType >
    void addInputElementaryTerm( const std::string &parameterName,
                                 const std::shared_ptr< ElementaryTerm< ValueType > > elemTerm ) {
        _inputElemTerms.insert( listElemTerms::value_type( parameterName, elemTerm ) );
    };

    /** @brief Add output elementary term */
    template < typename ValueType >
    void addOutputElementaryTerm( const std::string &parameterName,
                                  const std::shared_ptr< ElementaryTerm< ValueType > > elemTerm ) {
        _outputElemTerms.insert( listElemTerms::value_type( parameterName, elemTerm ) );
        _outputElemTermsExist.insert( listExists::value_type( parameterName, false ) );
    }

    /** @brief Clear all input fields */
    void clearInputFields() { _inputFields.clear(); };

    /** @brief Clear all input elementary term */
    void clearInputElemTerms() { _inputElemTerms.clear(); };

    /** @brief Clear all output fields */
    void clearOutputFields() {
        _outputFields.clear();
        _outputFieldsExist.clear();
    };

    /** @brief Clear all output elementary term */
    void clearOutputElemTerms() {
        _outputElemTerms.clear();
        _outputElemTermsExist.clear();
    };

    /** @brief Clear all outputs  */
    void clearOutputs() {
        clearOutputFields();
        clearOutputElemTerms();
    };

    /** @brief Clear all inputs  */
    void clearInputs() {
        clearInputFields();
        clearInputElemTerms();
    };

    /** @brief Is output if is elementary term */
    bool hasOutputElementaryTerm( const std::string &parameterName ) const;

    /** @brief Get output if is elementary term */
    ElementaryTermRealPtr getOutputElementaryTermReal( const std::string &parameterName ) const;
    ElementaryTermComplexPtr
    getOutputElementaryTermComplex( const std::string &parameterName ) const;

    /** @brief Compute option */
    void compute();

    /** @brief Add input fields for elementary characteristics */
    void addElementaryCharacteristicsField( const ElementaryCharacteristicsPtr elemChara );

    /** @brief Add input fields for orientation */
    void addOrientationField( const ElementaryCharacteristicsPtr elemChara );

    /** @brief Create and add input field for Fourier */
    void addFourierModeField( const ASTERINTEGER &nh );

    /** @brief Create and add input field for current time */
    void addTimeField( const std::string &parameterName, const ASTERDOUBLE time_value );

    /** @brief Create and add input field for current time (for thermics) */
    void addTimeField( const std::string &parameterName, const ASTERDOUBLE &time_value,
                       const ASTERDOUBLE &time_delta, const ASTERDOUBLE &time_theta );

    /** @brief Create and add input fields for XFEM */
    void addXFEMField( const XfemModelPtr xfemModel );

    /** @brief Create and add input fields for XFEM */
    void addXFEMField( const ModelPtr model );

    /** @brief Create and add input fields for HHO */
    void addHHOField( const HHOModelPtr hhoModel );

    /** @brief Create and add input fields for HHO */
    void addHHOField( const ModelPtr model );

    /** @brief Add input fields for non-linear behaviours */
    void addBehaviourField( const BehaviourPropertyPtr behaviour );

    /** @brief Get current FE descriptor */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() { return _FEDesc; };

    /** @brief Detect input fields with complex values */
    bool hasComplexInputFields( void ) const;
};

/**  @typedef CalculPtr */
using CalculPtr = std::shared_ptr< Calcul >;

#endif /* CALCUL_H_ */
