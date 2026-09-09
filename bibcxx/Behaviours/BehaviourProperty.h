#ifndef BEHAVIOURPROPERTY_H_
#define BEHAVIOURPROPERTY_H_

/**
 * @file BehaviourProperty.h
 * @brief Header for class BehaviourProperty
 * @author Mickaël Abas
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

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Materials/MaterialField.h"
#include "Modeling/Model.h"

/**
 * @class BehaviourProperty
 * @brief Class to define behaviour
 */
class BehaviourProperty : public DataStructure {
  private:
    /** @brief Mesh */
    BaseMeshPtr _mesh;

    /** @brief Model */
    ModelPtr _model;

    /** @brief Material field */
    MaterialFieldPtr _materialField;

    /** @brief Flag for initial state */
    bool _initialState;

    /** @brief Verbosity */
    bool _verbosity;

    /** @brief Map '.COMPOR' to define behaviours */
    ConstantFieldOnCellsChar16Ptr _COMPOR;

    /** @brief Mao '.MULCOM' for multi-behaviours (crystals) */
    ConstantFieldOnCellsChar16Ptr _MULCOM;

    /** @brief Map '.CARCRI' for parameters to integrate behaviours */
    ConstantFieldOnCellsRealPtr _CARCRI;

    /** Flag for annealing */
    bool _annealing;

  private:
    /** @brief Create objects (maps) */
    void createObjects();

    /** @brief Detect functionnalities */
    void detectFunctionnalities();

  public:
    /** @typedef Smart-pointer to behaviour class */
    typedef std::shared_ptr< BehaviourProperty > BehaviourPropertyPtr;

    /** @brief Constructor */
    BehaviourProperty( const std::string name );

    /** @brief Constructor */
    BehaviourProperty();

    /** @brief Constructor */
    BehaviourProperty( const std::string name, ModelPtr model, MaterialFieldPtr materialField );

    /** @brief Constructor */
    BehaviourProperty( ModelPtr model, MaterialFieldPtr materialField );

    /** @brief Build object */
    bool build();

    /** @brief Get model */
    ModelPtr getModel() const { return _model; }

    /** @brief Get material field */
    MaterialFieldPtr getMaterialField() const { return _materialField; }

    /** @brief Set flag for initial state */
    void setInitialState( const bool &value ) { _initialState = value; };

    /** @brief Set verbosity */
    void setVerbosity( const bool &value ) { _verbosity = value; };

    /** @brief Set behavior ConstantField '.COMPOR' */
    ConstantFieldOnCellsChar16Ptr getBehaviourField() const { return _COMPOR; };

    /** @brief Set behavior ConstantField '.MULT_COMP' */
    ConstantFieldOnCellsChar16Ptr getMultipleBehaviourField() const { return _MULCOM; };

    /** @brief Set behavior ConstantField '.CARCRI' */
    ConstantFieldOnCellsRealPtr getConvergenceCriteria() const { return _CARCRI; };

    /** @brief Detect non-linear behaviour */
    bool hasBehaviour( const std::string &behaviour ) const;

    VectorString variNameToCmp( const VectorLong &cells, const VectorString &variNames ) const;

    /** @brief Detect annealing as post-treatment */
    bool hasAnnealing() const { return _annealing; };
};

/** @typedef Smart-pointer to behaviour class */
using BehaviourPropertyPtr = std::shared_ptr< BehaviourProperty >;

#endif
