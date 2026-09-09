/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2026  EDF www.code-aster.org
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

#include "Discretization/DiscreteComputation.h"

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"

#include "Discretization/Calcul.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

ConstantFieldOnCellsRealPtr
DiscreteComputation::createTimeField( const ASTERDOUBLE time_value, const ASTERDOUBLE time_delta,
                                      const ASTERDOUBLE time_theta ) const {

    int sz = 1;
    if ( _phys_problem->getModel()->isThermal() )
        sz = 6;

    VectorString para_names = { "INST", "DELTAT", "THETA", "KHI", "R", "RHO" };
    VectorReal para_values = { time_value, time_delta, time_theta, 0., 0., 0. };

    VectorString reduced_names = VectorString( para_names.begin(), para_names.begin() + sz );
    VectorReal reduced_values = VectorReal( para_values.begin(), para_values.begin() + sz );

    // Get mesh
    auto mesh = _phys_problem->getMesh();

    // Create field
    auto field = std::make_shared< ConstantFieldOnCellsReal >( mesh );

    const std::string physicalName( "INST_R" );
    field->allocate( physicalName );
    ConstantFieldOnZone a( mesh );
    ConstantFieldValues< ASTERDOUBLE > b( reduced_names, reduced_values );
    field->setValueOnZone( a, b );

    return field;
}

CalculPtr DiscreteComputation::createCalculForNonLinear( const std::string option,
                                                         const ASTERDOUBLE &time_prev,
                                                         const ASTERDOUBLE &time_curr,
                                                         const FieldOnCellsRealPtr varc_prev,
                                                         const FieldOnCellsRealPtr varc_curr,
                                                         const VectorString &groupOfCells ) const {

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();
    auto currExternVarRefe = _phys_problem->getReferenceExternalStateVariables();
    AS_ASSERT( currMater );

    // No !
    if ( currModel->exists3DShell() ) {
        throw std::runtime_error( "COQUE_3D not implemented" );
    }
    if ( currModel->existsSTRX() ) {
        throw std::runtime_error( "Beams not implemented" );
    }

    // Prepare computing: the main object
    CalculPtr calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add external state variables
    if ( currMater->hasExternalStateVariable() ) {
        if ( !varc_prev || !varc_prev->exists() ) {
            raiseAsterError(
                "External state variables vector for beginning of time step is missing" );
        }
        if ( !varc_curr || !varc_curr->exists() ) {
            raiseAsterError( "External state variables vector for end of time step is missing" );
        }
        if ( currMater->hasExternalStateVariableWithReference() ) {
            AS_ASSERT( currExternVarRefe );
            calcul->addInputField( "PVARCRR", currExternVarRefe );
        }
        calcul->addInputField( "PVARCMR", varc_prev );
        calcul->addInputField( "PVARCPR", varc_curr );
    }

    // Add time fields
    calcul->addTimeField( "PINSTMR", time_prev );
    calcul->addTimeField( "PINSTPR", time_curr );

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
        calcul->addOrientationField( currElemChara );
    }
    calcul->addXFEMField( currModel );
    calcul->addBehaviourField( currBehaviour );

    return calcul;
};

ConstantFieldOnCellsLongPtr
DiscreteComputation::createDampingFluidField( const ASTERINTEGER damping,
                                              const ASTERINTEGER onde_flui ) const {

    VectorString para_names = { "X1", "X2" };
    VectorLong para_values = { damping, onde_flui };

    // Get mesh
    auto mesh = _phys_problem->getMesh();

    // Create field
    auto field = std::make_shared< ConstantFieldOnCellsLong >( mesh );

    const std::string physicalName( "NEUT_I" );
    field->allocate( physicalName );
    ConstantFieldOnZone a( mesh );
    ConstantFieldValues< ASTERINTEGER > b( para_names, para_values );
    field->setValueOnZone( a, b );

    return field;
};

ConstantFieldOnCellsLongPtr
DiscreteComputation::createWaveTypeFluidField( const ASTERINTEGER onde_flui ) const {

    VectorString para_names = { "X1" };
    VectorLong para_values = { onde_flui };

    // Get mesh
    auto mesh = _phys_problem->getMesh();

    // Create field
    auto field = std::make_shared< ConstantFieldOnCellsLong >( mesh );

    const std::string physicalName( "NEUT_I" );
    field->allocate( physicalName );
    ConstantFieldOnZone a( mesh );
    ConstantFieldValues< ASTERINTEGER > b( para_names, para_values );
    field->setValueOnZone( a, b );

    return field;
};
