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

#include "aster_fort_calcul.h"
#include "aster_fort_superv.h"

#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/ThermalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/HHO.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

ElementaryMatrixTemperatureRealPtr DiscreteComputation::getLinearConductivityMatrix(
    const ASTERDOUBLE time_curr, const ASTERINTEGER &modeFourier,
    const FieldOnCellsRealPtr varc_curr, const VectorString &groupOfCells,
    const bool &with_dual ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "RIGI_THER" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem->getModel(), option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Prepare computing
    CalculPtr calcul = std::make_shared< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );

    calcul->addHHOField( currModel );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            if ( !varc_curr || !varc_curr->exists() ) {
                raiseAsterError( "External state variables are needed but not given" );
            }
            calcul->addInputField( "PVARCPR", varc_curr );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
        calcul->addOrientationField( currElemChara );
    }

    calcul->addFourierModeField( modeFourier );
    calcul->addTimeField( "PINSTR", time_curr, 0.0, 0.0 );

    calcul->addXFEMField( currModel );

    calcul->addHHOField( currModel );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for rigidity
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
    };

    if ( with_dual ) {
        DiscreteComputation::baseDualLinearConductivityMatrix( calcul, elemMatr );
    }

    elemMatr->build();
    return elemMatr;
};

/**
 * @brief Compute elementary matrices for mass matrix (RIGI_THER_TANG)
 */
ElementaryMatrixTemperatureRealPtr DiscreteComputation::getTangentConductivityMatrix(
    const FieldOnNodesRealPtr temp_prev, const FieldOnNodesRealPtr temp_step,
    const FieldOnCellsRealPtr varc_curr, const VectorString &groupOfCells,
    const bool &with_dual ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "RIGI_THER_TANG" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem->getModel(), option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehav = _phys_problem->getBehaviourProperty();

    // Prepare computing
    auto calcul = std::make_shared< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            if ( !varc_curr || !varc_curr->exists() ) {
                raiseAsterError( "External state variables are needed but not given" );
            }
            calcul->addInputField( "PVARCPR", varc_curr );
        }
    }

    if ( currBehav ) {
        calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
        calcul->addOrientationField( currElemChara );
    }

    calcul->addXFEMField( currModel );

    calcul->addHHOField( currModel );

    // Current Thermal Field
    auto temp_curr = std::make_shared< FieldOnNodesReal >( *temp_prev + *temp_step );
    calcul->addInputField( "PTEMPEI", temp_curr );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PMATTSR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
        if ( calcul->hasOutputElementaryTerm( "PMATTSR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTSR" ) );
    };

    if ( with_dual ) {
        DiscreteComputation::baseDualLinearConductivityMatrix( calcul, elemMatr );
    }

    elemMatr->build();
    return elemMatr;
}

ElementaryMatrixTemperatureRealPtr
DiscreteComputation::getLinearCapacityMatrix( const ASTERDOUBLE time_curr,
                                              const FieldOnCellsRealPtr varc_curr,
                                              const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "MASS_THER" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem->getModel(), option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addHHOField( currModel );
    // Set to -1 because not used.
    calcul->addTimeField( "PINSTR", time_curr, 0.0, 0.0 );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            if ( !varc_curr || !varc_curr->exists() ) {
                raiseAsterError( "External state variables are needed but not given" );
            }
            calcul->addInputField( "PVARCPR", varc_curr );
        }
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
        calcul->addOrientationField( currElemChara );
    }

    calcul->addXFEMField( currModel );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
    };

    elemMatr->build();
    return elemMatr;
};

/**
 * @brief Compute elementary matrices for mass matrix (MASS_THER_TANG)
 */
ElementaryMatrixTemperatureRealPtr DiscreteComputation::getTangentCapacityMatrix(
    const FieldOnNodesRealPtr temp_prev, const FieldOnNodesRealPtr temp_step,
    const FieldOnCellsRealPtr varc_curr, const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "MASS_THER_TANG" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem->getModel(), option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehav = _phys_problem->getBehaviourProperty();

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( option );
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addHHOField( currModel );

    if ( currMater ) {
        calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

        if ( currMater->hasExternalStateVariable() ) {
            if ( !varc_curr || !varc_curr->exists() ) {
                raiseAsterError( "External state variables are needed but not given" );
            }
            calcul->addInputField( "PVARCPR", varc_curr );
        }
    }

    if ( currBehav ) {
        calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
    }

    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
        calcul->addOrientationField( currElemChara );
    }

    calcul->addXFEMField( currModel );

    // Current Thermal Field
    auto temp_curr = std::make_shared< FieldOnNodesReal >( *temp_prev + *temp_step );
    calcul->addInputField( "PTEMPEI", temp_curr );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) )
            elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
    };

    elemMatr->build();
    return elemMatr;
}

void DiscreteComputation::baseDualLinearConductivityMatrix(
    CalculPtr &calcul, ElementaryMatrixTemperatureRealPtr &elemMatr ) const {

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();

    // Select option
    calcul->setOption( "THER_DDLM_R" );

    auto impl = [calcul, elemMatr]( auto loads ) {
        for ( const auto &load : loads ) {
            auto FEDesc = load->getFiniteElementDescriptor();
            auto field = load->getMultiplicativeField();
            if ( field && field->exists() && FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( FEDesc );
                calcul->addInputField( "PDDLMUR", field );
                calcul->addOutputElementaryTerm( "PMATTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                    elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
                }
            }
        }
    };

    impl( listOfLoads->getThermalLoadsReal() );
    impl( listOfLoads->getThermalLoadsFunction() );

#ifdef ASTER_HAVE_MPI
    impl( listOfLoads->getParallelThermalLoadsReal() );
    impl( listOfLoads->getParallelThermalLoadsFunction() );
#endif
};

ElementaryMatrixTemperatureRealPtr DiscreteComputation::getDualLinearConductivityMatrix() const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    const std::string option( "THER_DDLM_R" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem->getModel(), option );

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );

    // Compute elementary matrices
    DiscreteComputation::baseDualLinearConductivityMatrix( calcul, elemMatr );

    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixTemperatureRealPtr
DiscreteComputation::getThermalExchangeMatrix( const ASTERDOUBLE &time_curr ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    const std::string option( "RIGI_THER" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem->getModel(), option );

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );

    auto currModel = _phys_problem->getModel();

    // Compute elementary matrices
    // Prepare loads
    const auto &_listOfLoads = _phys_problem->getListOfLoads();
    auto isXfem = _phys_problem->getModel()->existsXfem();

    auto therLoadReal = _listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        auto FEDesc = _phys_problem->getModel()->getFiniteElementDescriptor();
        auto load_FEDesc = load->getFiniteElementDescriptor();

        if ( load->hasLoadResult() ) {
            auto exchange_field = load->interpolateLoadResult( "COEF_H", time_curr );
            calcul->setOption( "RIGI_THER_ECHA_R" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PCOEFHR", exchange_field );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, 1.0 );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }

        if ( load->hasLoadField( "COEFH" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            calcul->setOption( "RIGI_THER_ECHA_R" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PCOEFHR", exchange_field );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, 1.0 );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }

        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "RIGI_THER_PARO_R" );
            calcul->clearInputs();
            calcul->clearOutputs();
            if ( isXfem ) {
                XfemModelPtr currXfemModel = _phys_problem->getModel()->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PHECHPR", wall_exchange_field );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, 1.0 );

            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }
    }

    auto therLoadFunc = _listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        auto FEDesc = _phys_problem->getModel()->getFiniteElementDescriptor();
        auto load_FEDesc = load->getFiniteElementDescriptor();

        if ( load->hasLoadField( "COEFH" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            calcul->setOption( "RIGI_THER_ECHA_F" );
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( FEDesc );
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PCOEFHF", exchange_field );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, 1.0 );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }

        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "RIGI_THER_PARO_F" );
            calcul->clearInputs();
            calcul->clearOutputs();
            if ( isXfem ) {
                XfemModelPtr currXfemModel = _phys_problem->getModel()->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addInputField( "PGEOMER",
                                   _phys_problem->getModel()->getMesh()->getCoordinates() );
            calcul->addInputField( "PHECHPF", wall_exchange_field );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, 1.0 );
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }
    }
    elemMatr->build();
    return elemMatr;
};

ElementaryMatrixTemperatureRealPtr DiscreteComputation::getThermalTangentNonLinearNeumannMatrix(
    const FieldOnNodesRealPtr temp_curr, const ASTERDOUBLE time_curr,
    const FieldOnCellsRealPtr varc_curr ) const {

    const std::string calcul_option( "MTAN_THER" );

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemMatr = std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem->getModel(),
                                                                         calcul_option );

    // Setup

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl = [&]( auto load, const std::string &option, const std::string &name,
                     const std::string &param, const FiniteElementDescriptorPtr FED ) {
        if ( load->hasLoadField( name ) ) {
            calcul->setOption( option );
            calcul->setFiniteElementDescriptor( FED );

            calcul->clearInputs();
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPEI", temp_curr );
            calcul->addHHOField( currModel );

            if ( currMater && currMater->hasExternalStateVariable() ) {
                if ( !varc_curr || !varc_curr->exists() ) {
                    raiseAsterError( "External state variables are needed but not given" );
                }
                calcul->addInputField( "PVARCPR", varc_curr );
            }

            calcul->addInputField( param, load->getConstantLoadField( name ) );

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                auto elem_term = calcul->getOutputElementaryTermReal( "PMATTTR" );
                elemMatr->addElementaryTerm( elem_term );
            }
        }
    };

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        impl( load, "MTAN_THER_RAYO_R", "RAYO", "PRAYONR", model_FEDesc );
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        impl( load, "MTAN_THER_FLUXNL", "FLUNL", "PFLUXNL", model_FEDesc );
        impl( load, "MTAN_THER_RAYO_F", "RAYO", "PRAYONF", model_FEDesc );
    }

    elemMatr->build();

    return elemMatr;
};

ElementaryMatrixTemperatureRealPtr DiscreteComputation::getThermalTangentNonLinearVolumetricMatrix(
    const FieldOnNodesRealPtr temp_curr, const ASTERDOUBLE time_curr ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "MTAN_THER_SOURNL" );

    auto elemMatr =
        std::make_shared< ElementaryMatrixTemperatureReal >( _phys_problem->getModel(), option );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();

    // Prepare computing
    CalculPtr calcul = std::make_shared< Calcul >( option );
    calcul->setModel( currModel );

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        if ( load->hasLoadField( "SOUNL" ) ) {
            calcul->clearInputs();
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPEI", temp_curr );

            calcul->addInputField( "PSOURNL", load->getConstantLoadField( "SOUNL" ) );

            calcul->addHHOField( currModel );

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PMATTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PMATTTR" ) ) {
                elemMatr->addElementaryTerm( calcul->getOutputElementaryTermReal( "PMATTTR" ) );
            }
        }
    }

    elemMatr->build();
    return elemMatr;
};
