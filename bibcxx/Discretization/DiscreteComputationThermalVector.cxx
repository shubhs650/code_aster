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

#include "DataFields/FieldOnCellsBuilder.h"
#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalNeumannForces( const ASTERDOUBLE time_curr,
                                              const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect =
        std::make_shared< ElementaryVectorTemperatureReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );
    auto isXfem = currModel->existsXfem();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        // Termes FLUX XYZ
        if ( load->hasLoadField( "FLURE" ) ) {
            auto flow_xyz_field = load->getConstantLoadField( "FLURE" );
            calcul->setOption( "CHAR_THER_FLUN_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PFLUXNR", flow_xyz_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes FLUX NORM
        if ( load->hasLoadField( "FLUR2" ) ) {
            auto flow_nor_field = load->getConstantLoadField( "FLUR2" );
            calcul->setOption( "CHAR_THER_FLUX_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PFLUXVR", flow_nor_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        // Termes FLUX XYZ
        if ( load->hasLoadField( "FLURE" ) ) {
            auto flow_xyz_field = load->getConstantLoadField( "FLURE" );
            calcul->setOption( "CHAR_THER_FLUN_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PFLUXNF", flow_xyz_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes FLUX NORM
        if ( load->hasLoadField( "FLUR2" ) ) {
            auto flow_nor_field = load->getConstantLoadField( "FLUR2" );
            calcul->setOption( "CHAR_THER_FLUX_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PFLUXVF", flow_nor_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalVolumetricForces( const ASTERDOUBLE time_curr,
                                                 const FieldOnCellsRealPtr varc_curr,
                                                 const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect =
        std::make_shared< ElementaryVectorTemperatureReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );
    auto isXfem = currModel->existsXfem();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {

        // Termes SOURCE
        if ( load->hasLoadField( "SOURE" ) ) {
            auto source_field = load->getConstantLoadField( "SOURE" );
            calcul->setOption( "CHAR_THER_SOUR_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PSOURCR", source_field );

            if ( currMater && currMater->hasExternalStateVariable() ) {
                if ( !varc_curr || !varc_curr->exists() ) {
                    raiseAsterError( "External state variables are needed but not given" );
                }
                calcul->addInputField( "PVARCPR", varc_curr );
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }
        // Termes SOURCE CALCULEE
        if ( load->hasLoadField( "SOURC" ) ) {
            auto computed_source_field = load->getLoadField( "SOURC" );
            calcul->setOption( "CHAR_THER_SOUR_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PSOURCR", computed_source_field );

            if ( currMater && currMater->hasExternalStateVariable() ) {
                if ( !varc_curr || !varc_curr->exists() ) {
                    raiseAsterError( "External state variables are needed but not given" );
                }
                calcul->addInputField( "PVARCPR", varc_curr );
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes PRE_GRAD_TEMP
        if ( load->hasLoadField( "GRAIN" ) ) {
            auto pregrad_field = load->getConstantLoadField( "GRAIN" );
            calcul->setOption( "CHAR_THER_GRAI_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
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
            calcul->addInputField( "PGRAINR", pregrad_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();

        if ( load->hasLoadField( "SOURE" ) ) {
            auto source_field = load->getConstantLoadField( "SOURE" );
            calcul->setOption( "CHAR_THER_SOUR_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PSOURCF", source_field );

            if ( currMater && currMater->hasExternalStateVariable() ) {
                if ( !varc_curr || !varc_curr->exists() ) {
                    raiseAsterError( "External state variables are needed but not given" );
                }
                calcul->addInputField( "PVARCPR", varc_curr );
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes PRE_GRAD_TEMP
        if ( load->hasLoadField( "GRAIN" ) ) {
            auto pregrad_field = load->getConstantLoadField( "GRAIN" );
            calcul->setOption( "CHAR_THER_GRAI_F" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );

            if ( currMater ) {
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

                if ( currMater->hasExternalStateVariable() ) {
                    if ( !varc_curr || !varc_curr->exists() ) {
                        raiseAsterError( "External state variables are needed but not given" );
                    }
                    calcul->addInputField( "PVARCPR", varc_curr );
                }
            }
            calcul->addInputField( "PGRAINF", pregrad_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalExchangeForces( const FieldOnNodesRealPtr temp_curr,
                                               const ASTERDOUBLE time_curr,
                                               const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect =
        std::make_shared< ElementaryVectorTemperatureReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );
    auto isXfem = currModel->existsXfem();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();

        // Termes ECHANGE
        if ( load->hasLoadField( "COEFH" ) && load->hasLoadField( "T_EXT" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            auto ext_temp_field = load->getConstantLoadField( "T_EXT" );

            calcul->setOption( "CHAR_THER_ECHA_R" );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PCOEFHR", exchange_field );
            calcul->addInputField( "PT_EXTR", ext_temp_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes ECHANGE_PAROI
        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "CHAR_THER_PARO_R" );
            calcul->clearInputs();
            if ( isXfem ) {
                XfemModelPtr currXfemModel = currModel->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( model_FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PHECHPR", wall_exchange_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        iload++;
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();
        // Termes ECHANGE
        if ( load->hasLoadField( "COEFH" ) && load->hasLoadField( "T_EXT" ) ) {
            auto exchange_field = load->getConstantLoadField( "COEFH" );
            auto ext_temp_field = load->getConstantLoadField( "T_EXT" );

            calcul->setOption( "CHAR_THER_ECHA_F" );
            calcul->clearInputs();
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->setFiniteElementDescriptor( model_FEDesc );
            calcul->addInputField( "PCOEFHF", exchange_field );
            calcul->addInputField( "PT_EXTF", ext_temp_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }

        // Termes ECHANGE_PAROI
        if ( load->hasLoadField( "HECHP" ) ) {
            auto wall_exchange_field = load->getConstantLoadField( "HECHP" );
            calcul->setOption( "CHAR_THER_PARO_F" );
            calcul->clearInputs();
            if ( isXfem ) {
                XfemModelPtr currXfemModel = currModel->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->setFiniteElementDescriptor( model_FEDesc );
            } else {
                calcul->setFiniteElementDescriptor( load_FEDesc );
            }
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addHHOField( currModel );
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PHECHPF", wall_exchange_field );
            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             iload );
            }
        }
        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalNonLinearNeumannForces( const FieldOnNodesRealPtr temp_curr,
                                                       const ASTERDOUBLE time_curr,
                                                       const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect =
        std::make_shared< ElementaryVectorTemperatureReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto impl = [&]( auto load, const ASTERINTEGER &load_i, const std::string &option,
                     const std::string &name, const std::string &param,
                     const FiniteElementDescriptorPtr FED ) {
        if ( load->hasLoadField( name ) ) {
            calcul->setOption( option );
            calcul->setFiniteElementDescriptor( FED );

            calcul->clearInputs();
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );
            calcul->addHHOField( currModel );

            calcul->addInputField( param, load->getConstantLoadField( name ) );

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             load_i );
            }
        }
    };

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        impl( load, iload, "CHAR_THER_RAYO_R", "RAYO", "PRAYONR", model_FEDesc );

        iload++;
    }

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        impl( load, iload, "CHAR_THER_FLUNL", "FLUNL", "PFLUXNL", model_FEDesc );
        impl( load, iload, "CHAR_THER_RAYO_F", "RAYO", "PRAYONF", model_FEDesc );

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalNonLinearVolumetricForces( const FieldOnNodesRealPtr temp_curr,
                                                          const ASTERDOUBLE time_curr,
                                                          const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect =
        std::make_shared< ElementaryVectorTemperatureReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto impl = [&]( auto load, const ASTERINTEGER &load_i, const std::string &option,
                     const std::string &name, const std::string &param,
                     const FiniteElementDescriptorPtr FED ) {
        if ( load->hasLoadField( name ) ) {
            calcul->setOption( option );
            calcul->setFiniteElementDescriptor( FED );

            calcul->clearInputs();
            calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PTEMPER", temp_curr );

            calcul->addHHOField( currModel );

            calcul->addInputField( param, load->getConstantLoadField( name ) );

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                             load_i );
            }
        }
    };

    auto therLoadFunc = listOfLoads->getThermalLoadsFunction();
    for ( const auto &load : therLoadFunc ) {
        impl( load, iload, "CHAR_THER_SOURNL", "SOUNL", "PSOURNL", model_FEDesc );

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getTransientThermalLoadForces( const ASTERDOUBLE time_curr,
                                                    const FieldOnNodesRealPtr temp_prev,
                                                    const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect =
        std::make_shared< ElementaryVectorTemperatureReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );

    auto therLoadReal = listOfLoads->getThermalLoadsReal();
    for ( const auto &load : therLoadReal ) {
        if ( load->hasLoadResult() ) {
            if ( load->canInterpolateLoadResult( "FLUN" ) ) {
                FieldOnCellsRealPtr evol_flow_xyz_field =
                    load->interpolateLoadResult( "FLUN", time_curr );

                calcul->setOption( "CHAR_THER_FLUN_R" );
                calcul->setFiniteElementDescriptor( model_FEDesc );
                calcul->clearInputs();
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
                calcul->addInputField( "PFLUXNR", evol_flow_xyz_field );
                calcul->clearOutputs();
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                }

            } else {
                AS_ASSERT( temp_prev && temp_prev->exists() );
                FieldOnCellsRealPtr evol_exchange_field =
                    load->interpolateLoadResult( "COEF_H", time_curr );
                FieldOnCellsRealPtr evol_ext_temp_field =
                    load->interpolateLoadResult( "T_EXT", time_curr );

                calcul->setOption( "CHAR_THER_ECHA_R" );
                calcul->setFiniteElementDescriptor( model_FEDesc );
                calcul->clearInputs();
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addInputField( "PTEMPER", temp_prev );
                calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
                calcul->addInputField( "PCOEFHR", evol_exchange_field );
                calcul->addInputField( "PT_EXTR", evol_ext_temp_field );
                calcul->clearOutputs();
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                }
            }
        }

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

/** @brief Compute AFFE_CHAR_THER TEMP_IMPO */
std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getThermalImposedDualBC( const ASTERDOUBLE time_curr,
                                              const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    auto elemVect =
        std::make_shared< ElementaryVectorTemperatureReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_THER" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl = [&]( auto loads, bool real ) {
        std::string name;
        if ( real ) {
            calcul->setOption( "THER_DDLI_R" );
            name = "PDDLIMR";
        } else {
            calcul->setOption( "THER_DDLI_F" );
            name = "PDDLIMF";
        }
        for ( const auto &load : loads ) {
            auto load_FEDesc = load->getFiniteElementDescriptor();
            auto impo_field = load->getImposedField();
            if ( impo_field && impo_field->exists() && load_FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( load_FEDesc );
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addInputField( name, impo_field );
                calcul->addTimeField( "PINSTR", time_curr, 0.0, -1.0 );
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                }
            }
            iload++;
        }
    };

    impl( listOfLoads->getThermalLoadsReal(), true );
    impl( listOfLoads->getThermalLoadsFunction(), false );

#ifdef ASTER_HAVE_MPI
    impl( listOfLoads->getParallelThermalLoadsReal(), true );
    impl( listOfLoads->getParallelThermalLoadsFunction(), false );
#endif

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
}

/** @brief Compute CHAR_THER_EVOL */
FieldOnNodesRealPtr DiscreteComputation::getTransientThermalForces(
    const ASTERDOUBLE time_curr, const ASTERDOUBLE time_step, const ASTERDOUBLE theta,
    const FieldOnNodesRealPtr previousPrimalField, const FieldOnCellsRealPtr varc_curr ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    AS_ASSERT( previousPrimalField && previousPrimalField->exists() );

    auto elemVect = std::make_shared< ElementaryVectorReal >( _phys_problem->getModel() );

    // Setup
    const std::string calcul_option( "CHAR_THER_EVOL" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setModel( currModel );
    calcul->clearInputs();
    calcul->clearOutputs();

    // Add input fields
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addTimeField( "PINSTR", time_curr, time_step, theta );
    calcul->addInputField( "PTEMPER", previousPrimalField );
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
    calcul->addXFEMField( currModel );

    calcul->addHHOField( currModel );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PVECTTR", std::make_shared< ElementaryTermReal >() );

    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ) );
    };

    elemVect->build();
    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}

/**
 * @brief Compute elementary forces for internal forces (RAPH_THER)
 */
std::tuple< ASTERINTEGER, FieldOnCellsRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getInternalThermalForces( const FieldOnNodesRealPtr temp_prev,
                                               const FieldOnNodesRealPtr temp_step,
                                               const FieldOnCellsRealPtr varc_curr,
                                               const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "RAPH_THER" );

    auto elemVect = std::make_shared< ElementaryVectorReal >( _phys_problem->getModel() );

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
    auto FEDesc = calcul->getFiniteElementDescriptor();

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

    // Add Thermal Field
    auto temp_curr = std::make_shared< FieldOnNodesReal >( *temp_prev + *temp_step );
    calcul->addInputField( "PTEMPEI", temp_curr );

    // TODO:
    // calcul->addInputField( "PTMPCHF", dry_curr );

    // Create output fields
    auto flux_curr = FieldOnCellsPtrBuilder< ASTERDOUBLE >( FEDesc, "ELGA", "FLUX_R" );
    calcul->addOutputField( "PFLUXPR", flux_curr );

    // Add output elementary terms
    calcul->addOutputElementaryTerm( "PRESIDU", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PRESIDU" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PRESIDU" ) );
    };

    elemVect->build();
    auto internalForces = elemVect->assemble( _phys_problem->getDOFNumbering() );

    return std::make_tuple( 0, flux_curr, internalForces );
}

/**
 * @brief Compute elementary forces for capacitiy forces (MASS_THER_RESI)
 */
FieldOnNodesRealPtr DiscreteComputation::getNonLinearCapacityForces(
    const FieldOnNodesRealPtr temp_prev, const FieldOnNodesRealPtr temp_step,
    const FieldOnCellsRealPtr varc_curr, const VectorString &groupOfCells ) const {
    AS_ASSERT( _phys_problem->getModel()->isThermal() );
    const std::string option( "MASS_THER_RESI" );

    auto elemVect = std::make_shared< ElementaryVectorReal >( _phys_problem->getModel() );

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
    calcul->addOutputElementaryTerm( "PRESIDU", std::make_shared< ElementaryTermReal >() );

    // Compute elementary matrices for mass
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PRESIDU" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PRESIDU" ) );
    };

    elemVect->build();
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
}

FieldOnNodesRealPtr DiscreteComputation::getDualTemperature( FieldOnNodesRealPtr temp_curr,
                                                             ASTERDOUBLE scaling ) const {

    auto elemVect = std::make_shared< ElementaryVectorReal >( _phys_problem->getModel() );

    if ( !_phys_problem->isThermal() ) {
        AS_ABORT( "Not implemented" );
    };

    // Prepare loads
    const auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( "THER_BU_R" );

    // ConstantField for scaling
    auto cf_scaling = std::make_shared< ConstantFieldOnCellsReal >( _phys_problem->getMesh() );
    cf_scaling->allocate( "NEUT_R" );
    ConstantFieldOnZone a( _phys_problem->getMesh() );
    ConstantFieldValues< ASTERDOUBLE > b( { "X1" }, { scaling } );
    cf_scaling->setValueOnZone( a, b );

    auto impl = [&]( auto loads ) {
        for ( const auto &load : loads ) {
            auto FEDesc = load->getFiniteElementDescriptor();
            auto field = load->getMultiplicativeField();
            if ( field && field->exists() && FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( FEDesc );
                calcul->addInputField( "PDDLMUR", field );
                calcul->addInputField( "PDDLIMR", temp_curr );
                calcul->addInputField( "PALPHAR", cf_scaling );
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ) );
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

    // Assemble
    elemVect->build();

    FieldOnNodesRealPtr buth;
    if ( elemVect->hasElementaryTerm() ) {
        buth = elemVect->assemble( _phys_problem->getDOFNumbering() );
    } else {
        buth = std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
        buth->setValues( 0.0 );
        buth->build();
    }

    if ( _phys_problem->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( buth->getName() );

    return buth;
};

FieldOnNodesRealPtr DiscreteComputation::dualThermalVector( FieldOnNodesRealPtr lagr_curr ) const {

    AS_ASSERT( _phys_problem->getModel()->isThermal() );

    // Init
    ASTERINTEGER iload = 1;

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currMater = _phys_problem->getMaterialField();
    auto currListOfLoads = _phys_problem->getListOfLoads();

    // Select option to compute
    std::string option = "THER_BTLA_R";

    // Create elementary vector
    auto elemVect = std::make_shared< ElementaryVectorReal >( currModel );

    // Setup
    const std::string calcul_option( "CHAR_THER" );

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );

    auto impl_load = [&]( auto loads ) {
        for ( const auto &load : loads ) {
            auto load_FEDesc = load->getFiniteElementDescriptor();
            auto mult_field = load->getMultiplicativeField();
            if ( mult_field && mult_field->exists() && load_FEDesc ) {

                // Add input fields
                calcul->clearInputs();
                calcul->addInputField( "PLAGRAR", lagr_curr );
                calcul->addInputField( "PDDLMUR", mult_field );

                // Add output fields
                calcul->clearOutputs();
                calcul->addOutputElementaryTerm( "PVECTTR",
                                                 std::make_shared< ElementaryTermReal >() );

                // Compute
                calcul->setFiniteElementDescriptor( load_FEDesc );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTTR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTTR" ),
                                                 iload );
                }
            }
            iload++;
        }
    };

    impl_load( currListOfLoads->getThermalLoadsReal() );
    impl_load( currListOfLoads->getThermalLoadsFunction() );

#ifdef ASTER_HAVE_MPI
    impl_load( currListOfLoads->getParallelThermalLoadsReal() );
    impl_load( currListOfLoads->getParallelThermalLoadsFunction() );
#endif

    // Construct elementary vector
    auto FEDs = _phys_problem->getListOfLoads()->getFiniteElementDescriptors();
    FEDs.push_back( _phys_problem->getModel()->getFiniteElementDescriptor() );
    elemVect->build( FEDs );

    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
};
