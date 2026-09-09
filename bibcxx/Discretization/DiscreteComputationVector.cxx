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
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"
#include "Supervis/CommandSyntax.h"
#include "Utilities/Tools.h"

FieldOnNodesRealPtr DiscreteComputation::getDualForces( FieldOnNodesRealPtr lagr_curr ) const {

    if ( _phys_problem->getModel()->isMechanical() ) {
        return dualMechanicalVector( lagr_curr );
    } else if ( _phys_problem->getModel()->isThermal() ) {
        return dualThermalVector( lagr_curr );
    } else {
        AS_ABORT( "Should not be here" );
    }
};

FieldOnNodesRealPtr DiscreteComputation::getDualPrimal( FieldOnNodesRealPtr primal_curr,
                                                        ASTERDOUBLE scaling ) const {
    if ( _phys_problem->isMechanical() ) {
        return this->getDualDisplacement( primal_curr, scaling );
    } else if ( _phys_problem->isThermal() ) {
        return this->getDualTemperature( primal_curr, scaling );
    } else {

        AS_ABORT( "Not implemented" );
    }
};

template < typename T >
std::shared_ptr< FieldOnNodes< T > >
DiscreteComputation::_getDirichletBC( const ASTERDOUBLE time_curr ) const {

    auto dofNume = _phys_problem->getDOFNumbering();
    auto vectAsse = std::make_shared< FieldOnNodes< T > >( dofNume );

    // Prepare loads
    const auto &_listOfLoads = _phys_problem->getListOfLoads();
    if ( !_listOfLoads->isBuilt() )
        _listOfLoads->build( _phys_problem->getModel() );

    JeveuxVectorChar24 listOfLoadsList = _listOfLoads->getLoadNames();
    JeveuxVectorLong listOfLoadsInfo = _listOfLoads->getInformationVector();
    JeveuxVectorChar24 listOfLoadsFunc = _listOfLoads->getListOfFunctions();
    std::string nameLcha = ljust( listOfLoadsList->getName(), 24 );
    std::string nameInfc = ljust( listOfLoadsInfo->getName(), 24 );
    std::string nameFcha = ljust( listOfLoadsFunc->getName(), 24 );

    // Get JEVEUX names of objects to call Fortran
    std::string vectAsseName = vectAsse->getName();
    std::string dofNumName = dofNume->getName();
    std::string loadName = _listOfLoads->getName();
    std::string modelName = _phys_problem->getModel()->getName();
    std::string base( "G" );

    // Wrapper FORTRAN
    CALLO_ASCAVC_WRAP( modelName, loadName, dofNumName, &time_curr, vectAsseName, base );

    // Construct vect_asse object
    vectAsse->build();

    // Assemble
    return vectAsse;
};

FieldOnNodesRealPtr DiscreteComputation::_getDirichletDifferentialBC() const {

    auto dofNume = _phys_problem->getDOFNumbering();
    auto equaNume = dofNume->getEquationNumbering();
    equaNume->updateValuePointers();
    const auto &_listOfLoads = _phys_problem->getListOfLoads();
    if ( !_listOfLoads->isBuilt() )
        _listOfLoads->build( _phys_problem->getModel() );
    auto didi_depl = _listOfLoads->getDifferentialDisplacement();

    auto vectAsse = std::make_shared< FieldOnNodesReal >( dofNume );
    vectAsse->setValues( 0.0 );

    ASTERINTEGER iload = 1;
    auto types = _listOfLoads->getListOfDiriTyp();
    for ( const auto &diriBC : _listOfLoads->getDirichletBCs() ) {
        if ( types[iload - 1] == "DIDI" ) {
            // add didi_depl on char_cine dofs (voir numchc.F90)
            auto intParam = diriBC->getIntParam();
            intParam->updateValuePointer();
            ASTERINTEGER nbloc = ( *intParam )[0];
            for ( ASTERINTEGER ibloc = 0; ibloc < nbloc; ibloc++ ) {
                ASTERINTEGER node = ( *intParam )[1 + 3 * ( ibloc )] - 1;
                ASTERINTEGER cmp = ( *intParam )[1 + 3 * ( ibloc ) + 1] - 1;
                ASTERINTEGER dof = equaNume->getDOFsFromNode( node )[cmp];
                ( *vectAsse )[dof] = ( *didi_depl )[dof];
            }
        }
        iload++;
    }

    return vectAsse;
}

FieldOnNodesRealPtr
DiscreteComputation::getMechanicalDirichletBC( const ASTERDOUBLE time_curr ) const {
    auto vectAsse = this->_getDirichletBC< ASTERDOUBLE >( time_curr );
    if ( _phys_problem->getListOfLoads()->hasDifferentialDirichletBC() )
        ( *vectAsse ) += ( *this->_getDirichletDifferentialBC() );
    return vectAsse;
}

FieldOnNodesRealPtr
DiscreteComputation::getThermalDirichletBC( const ASTERDOUBLE time_curr ) const {
    return this->_getDirichletBC< ASTERDOUBLE >( time_curr );
}

FieldOnNodesComplexPtr
DiscreteComputation::getAcousticDirichletBC( const ASTERDOUBLE time_curr ) const {
    CommandSyntax cmdSt( "CALC_CHAR_CINE" );
    SyntaxMapContainer dict;
    ListSyntaxMapContainer listeExcit;
    listeExcit.push_back( dict );
    SyntaxMapContainer dict2;

    dict.container["EXCIT"] = listeExcit;
    cmdSt.define( dict );

    return this->_getDirichletBC< ASTERCOMPLEX >( time_curr );
}

FieldOnNodesRealPtr
DiscreteComputation::getIncrementalDirichletBC( const ASTERDOUBLE &time_curr,
                                                const FieldOnNodesRealPtr disp_curr ) const {
    auto dofNume = _phys_problem->getDOFNumbering();

    const auto listOfLoads = _phys_problem->getListOfLoads();
    if ( listOfLoads->hasDirichletBC() ) {
        auto diri_curr = this->_getDirichletBC< ASTERDOUBLE >( time_curr );
        if ( listOfLoads->hasDifferentialDirichletBC() )
            ( *diri_curr ) += ( *this->_getDirichletDifferentialBC() );

        auto diri_impo = *( diri_curr ) - *( disp_curr );
        diri_impo.updateValuePointers();
        // Set to zero terms not imposed
        auto eliminatedDofs = _phys_problem->getDirichletBCDOFs();
        auto nbElimination = eliminatedDofs.size();

        for ( ASTERINTEGER ieq = 0; ieq < nbElimination; ieq++ ) {
            if ( eliminatedDofs[ieq] == 0 )
                diri_impo[ieq] = 0.0;
        }

        return std::make_shared< FieldOnNodesReal >( diri_impo );
    }

    // Construct vect_asse object
    FieldOnNodesRealPtr vectAsse = std::make_shared< FieldOnNodesReal >( dofNume );
    vectAsse->setValues( 0.0 );

    return vectAsse;
};

std::variant< ElementaryVectorRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getExternalStateVariablesForces(
    const ASTERDOUBLE time_curr, const FieldOnCellsRealPtr varc_curr,
    const FieldOnCellsRealPtr varc_prev, const FieldOnCellsRealPtr vari_curr,
    const FieldOnCellsRealPtr stress_prev, const ASTERINTEGER mode_fourier, const bool assembly,
    const FieldOnCellsLongPtr maskField ) const {

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currExternVarRefe = _phys_problem->getReferenceExternalStateVariables();
    auto currBehav = _phys_problem->getBehaviourProperty();

    // Some checks
    AS_ASSERT( currMater );
    AS_ASSERT( currMater->hasExternalStateVariableForLoad() );
    if ( currMater->hasExternalStateVariableWithReference() ) {
        AS_ASSERT( currExternVarRefe );
    }

    if ( !varc_curr || !varc_curr->exists() ) {
        raiseAsterError( "External state variables are needed but not given" );
    }

    // Main object
    CalculPtr calcul = std::make_unique< Calcul >( "CHAR_VARC" );

    // Create specific output field for XFEM
    FieldOnCellsRealPtr sigmXfem;
    if ( currModel->existsXfem() ) {
        const std::string option = "SIEF_ELGA";
        const std::string paraName = "PCONTRR";
        sigmXfem = std::make_shared< FieldOnCellsReal >( currModel, option, paraName );
    }

    // Create elementary vectors
    auto elemVect = std::make_shared< ElementaryVectorReal >( _phys_problem->getModel() );

    int nbExternVar = static_cast< int >( externVarEnumInt::NumberOfExternVarTypes );
    for ( auto iExternVar = 0; iExternVar < nbExternVar; iExternVar++ ) {
        externVarEnumInt numeExternVar = static_cast< externVarEnumInt >( iExternVar );
        if ( currMater->hasExternalStateVariable( numeExternVar ) &&
             ExternalVariableTraits::externVarHasStrain( numeExternVar ) ) {
            const auto option = ExternalVariableTraits::getExternVarOption( numeExternVar );
            calcul->setOption( option );
            calcul->setModel( currModel );
            calcul->clearInputs();
            calcul->clearOutputs();

            // Add input fields
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
            if ( currBehav ) {
                calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
            }
            calcul->addFourierModeField( mode_fourier );
            calcul->addInputField( "PVARCPR", varc_curr );
            if ( varc_prev ) {
                calcul->addInputField( "PVARCMR", varc_prev );
            }
            if ( vari_curr ) {
                calcul->addInputField( "PVARIPR", vari_curr );
            }
            if ( stress_prev ) {
                calcul->addInputField( "PCONTMR", stress_prev );
            }
            if ( currElemChara ) {
                calcul->addElementaryCharacteristicsField( currElemChara );
                calcul->addOrientationField( currElemChara );
            }

            calcul->addTimeField( "PINSTR", time_curr );
            if ( currExternVarRefe ) {
                calcul->addInputField( "PVARCRR", currExternVarRefe );
            }
            if ( currModel->existsXfem() ) {
                XfemModelPtr currXfemModel = currModel->getXfemModel();
                calcul->addXFEMField( currXfemModel );
                calcul->addOutputField( "PCONTRT", sigmXfem );
            }
            calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ) );
            }
        }
    }
    // Build elementary vectors
    elemVect->build();

    if ( assembly ) { // Assemble vector
        if ( maskField ) {
            return elemVect->assembleWithMask( _phys_problem->getDOFNumbering(), maskField, 1 );
        } else {
            return elemVect->assemble( _phys_problem->getDOFNumbering() );
        }
    }

    return elemVect;
}
