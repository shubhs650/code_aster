/**
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

#pragma once

#include "astercxx.h"

#include "aster_fort_ds.h"

#include "Behaviours/BehaviourProperty.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/SimpleFieldOnCells.h"
#include "Discretization/Calcul.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Modeling/Model.h"

/**
 * There is a separate file to construct some object to avoid circular inclusions.
 * Theses new methods are binding as contructor.
 */

/**
 * @brief Constructor for empty FieldOnCells with dynamic components
 * @param model model
 * @param behaviour Description of behaviour (for size of dynamic components as VARI_ELGA)
 * @param carael Description of elementary characteristics (for size of dynamic components as
 * VARI_ELGA)
 * @param typcham Type de champ à calculer
 */
template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
FieldOnCellsPtrBuilder( const FiniteElementDescriptorPtr FEDesc, const std::string &loc,
                        const std::string &quantity, const BehaviourPropertyPtr behaviour = nullptr,
                        const ElementaryCharacteristicsPtr carael = nullptr ) {
    auto cham_elem = std::make_shared< FieldOnCells< ValueType > >( FEDesc );
    std::string option;
    std::string nompar;
    bool hastrx = ( loc == "ELGA" && quantity == "STRX_R" );

    if ( loc == "ELGA" ) {
        option = "TOU_INI_ELGA";
    } else if ( loc == "ELNO" ) {
        option = "TOU_INI_ELNO";
    } else if ( loc == "ELEM" ) {
        option = "TOU_INI_ELEM";
    } else {
        option = loc;
    };

    if ( quantity[0] != 'P' || quantity == "PRES_R" ) {
        nompar = "P" + quantity;
    } else {
        nompar = quantity;
    };

    ASTERINTEGER iret = 0;

    std::string dcel = " ";

    if ( !hastrx && ( behaviour || carael ) ) {
        std::string carele = " ", comporName = " ";

        if ( behaviour ) {
            comporName = behaviour->getBehaviourField()->getName();
        }

        if ( carael ) {
            carele = carael->getName();
        }

        auto _DCEL =
            std::make_shared< SimpleFieldOnCellsLong >( cham_elem->getName(), FEDesc->getMesh() );
        CALLO_CESVAR( carele, comporName, cham_elem->getDescription()->getName(),
                      _DCEL->getName() );
        dcel = _DCEL->getName();
        cham_elem->setExtentedInformations( _DCEL );
    }

    if ( hastrx ) {
        AS_ASSERT( carael );
        CalculPtr calcul = std::make_unique< Calcul >( "INI_STRX" );
        calcul->setModel( carael->getModel() );
        calcul->addElementaryCharacteristicsField( carael );
        calcul->addOrientationField( carael );
        calcul->addOutputField( "PSTRX_R", cham_elem );
        calcul->compute();
    } else {
        CALLO_ALCHML( cham_elem->getDescription()->getName(), option, nompar,
                      JeveuxMemoryTypesNames[Permanent], cham_elem->getName(), &iret, dcel );
        AS_ASSERT( iret == 0 );
    }
    cham_elem->updateValuePointers();

    return cham_elem;
};

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
FieldOnCellsPtrBuilder( const FiniteElementDescriptorPtr FEDesc, const std::string &loc,
                        const std::string &quantity, const BehaviourPropertyPtr behaviour ) {
    return FieldOnCellsPtrBuilder< ValueType >( FEDesc, loc, quantity, behaviour, nullptr );
};

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
FieldOnCellsPtrBuilder( const FiniteElementDescriptorPtr FEDesc, const std::string &loc,
                        const std::string &quantity, const ElementaryCharacteristicsPtr carael ) {
    return FieldOnCellsPtrBuilder< ValueType >( FEDesc, loc, quantity, nullptr, carael );
};

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
FieldOnCellsPtrBuilder( const ModelPtr model, const std::string &loc, const std::string &quantity,
                        const BehaviourPropertyPtr behaviour = nullptr,
                        const ElementaryCharacteristicsPtr carael = nullptr ) {
    return FieldOnCellsPtrBuilder< ValueType >( model->getFiniteElementDescriptor(), loc, quantity,
                                                behaviour, carael );
};

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
FieldOnCellsPtrBuilder( const ModelPtr model, const std::string &loc, const std::string &quantity,
                        const BehaviourPropertyPtr behaviour ) {
    return FieldOnCellsPtrBuilder< ValueType >( model, loc, quantity, behaviour, nullptr );
};

template < typename ValueType >
std::shared_ptr< FieldOnCells< ValueType > >
FieldOnCellsPtrBuilder( const ModelPtr model, const std::string &loc, const std::string &quantity,
                        const ElementaryCharacteristicsPtr carael ) {
    return FieldOnCellsPtrBuilder< ValueType >( model, loc, quantity, nullptr, carael );
};
