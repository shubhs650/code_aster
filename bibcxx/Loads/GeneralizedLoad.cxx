/**
 * @file GeneralizedLoad.cxx
 * @brief Implementation de GeneralizedLoad
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

#include "Loads/GeneralizedLoad.h"

#include "LinearAlgebra/GeneralizedAssemblyVector.h"

const GeneralizedAssemblyVectorRealPtr GeneralizedLoad::getAssemblyVector() const {
    GeneralizedAssemblyVectorRealPtr assemblyVector( new GeneralizedAssemblyVectorReal() );
    assemblyVector->setGeneralizedDOFNumbering( _nume );
    auto basis = _nume->getModalBasisFromModeResult();
    assemblyVector->setModalBasis( basis );
    assemblyVector->allocate();
    return assemblyVector;
};

void GeneralizedLoad::setMPCs( std::vector< std::tuple< VectorInt, VectorReal, double > > &mpcs ) {
    for ( auto [nume_mode, coef_mult, coef_impo] : mpcs ) {
        if ( nume_mode.size() != coef_mult.size() )
            throw std::runtime_error( "Inconsistent size between NUME_MODE and COEF_MULT" );
    }
    _mpcs = mpcs;
}
