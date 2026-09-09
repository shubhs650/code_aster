/**
 * @file TimesList.cxx
 * @brief Implementation de TimesList
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

#include "Steppers/TimesList.h"

bool TimesList::setValues( const VectorReal &values ) {
    _values->clear();
    _values->reserve( values.size() );
    _values->updateValuePointer();

    ASTERDOUBLE previous = -1.e300, minstep = 1.e300, maxstep = 0.;
    for ( VectorRealCIter tmp = values.begin(); tmp != values.end(); ++tmp ) {
        _values->push_back( *tmp );
        const ASTERDOUBLE &curVal = *tmp;
        if ( previous >= curVal )
            throw std::runtime_error( "Time function not strictly increasing" );
        if ( *tmp - previous < minstep )
            minstep = *tmp - previous;
        previous = *tmp;
    }

    _infor->updateValuePointer();
    // 0: MANUEL / AUTO
    // 1: PAS_MINI / 1.d-12
    ( *_infor )[2] = values.back() - values.front();
    // 3: NB_PAS_MAXI
    ( *_infor )[4] = minstep;
    ( *_infor )[7] = values.size();

    return true;
};
