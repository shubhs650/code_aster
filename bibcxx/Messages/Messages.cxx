/**
 * @file Message.cxx
 * @brief Fichier entete de la class Messages
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

#include "Messages/Messages.h"

#include "aster_fort_utils.h"

#include <algorithm>

void UTMESS( const std::string &typm, const std::string &idmess, const VectorString &vk,
             const VectorLong &vi, const VectorReal &vr ) {
    if ( typm == "A" || typm == "I" ) {
        std::string typm2( typm ), idmess2( idmess );
        ASTERINTEGER nk( vk.size() ), ni( vi.size() ), nr( vr.size() );
        nk = std::max( nk, (ASTERINTEGER)1 );
        ni = std::max( ni, (ASTERINTEGER)1 );
        nr = std::max( nr, (ASTERINTEGER)1 );
        ASTERINTEGER nexc = 0;
        char *valk;
        char *fname;
        fname = MakeBlankFStr( 1 );
        valk = MakeTabFStr( nk, VALK_SIZE );
        for ( int i = 0; i < vk.size(); ++i ) {
            SetTabFStr( valk, i, vk[i].data(), VALK_SIZE );
        }
        ASTERINTEGER vali[ni];
        for ( int i = 0; i < vi.size(); ++i ) {
            vali[i] = vi[i];
        }
        ASTERDOUBLE valr[nr];
        for ( int i = 0; i < vr.size(); ++i ) {
            valr[i] = vr[i];
        }
        CALL_UTMESS_CORE( typm2.data(), idmess2.data(), &nk, valk, &ni, vali, &nr, valr, &nexc,
                          fname );
        FreeStr( valk );
        FreeStr( fname );
    } else {
        raiseAsterError( idmess, vk, vi, vr );
    }
}

void UTMESS( char *typm, char *idmess ) { UTMESS( std::string( typm ), std::string( idmess ) ); }
void UTMESS( const char *typm, const char *idmess ) { UTMESS( (char *)typm, (char *)idmess ); }
