/**
 * @file ParallelDOFNumbering.cxx
 * @brief Implementation de ParallelDOFNumbering
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

#include "Meshes/Joints.h"

#ifdef ASTER_HAVE_MPI

#include "aster_mpi.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/Exceptions.h"

#include <algorithm>
#include <chrono>

struct HashStruct {
    std::string nameRank;
    std::string timeString;
    bool operator==( const HashStruct & ) const = default;
};

template <>
struct std::hash< HashStruct > {
    std::size_t operator()( const HashStruct &s ) const noexcept {
        std::size_t h1 = std::hash< std::string > {}( s.nameRank );
        std::size_t h2 = std::hash< std::string > {}( s.timeString );
        return h1 ^ ( h2 << 1 );
    }
};

Joints::Joints() : Joints( DataStructureNaming::getNewName() ) {};

Joints::Joints( const std::string name )
    : DataStructure( name, 19, "DOMJOINTS" ),
      _layer( JeveuxVectorLong( getName() + ".NBLG" ) ),
      _domj( JeveuxVectorLong( getName() + ".DOMJ" ) ),
      _send( JeveuxCollectionLong( getName() + ".SEND" ) ),
      _recv( JeveuxCollectionLong( getName() + ".RECV" ) ),
      _procGroupIds( JeveuxVectorShort( getName() + ".PGID" ) ),
      _groupComm( JeveuxVectorLong( getName() + ".GCOM" ) ),
      _hashContainer( JeveuxVectorChar24( getName() + ".HASH" ) ) {
    buildGroup();
};

void Joints::buildGroup() {
    _procGroupIds->deallocate();
    _groupComm->deallocate();
    _hashContainer->deallocate();
    const auto aster_comm = aster_get_comm_world();
    const MPI_Comm comm = aster_comm->id;
    int rank, size;
    aster_get_mpi_info( aster_comm, &rank, &size );
    VectorInt toGather, vecOut, idsInParent, groupProcIds;
    if ( _domj->exists() )
        toGather.push_back( rank );
    else
        toGather.push_back( -1 );
    AsterMPI::all_gather( toGather, vecOut );
    int cmpt = 0;
    for ( auto id : vecOut ) {
        if ( id != -1 ) {
            idsInParent.push_back( cmpt );
            groupProcIds.push_back( id );
            ++cmpt;
        } else {
            idsInParent.push_back( -1 );
        }
    }
    *_procGroupIds = idsInParent;

    if ( groupProcIds.size() == 0 ) {
        groupProcIds.push_back( rank );
    }
    if ( _procGroupIds->exists() ) {
        _jointsMPIGroup = MPIGroupPtr( new MPIGroup() );
        _jointsMPIGroup->buildFromProcsVector( comm, groupProcIds );
        _groupComm->allocate( 1 );
        std::string uniqueHash( 24, ' ' );
        const auto grpSize = groupProcIds.size();
        const auto minRank = *std::min_element( groupProcIds.begin(), groupProcIds.end() );

        const auto now = std::chrono::system_clock::now().time_since_epoch();
        const auto milliseconds =
            std::chrono::duration_cast< std::chrono::milliseconds >( now ).count();
        const auto timeStr( std::to_string( milliseconds ) );
        HashStruct obj { getName() + std::to_string( rank ), timeStr };
        // build unique identifier
        uniqueHash = std::to_string( std::hash< HashStruct > {}( obj ) ) + std::string( 24, ' ' );
        uniqueHash = uniqueHash.substr( 0, 24 );

        const auto rankSearch = std::find( groupProcIds.begin(), groupProcIds.end(), rank );
        bool rankFind = ( rankSearch != groupProcIds.end() ) ? true : false;
        if ( grpSize > 1 && rankFind ) {
            const aster_comm_t *curACommT = _jointsMPIGroup->getAsterCommunicator();
            aster_comm_t *curAComm = const_cast< aster_comm_t * >( curACommT );
            AsterMPI::bcast( uniqueHash, 0, curAComm );
        }
        _hashContainer->allocate( 1 );
        ( *_hashContainer )[0] = uniqueHash;
        const auto tmp = MPI_Comm_c2f( _jointsMPIGroup->getCommunicator() );
        ( *_groupComm )[0] = tmp;
    }
};

void Joints::setOppositeDomains( const VectorLong &oppositeDomains ) {
    if ( oppositeDomains.size() != 0 )
        ( *_domj ) = oppositeDomains;
    buildGroup();
};

const JeveuxVectorLong &Joints::getOppositeDomains() const {
    if ( _domj->exists() ) {
        _domj->updateValuePointer();
    }

    return _domj;
}

bool Joints::build() {
    _send->build();
    _recv->build();
    buildGroup();

    return true;
}

void Joints::setSendedElements( const VectorOfVectorsLong &send ) {
    if ( send.size() > 0 ) {
        _send->allocate( send.size() );

        ASTERINTEGER i = 1;
        for ( auto &send_i : send ) {
            if ( send_i.size() > 0 ) {
                auto obj = _send->allocateObject( i, send_i );
            }
            i++;
        }
    }
};

void Joints::setReceivedElements( const VectorOfVectorsLong &recv ) {
    if ( recv.size() > 0 ) {
        _recv->allocate( recv.size() );

        ASTERINTEGER i = 1;
        for ( auto &recv_i : recv ) {
            if ( recv_i.size() > 0 ) {
                auto obj = _recv->allocateObject( i, recv_i );
            }
            i++;
        }
    }
};

VectorLong Joints::getSendedElements( const ASTERINTEGER &id ) const {
    const auto &obj = ( *_send )[id + 1];
    obj->updateValuePointer();
    return obj->toVector();
};

VectorLong Joints::getReceivedElements( const ASTERINTEGER &id ) const {
    const auto &obj = ( *_recv )[id + 1];
    obj->updateValuePointer();
    return obj->toVector();
};

void Joints::setNumberOfGhostLayer( const ASTERINTEGER &nb_layer ) {
    if ( !_layer->exists() ) {
        _layer->allocate( 1 );
    }
    _layer->updateValuePointer();

    if ( nb_layer < 1 ) {
        raiseAsterError( "ValueError: ghost layer number must be at least 1." );
    }

    ( *_layer )[0] = nb_layer;
};

ASTERINTEGER Joints::getNumberOfGhostLayer( void ) const {
    _layer->updateValuePointer();
    return ( *_layer )[0];
};

#endif /* ASTER_HAVE_MPI */
