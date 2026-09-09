#ifndef ASTERMPI_H_
#define ASTERMPI_H_

/**
 * @file AsterMPI.h
 * @brief Fichier entete contenant des utilitaires de manipulation de containers
 * STL en parallèle
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

#include "aster_mpi.h"

#include <algorithm>
#include <numeric>
#include <set>

#ifdef ASTER_HAVE_MPI

#include "astercxx.h"

#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/JeveuxVector.h"

#endif /* ASTER_HAVE_MPI */

/** @brief Get MPI number of procs */
int getMPISize( aster_comm_t *comm = aster_get_current_comm() );

/** @brief Get MPI rank */
int getMPIRank( aster_comm_t *comm = aster_get_current_comm() );

#ifdef ASTER_HAVE_MPI

/* this code is inspired by the MPI class of the dolfin project :
 * fenicsproject.org */

class AsterMPI {
  private:
    template < typename T >
    struct dependent_false : std::false_type {};
    template < typename T >
    static MPI_Datatype mpi_type() {
        static_assert( dependent_false< T >::value, "Unknown MPI type" );
        throw std::runtime_error( "Unknown MPI type" );
        return MPI_CHAR;
    }

  public:
    /// Gather values from each process (variable count per process)
    template < typename T >
    static void all_gather( const std::vector< T > &in_values,
                            std::vector< std::vector< T > > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values from each process (variable count per process)
    template < typename T >
    static void all_gather( const std::vector< T > &in_values, std::vector< T > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values from each process (variable count per process)
    template < typename T >
    static void all_gather( const std::set< T > &in_values, std::set< T > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values from each process (variable count per process)
    template < typename T1, typename T2 >
    static void all_gather( const std::map< T1, T2 > &in_values, std::map< T1, T2 > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values from each process (variable count per process)
    template < typename T >
    static void all_gather( const std::vector< T > &in_values, JeveuxVector< T > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values, one primitive from each process (MPI_Allgather)
    template < typename T >
    static void all_gather( const T in_value, std::vector< T > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values, one primitive from each process (MPI_Allgather).
    /// Specialization for std::string
    static void all_gather( const std::string &in_values, VectorString &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values, one primitive from each process (MPI_Allgather).
    /// Specialization for VectorString
    static void all_gather( const VectorString &in_values, VectorString &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Gather values, one primitive from each process (MPI_Allgather).
    /// Specialization for JeveuxString
    template < int length >
    static void all_gather( const std::vector< JeveuxString< length > > &in_values,
                            std::vector< JeveuxString< length > > &out_values,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// AllReduce values, one value from each process (MPI_AllReduce).
    template < typename T >
    static void all_reduce( const T in_value, T &out_value, MPI_Op op,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );
    template < typename T >
    static T all_reduce( const T in_value, MPI_Op op,
                         aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// AllReduce a vector (MPI_AllReduce).
    template < typename T >
    static void all_reduce( const std::vector< T > in_value, std::vector< T > &out_value, MPI_Op op,
                            aster_comm_t *_commCurrent = aster_get_current_comm() );
    template < typename T >
    static void all_reduce( const JeveuxVector< T > in_value, JeveuxVector< T > &out_value,
                            MPI_Op op, aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// AllReduce values for given OP, one value from each process
    /// (MPI_AllReduce).
    template < typename T >
    static T max( const T in_value, aster_comm_t *_commCurrent = aster_get_current_comm() );
    template < typename T >
    static T min( const T in_value, aster_comm_t *_commCurrent = aster_get_current_comm() );
    template < typename T >
    static T sum( const T in_value, aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Broadcast one value from root
    template < typename T >
    static void bcast( T &value, int root, aster_comm_t *_commCurrent = aster_get_current_comm() );

    static void bcast( std::string &value, int root,
                       aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Broadcast a vector from root
    template < typename T >
    static void bcast( std::vector< T > &value, int root,
                       aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Barrier
    static void barrier( aster_comm_t *_commCurrent = aster_get_current_comm() ) {
        aster_set_mpi_barrier( _commCurrent );
    };

    template < typename T >
    static void send( const std::vector< T > &in_values, int dest, int tag,
                      aster_comm_t *_commCurrent = aster_get_current_comm() );

    template < typename T >
    static void receive( const std::vector< T > &out_values, int source, int tag,
                         aster_comm_t *_commCurrent = aster_get_current_comm() );

    template < typename T >
    static void send_receive( const std::vector< T > &in_values, std::vector< T > &out_values,
                              int peer, int tag,
                              aster_comm_t *_commCurrent = aster_get_current_comm() );

    /// Get Current Comm
    aster_comm_t *getCurrentCommunicator() { return aster_get_current_comm(); }

    /// Set Current Comm
    void setCurrentCommunicator( aster_comm_t *commCurrent ) {
        aster_set_current_comm( commCurrent );
    };

    void freeCurrentCommunicator() {
        auto commCurr = getCurrentCommunicator();
        auto parent = commCurr->parent;
        setCurrentCommunicator( parent );
        aster_free_comm( commCurr );
    };

    void freeCommunicator( aster_comm_t *commCurr ) { aster_free_comm( commCurr ); }

    // Split communicator
    aster_comm_t *splitCommunicator( int color,
                                     aster_comm_t *_commCurrent = aster_get_current_comm() );
};

//---------------------------------------------------------------------------
template <>
inline MPI_Datatype AsterMPI::mpi_type< char >() {
    return MPI_CHAR;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< float >() {
    return MPI_FLOAT;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< double >() {
    return MPI_DOUBLE;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< short int >() {
    return MPI_SHORT;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< int >() {
    return MPI_INT;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< long int >() {
    return MPI_LONG;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< unsigned int >() {
    return MPI_UNSIGNED;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< unsigned long int >() {
    return MPI_UNSIGNED_LONG;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< long long >() {
    return MPI_LONG_LONG;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< bool >() {
    return MPI_CXX_BOOL;
}
template <>
inline MPI_Datatype AsterMPI::mpi_type< ASTERCOMPLEX >() {
    return MPI_C_DOUBLE_COMPLEX;
}
//---------------------------------------------------------------------------
inline void AsterMPI::all_gather( const std::string &in_values, VectorString &out_values,
                                  aster_comm_t *_commCurrent ) {
    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts;
    AsterMPI::all_gather( int( in_values.size() ), pcounts, _commCurrent );

    // Build offsets
    VectorInt offsets( comm_size + 1, 0 );
    for ( std::size_t i = 1; i <= comm_size; ++i )
        offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather
    const std::size_t n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );
    std::vector< char > _out( n );
    aster_mpi_allgatherv( const_cast< char * >( in_values.data() ), in_values.size(), MPI_CHAR,
                          _out.data(), pcounts.data(), offsets.data(), MPI_CHAR, _commCurrent );

    // Rebuild
    out_values.resize( comm_size );
    for ( std::size_t p = 0; p < comm_size; ++p ) {
        out_values[p] = std::string( _out.begin() + offsets[p], _out.begin() + offsets[p + 1] );
    }
}
//---------------------------------------------------------------------------
inline void AsterMPI::all_gather( const VectorString &in_values, VectorString &out_values,
                                  aster_comm_t *_commCurrent ) {
    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts_lc, pcounts_gl;
    pcounts_lc.reserve( in_values.size() );
    for ( auto str : in_values ) {
        pcounts_lc.push_back( str.size() );
    }
    AsterMPI::all_gather( pcounts_lc, pcounts_gl, _commCurrent );

    const std::size_t local_size = std::accumulate( pcounts_lc.begin(), pcounts_lc.end(), 0 );
    const std::size_t global_size = std::accumulate( pcounts_gl.begin(), pcounts_gl.end(), 0 );

    // Gather
    std::vector< char > _in, _out;
    _in.reserve( local_size );
    for ( auto &str : in_values ) {
        for ( auto &carac : str )
            _in.push_back( carac );
    }
    AS_ASSERT( _in.size() == local_size );

    AsterMPI::all_gather( _in, _out, _commCurrent );

    // Rebuild
    VectorInt offsets_gl( pcounts_gl.size() + 1, 0 );
    for ( std::size_t i = 1; i <= pcounts_gl.size(); ++i )
        offsets_gl[i] = offsets_gl[i - 1] + pcounts_gl[i - 1];

    out_values.clear();
    out_values.resize( pcounts_gl.size() );
    for ( std::size_t p = 0; p < pcounts_gl.size(); ++p ) {
        out_values[p] =
            std::string( _out.begin() + offsets_gl[p], _out.begin() + offsets_gl[p + 1] );
    }
}
//---------------------------------------------------------------------------
template < int length >
inline void AsterMPI::all_gather( const std::vector< JeveuxString< length > > &in_values,
                                  std::vector< JeveuxString< length > > &out_values,
                                  aster_comm_t *_commCurrent ) {
    // Get number of procs
    const std::size_t comm_size = getMPISize();
    const std::size_t comm_rank = getMPIRank();

    typedef JeveuxString< length > JeveuxChar;
    VectorInt sizes2;
    AsterMPI::all_gather( int( in_values.size() ), sizes2, _commCurrent );

    out_values.clear();
    for ( int rank = 0; rank < comm_size; ++rank ) {
        JeveuxChar *retour = new JeveuxChar[sizes2[rank]];

        if ( rank == comm_rank )
            for ( int position = 0; position < sizes2[rank]; ++position )
                retour[position] = in_values[position];

        aster_mpi_bcast( retour, length * sizes2[rank], MPI_CHAR, rank, _commCurrent );
        for ( int position = 0; position < sizes2[rank]; ++position )
            out_values.push_back( JeveuxChar( retour[position] ) );
        delete[] retour;
    }
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const std::vector< T > &in_values,
                           std::vector< std::vector< T > > &out_values,
                           aster_comm_t *_commCurrent ) {

    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    AsterMPI::all_gather( local_size, pcounts, _commCurrent );
    assert( pcounts.size() == comm_size );

    // Build offsets
    VectorInt offsets( comm_size + 1, 0 );
    for ( std::size_t i = 1; i <= comm_size; ++i )
        offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );
    std::vector< T > recvbuf( n );
    aster_mpi_allgatherv( const_cast< T * >( in_values.data() ), in_values.size(), mpi_type< T >(),
                          recvbuf.data(), pcounts.data(), offsets.data(), mpi_type< T >(),
                          _commCurrent );

    // Repack data
    out_values.resize( comm_size );
    for ( std::size_t p = 0; p < comm_size; ++p ) {
        out_values[p].resize( pcounts[p] );
        for ( int i = 0; i < pcounts[p]; ++i )
            out_values[p][i] = recvbuf[offsets[p] + i];
    }
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const std::vector< T > &in_values, std::vector< T > &out_values,
                           aster_comm_t *_commCurrent ) {

    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    AsterMPI::all_gather( local_size, pcounts, _commCurrent );
    assert( pcounts.size() == comm_size );

    // Build offsets
    VectorInt offsets( comm_size + 1, 0 );
    for ( std::size_t i = 1; i <= comm_size; ++i )
        offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );
    out_values.resize( n );
    aster_mpi_allgatherv( const_cast< T * >( in_values.data() ), in_values.size(), mpi_type< T >(),
                          out_values.data(), pcounts.data(), offsets.data(), mpi_type< T >(),
                          _commCurrent );
}

//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const std::set< T > &in_values, std::set< T > &out_values,
                           aster_comm_t *_commCurrent ) {
    std::vector< T > vect_in( in_values.begin(), in_values.end() ), vect_out;

    AsterMPI::all_gather( vect_in, vect_out, _commCurrent );

    std::copy( vect_out.begin(), vect_out.end(), std::inserter( out_values, out_values.end() ) );
}

//---------------------------------------------------------------------------
template < typename T1, typename T2 >
void AsterMPI::all_gather( const std::map< T1, T2 > &in_values, std::map< T1, T2 > &out_values,
                           aster_comm_t *_commCurrent ) {

    std::vector< T1 > key_in, key_out;
    std::vector< T2 > val_in, val_out;

    key_in.reserve( in_values.size() );
    val_in.reserve( in_values.size() );

    for ( const auto &[key, value] : in_values ) {
        key_in.push_back( key );
        val_in.push_back( value );
    }

    AsterMPI::all_gather( key_in, key_out, _commCurrent );
    AsterMPI::all_gather( val_in, val_out, _commCurrent );

    out_values.clear();
    auto nb_val = key_out.size();
    for ( auto i_val = 0; i_val < nb_val; i_val++ ) {
        out_values[key_out[i_val]] = val_out[i_val];
    }
}

//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const std::vector< T > &in_values, JeveuxVector< T > &out_values,
                           aster_comm_t *_commCurrent ) {

    // Get number of procs
    const std::size_t comm_size = getMPISize();

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    AsterMPI::all_gather( local_size, pcounts, _commCurrent );
    assert( pcounts.size() == comm_size );

    // Build offsets
    VectorInt offsets( comm_size + 1, 0 );
    for ( std::size_t i = 1; i <= comm_size; ++i )
        offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );

    if ( out_values->size() < n ) {
        if ( out_values.exists() )
            out_values->deallocate();
        out_values->allocate( n );
    }
    out_values->updateValuePointer();

    aster_mpi_allgatherv( const_cast< T * >( in_values.data() ), in_values.size(), mpi_type< T >(),
                          out_values->getDataPtr(), pcounts.data(), offsets.data(), mpi_type< T >(),
                          _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_gather( const T in_value, std::vector< T > &out_values,
                           aster_comm_t *_commCurrent ) {
    out_values.resize( getMPISize() );
    aster_mpi_allgather( const_cast< T * >( &in_value ), 1, mpi_type< T >(), out_values.data(), 1,
                         mpi_type< T >(), _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_reduce( const T in_value, T &out_value, MPI_Op op, aster_comm_t *_commCurrent ) {
    aster_mpi_allreduce( const_cast< T * >( &in_value ), static_cast< T * >( &out_value ), 1,
                         mpi_type< T >(), op, _commCurrent );
}
template < typename T >
T AsterMPI::all_reduce( const T in_value, MPI_Op op, aster_comm_t *_commCurrent ) {
    T out_value;
    AsterMPI::all_reduce( in_value, out_value, op, _commCurrent );

    return out_value;
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_reduce( const std::vector< T > in_value, std::vector< T > &out_value, MPI_Op op,
                           aster_comm_t *_commCurrent ) {
    if ( in_value.size() != out_value.size() ) {
        throw std::runtime_error( "Inconsistent sizes" );
    }
    aster_mpi_allreduce( const_cast< T * >( in_value.data() ),
                         static_cast< T * >( out_value.data() ), in_value.size(), mpi_type< T >(),
                         op, _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::all_reduce( const JeveuxVector< T > in_value, JeveuxVector< T > &out_value,
                           MPI_Op op, aster_comm_t *_commCurrent ) {
    if ( in_value->size() != out_value->size() ) {
        throw std::runtime_error( "Inconsistent sizes" );
    }
    in_value->updateValuePointer();
    out_value->updateValuePointer();
    aster_mpi_allreduce( const_cast< T * >( in_value->getDataPtr() ),
                         static_cast< T * >( out_value->getDataPtr() ), in_value->size(),
                         mpi_type< T >(), op, _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
T AsterMPI::max( const T in_value, aster_comm_t *_commCurrent ) {
    return AsterMPI::all_reduce( in_value, MPI_MAX, _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
T AsterMPI::min( const T in_value, aster_comm_t *_commCurrent ) {
    return AsterMPI::all_reduce( in_value, MPI_MIN, _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
T AsterMPI::sum( const T in_value, aster_comm_t *_commCurrent ) {
    return AsterMPI::all_reduce( in_value, MPI_SUM, _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::bcast( T &value, int root, aster_comm_t *_commCurrent ) {
    aster_mpi_bcast( const_cast< T * >( &value ), 1, mpi_type< T >(), root, _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::bcast( std::vector< T > &value, int root, aster_comm_t *_commCurrent ) {
    aster_mpi_bcast( const_cast< T * >( &value[0] ), value.size(), mpi_type< T >(), root,
                     _commCurrent );
}
//---------------------------------------------------------------------------
template < typename T >
void AsterMPI::send( const std::vector< T > &in_values, int dest, int tag,
                     aster_comm_t *_commCurrent ) {
    const T &firstValue = in_values[0];
    T *firstValueP = const_cast< T * >( &firstValue );
    aster_mpi_send( firstValueP, in_values.size(), mpi_type< T >(), dest, tag, _commCurrent );
}

template < typename T >
void AsterMPI::receive( const std::vector< T > &out_values, int source, int tag,
                        aster_comm_t *_commCurrent ) {
    const T &firstValue = out_values[0];
    T *firstValueP = const_cast< T * >( &firstValue );
    aster_mpi_recv( firstValueP, out_values.size(), mpi_type< T >(), source, tag, _commCurrent );
}

template < typename T >
void AsterMPI::send_receive( const std::vector< T > &in_values, std::vector< T > &out_values,
                             int peer, int tag, aster_comm_t *_commCurrent ) {

    int size_in = in_values.size(), size_out = 0;
    aster_mpi_sendrecv( const_cast< int * >( &size_in ), 1, mpi_type< int >(), peer, tag,
                        static_cast< int * >( &size_out ), 1, mpi_type< int >(), peer, tag,
                        _commCurrent );

    out_values.resize( size_out );

    aster_mpi_sendrecv( const_cast< T * >( in_values.data() ), in_values.size(), mpi_type< T >(),
                        peer, tag, out_values.data(), size_out, mpi_type< T >(), peer, tag,
                        _commCurrent );
};

#endif /* ASTER_HAVE_MPI */

#endif /* ASTERMPI_H_ */
