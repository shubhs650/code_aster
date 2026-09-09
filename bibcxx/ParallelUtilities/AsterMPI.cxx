/**
 * @file AsterMPI.cxx
 * @brief Implementation de AsterMPI
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

#include "ParallelUtilities/AsterMPI.h"

int getMPISize( aster_comm_t *comm ) {
#ifdef ASTER_HAVE_MPI
    int rank = -1, nbProcs = -1;
    aster_get_mpi_info( comm, &rank, &nbProcs );
    if ( rank == -1 || nbProcs == -1 )
        throw std::runtime_error( "Error with MPI Infos" );
#else
    int nbProcs = 1;
#endif
    return nbProcs;
};

int getMPIRank( aster_comm_t *comm ) {
#ifdef ASTER_HAVE_MPI
    int rank = -1, nbProcs = -1;
    aster_get_mpi_info( comm, &rank, &nbProcs );
    if ( rank == -1 || nbProcs == -1 )
        throw std::runtime_error( "Error with MPI Infos" );
#else
    int rank = 0;
#endif
    return rank;
};

#ifdef ASTER_HAVE_MPI

aster_comm_t *AsterMPI::splitCommunicator( int color, aster_comm_t *_commCurrent ) {
    return aster_split_comm( _commCurrent, color, getMPIRank(),
                             std::string( "COMM_SPLIT" ).data() );
};

void AsterMPI::bcast( std::string &value, int root, aster_comm_t *_commCurrent ) {
    aster_mpi_bcast( value.data(), value.size(), MPI_CHAR, root, _commCurrent );
}

#endif
