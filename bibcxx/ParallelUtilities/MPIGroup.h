#ifndef MPIGROUP_H_
#define MPIGROUP_H_

/**
 * @file MPIGroup.h
 * @brief Header of connection graph
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

#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#include "aster_mpi.h"

#include "MemoryManager/JeveuxVector.h"

/**
 * @class MPIGroup
 * @brief Class describing a mpi group with jeveux vectors
 * @author Nicolas Sellenet
 */
class MPIGroup {
    /** @brief MPI parent group */
    MPI_Group _parentGroup;
    /** @brief MPI parent communicator */
    MPI_Comm _parentComm;
    /** @brief MPI current group */
    MPI_Group _currentGroup;
    /** @brief MPI current communicator */
    MPI_Comm _groupComm;
    /** @brief booleans to manage groups */
    bool _groupCreated = false;
    bool _newGroupCreated = false;
    aster_comm_t _asterComm;
    /** @brief booleans to manage groups */
    int _groupSize = 0;

  public:
    /**
     * @brief Constructor
     */
    MPIGroup() {};

    ~MPIGroup();

    /** @brief Build group from parent communicator and proc id list */
    void buildFromProcsVector( const MPI_Comm &parentComm, const VectorInt &procIdVector );

    /** @brief Get aster_comm_t communicator */
    const aster_comm_t *getAsterCommunicator() const;

    /** @brief Get communicator */
    MPI_Comm getCommunicator() const;

    /** @brief Get communicator */
    int getGroupSize() const { return _groupSize; };
};

using MPIGroupPtr = std::shared_ptr< MPIGroup >;

#endif /* ASTER_HAVE_MPI */

#endif /* MPIGROUP_H_ */
