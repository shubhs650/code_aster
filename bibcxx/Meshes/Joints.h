/**
 * @file Joints.h
 * @brief Fichier entete de la classe Joints
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

#pragma once

#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "ParallelUtilities/MPIGroup.h"

/**
 * @class Joints
 * @brief Cette classe decrit un maillage Aster parallèle
 * @author Nicolas Sellenet
 */
class Joints : public DataStructure {
  private:
    /** @brief Number of ghost layers */
    JeveuxVectorLong _layer;
    /** @brief List of connected domains */
    JeveuxVectorLong _domj;
    /** @brief two-side index of nodes to send */
    JeveuxCollectionLong _send;
    /** @brief two-side index of nodes to receive */
    JeveuxCollectionLong _recv;
    /** @brief vector of procs involved in MPI group (id given in parent group) */
    JeveuxVectorShort _procGroupIds;
    /** @brief MPI group communicator */
    JeveuxVectorLong _groupComm;
    /** @brief MPI group for joints */
    MPIGroupPtr _jointsMPIGroup;
    /** @brief Joints hash container */
    JeveuxVectorChar24 _hashContainer;

    void buildGroup();

  public:
    /**
     * @typedef JointsPtr
     * @brief Pointeur intelligent vers un Joints
     */
    typedef std::shared_ptr< Joints > JointsPtr;

    /**
     * @brief Constructeur
     */
    Joints();

    Joints( const std::string name );

    void setOppositeDomains( const VectorLong &oppositeDomains );

    const JeveuxVectorLong &getOppositeDomains() const;

    VectorLong getSendedElements( const ASTERINTEGER &id ) const;

    VectorLong getReceivedElements( const ASTERINTEGER &id ) const;

    void setSendedElements( const VectorOfVectorsLong &send );

    void setReceivedElements( const VectorOfVectorsLong &recv );

    void setNumberOfGhostLayer( const ASTERINTEGER &nb_layer );

    ASTERINTEGER getNumberOfGhostLayer( void ) const;

    bool build();
};

/**
 * @typedef JointsPtr
 * @brief Pointeur intelligent vers un Joints
 */
typedef std::shared_ptr< Joints > JointsPtr;

#endif /* ASTER_HAVE_MPI */
