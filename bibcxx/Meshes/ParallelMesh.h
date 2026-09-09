#ifndef PARALLELMESH_H_
#define PARALLELMESH_H_

/**
 * @file ParallelMesh.h
 * @brief Fichier entete de la classe ParallelMesh
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

#include <set>

#ifdef ASTER_HAVE_MPI

#include "MemoryManager/NamesMap.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/Joints.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ParallelMesh
 * @brief Cette classe decrit un maillage Aster parallèle
 * @author Nicolas Sellenet
 */
class ParallelMesh : public BaseMesh {
  private:
    typedef std::set< std::string > SetOfString;
    typedef SetOfString::iterator SetOfStringIter;
    typedef SetOfString::const_iterator SetOfStringCIter;

  protected:
    /** @brief All groups of nodes (parallel mesh) */
    NamesMapChar24 _globalGroupOfNodes;
    /** @brief Set of all groups of nodes (parallel mesh) */
    SetOfString _setOfAllGON;
    /** @brief All groups of cells (parallel mesh) */
    NamesMapChar24 _globalGroupOfCells;
    /** @brief Set of all groups of cells (parallel mesh) */
    SetOfString _setOfAllGOE;
    /** @brief Identify outer nodes */
    JeveuxVectorLong _nodesOwner;
    /** @brief Identify outer cells */
    JeveuxVectorLong _cellsOwners;
    /** @brief Global node numbering */
    JeveuxVectorLong _globalNodeIds;
    /** @brief Global cell numbering */
    JeveuxVectorLong _globalCellIds;
    /** @brief Objet Jeveux '.LASTGHOLAYER' */
    JeveuxVectorShort _lastGhostsLayer;
    /** @brief List of joints */
    JointsPtr _joints;
    /** @brief Global to local node number */
    std::shared_ptr< MapLong > _global2localNodeIdsPtr;

    void _buildGlobal2LocalNodeIdsMapPtr();

  public:
    /**
     * @typedef ParallelMeshPtr
     * @brief Pointeur intelligent vers un ParallelMesh
     */
    typedef std::shared_ptr< ParallelMesh > ParallelMeshPtr;

    /**
     * @brief Constructeur
     */
    ParallelMesh() : ParallelMesh( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     */
    ParallelMesh( const std::string &name );

    /**
     * @brief Return the MED type for each cell
     */
    const VectorLong getAllMedCellsTypes() const;

    /**
     * @brief Get the JeveuxVector for outer subdomain nodes
     * @return _nodesOwner
     */
    const JeveuxVectorLong getNodesOwner() const { return _nodesOwner; };

    /**
     * @brief Get the JeveuxVector for outer subdomain cells
     * @return _cellsOwners
     */
    const JeveuxVectorLong getCellsOwner() const { return _cellsOwners; };

    bool isQuadratic( const bool local = false ) const;

    bool hasGroupOfCells( const std::string &name, const bool local = false ) const;

    bool hasGroupOfNodes( const std::string &name, const bool local = false ) const;

    VectorString getGroupsOfCells( const bool local = false ) const;

    VectorString getGroupsOfNodes( const bool local = false ) const;

    /**
     * @brief Create a group of cells
     */
    void setGroupOfCells( const std::string &name, const VectorLong &cell_ids );

    /**
     * @brief Create a group of nodes
     */
    void setGroupOfNodes( const std::string &name, const VectorLong &node_ids,
                          const bool localNumbering = false );

    /**
     * @brief Steal from input vector ids in local numbering of ghost nodes on the last layer
     */
    void setLastGhostsLayer( VectorInt &node_ids ) {
        // copy data to JeveuxVector _lastGhostsLayer
        *_lastGhostsLayer = node_ids;
        // reset input vector to empty state
        node_ids = VectorInt();
    };

    VectorLong getCells( const std::string name ) const;

    VectorLong getCells( const VectorString &names = {} ) const;

    /**
     * @brief Return ids in local numbering of ghost nodes of the last layer
     */
    JeveuxVectorShort getLastGhostsLayer() const {
        _lastGhostsLayer->updateValuePointer();
        return _lastGhostsLayer;
    };

    /**
     * @brief Return list of nodes
     * @param name name of group (if empty all the nodes)
     * @param localNumbering node id in local or global numbering
     * @param same_rank keep or not the nodes owned by the current domain
     * @return list of Nodes
     */

    VectorLong getNodes( const VectorString &names, const bool localNumbering = false,
                         const ASTERINTEGER same_rank = PythonBool::None ) const;

    VectorLong getNodes( const std::string name = std::string(), const bool localNumbering = false,
                         const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Get inner nodes
     * @return list of node ids
     */
    VectorLong getInnerNodes() const {
        return this->getNodes( std::string(), true, PythonBool::True );
    };

    /**
     * @brief Get outer nodes
     * @return list of node ids
     */
    VectorLong getOuterNodes() const {
        return this->getNodes( std::string(), true, PythonBool::False );
    };

    /**
     * @brief Get inner nodes
     * @return list of node ids
     */
    VectorLong getInnerCells() const;

    /**
     * @brief Get outer nodes
     * @return list of node ids
     */
    VectorLong getOuterCells() const;

    /**
     * @brief Get the mapping between local and global numbering of nodes
     * @return JeveuxVector of the indirection
     */
    const JeveuxVectorLong getLocalToGlobalNodeIds() const;

    std::shared_ptr< const MapLong > getGlobalToLocalNodeIds() const {
        return _global2localNodeIdsPtr;
    };

    /**
     * @brief Get the mapping between local and global numbering of cells
     * @return JeveuxVector of the indirection
     */
    const JeveuxVectorLong getLocalToGlobalCellIds() const;

    /**
     * @brief Returns the nodes indexes of a group of cells
     * @param name name of group of cells
     * @param local node id in local or global numbering
     * @param same_rank keep or not the nodes owned by the current domain
     * @return list of Nodes
     */
    VectorLong getNodesFromCells( const std::string name, const bool localNumbering = false,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    VectorLong getNodesFromCells( const VectorString &names, const bool localNumbering = false,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    VectorLong getNodesFromCells( const VectorLong &cells, const bool localNumbering = false,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    /** @brief Returns the number of joints */
    const JeveuxVectorLong &getOppositeDomains() const { return _joints->getOppositeDomains(); };

    /** @brief Returns a joint */
    VectorLong getSendJoint( const int &id ) const;
    VectorLong getReceiveJoint( const int &id ) const;

    /**
     * @brief Fonction permettant de savoir si un maillage est parallel
     * @return retourne true si le maillage est parallel
     */
    virtual bool isParallel() const { return true; };

    /**
     * @brief Get the local number of a node from its local number
     * @return Return the local number if the node if present on the subdomain ;
     * otherwise raise an exception
     */
    ASTERINTEGER getGlobalToLocalNodeId( const ASTERINTEGER &nodeId,
                                         const bool &stop = true ) const;

    bool updateGlobalGroupOfNodes( void );

    bool updateGlobalGroupOfCells( void );

    bool build();

    ParallelMeshPtr convertToLinear( const ASTERINTEGER info = 1 );

    ParallelMeshPtr convertToQuadratic( const ASTERINTEGER info = 1 );

    ParallelMeshPtr convertToBiQuadratic( const ASTERINTEGER info = 1 );

    ParallelMeshPtr convertToCubic( const ASTERINTEGER info = 1 );

    ParallelMeshPtr fix( const bool remove_orphan, const bool positive_measure,
                         const bool outward_normal, const bool double_nodes,
                         const bool double_cells, const ASTERDOUBLE tole,
                         const ASTERINTEGER info = 1 );

    /* Mesh builder functions */
    void create_joints( const VectorLong &domains, const VectorLong &globalNodeIds,
                        const VectorLong &nodesOwner, const VectorLong &globalCellIds,
                        const VectorOfVectorsLong &joints, const ASTERINTEGER &nbLayer );

    void endDefinition();

    VectorOfVectorsLong getNodesRanks() const;

    VectorOfVectorsLong getCellsRanks() const;

    /**
     * @brief Impression du maillage au format MED
     * @param fileName Nom du fichier MED à imprimer
     * @param local booleen pour demander une impression locale
     * @param version version de med de l'impression
     * @param mode mode d'ecriture 'a' ou 'w'
     * @return true
     */
    bool printMedFile( const std::filesystem::path &fileName, bool local = true,
                       std::array< int, 3 > version = { 0, 0, 0 },
                       const std::string &mode = "a" ) const;
};

/**
 * @typedef ParallelMeshPtr
 * @brief Pointeur intelligent vers un ParallelMesh
 */
typedef std::shared_ptr< ParallelMesh > ParallelMeshPtr;

#endif /* ASTER_HAVE_MPI */

#endif /* PARALLELMESH_H_ */
