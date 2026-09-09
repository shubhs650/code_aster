#ifndef CONNECTIONMESH_H_
#define CONNECTIONMESH_H_

/**
 * @file ConnectionMesh.h
 * @brief Fichier entete de la classe
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

#include "Meshes/BaseMesh.h"
#include "Meshes/ParallelMesh.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ConnectionMesh
 * @brief Cette classe decrit un maillage partiel reconstruit a partir d'une liste de groupe de
 * noeuds
 * @author Nicolas Sellenet
 */
class ConnectionMesh : public BaseMesh {
  private:
    typedef JeveuxCollection< ASTERINTEGER, NamesMapChar24 > JeveuxCollectionLongNamePtr;
    /** @brief Base ParallelMesh */
    ParallelMeshPtr _pMesh;
    /** @brief id of node in local numbering */
    JeveuxVectorLong _nodesLocalNumbering;
    /** @brief id of node in global numbering */
    JeveuxVectorLong _nodesGlobalNumbering;
    /** @brief number of owner proc for each nodes */
    JeveuxVectorLong _nodesOwner;
    /** @brief id of cell in local numbering */
    JeveuxVectorLong _cellsLocalNumbering;
    /** @brief number of owner proc for each cells */
    JeveuxVectorLong _cellsOwner;

    VectorLong getCellsGlobalNumbering( const JeveuxVectorLong &rankOfCells ) const;

  public:
    /**
     * @typedef ConnectionMeshPtr
     * @brief Pointeur intelligent vers un ConnectionMesh
     */
    typedef std::shared_ptr< ConnectionMesh > ConnectionMeshPtr;

    /**
     * @brief Constructeur
     */
    ConnectionMesh( const ParallelMeshPtr &mesh, const VectorString &groupsOfNodes,
                    const VectorString &groupsOfCells )
        : ConnectionMesh( ResultNaming::getNewResultName(), mesh, groupsOfNodes, groupsOfCells ) {};

    /**
     * @brief Constructeur
     */
    ConnectionMesh( const std::string &name, const ParallelMeshPtr &mesh,
                    const VectorString &groupsOfNodes, const VectorString &groupsOfCells );

    const JeveuxVectorLong &getNodesGlobalNumbering() const { return _nodesGlobalNumbering; };

    const JeveuxVectorLong &getNodesLocalNumbering() const { return _nodesLocalNumbering; };

    const JeveuxVectorLong getNodesOwner() const { return _nodesOwner; };

    const JeveuxVectorLong &getCellsLocalNumbering() const { return _cellsLocalNumbering; };

    const JeveuxVectorLong getCellsOwner() const { return _cellsOwner; };

    const ParallelMeshPtr &getParallelMesh() const { return _pMesh; };

    VectorString getGroupsOfCells( const bool local = false ) const;

    VectorString getGroupsOfNodes( const bool local = false ) const;

    bool hasGroupOfCells( const std::string &name, const bool local = false ) const;

    bool hasGroupOfNodes( const std::string &name, const bool local = false ) const;

    VectorLong getCells( const std::string name = "" ) const;

    VectorLong getCells( const VectorString &names = {} ) const;

    /**
     * @brief Returns the nodes indexes of a group of cells
     * @param name name of group of cells
     * @param local node id in local or global numbering
     * @param same_rank keep or not the nodes owned by the current domain
     * @return list of nodes indexes
     */
    VectorLong getNodesFromCells( const std::string, const bool localNumbering = false,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    VectorLong getNodesFromCells( const VectorString &, const bool localNumbering = false,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    VectorLong getNodesFromCells( const VectorLong &cells, const bool localNumbering = false,
                                  const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Fonction permettant de savoir si un maillage est partiel
     * @return retourne true si le maillage est partiel
     */
    virtual bool isConnection() const { return true; };

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
 * @typedef ConnectionMeshPtr
 * @brief Pointeur intelligent vers un ConnectionMesh
 */
typedef std::shared_ptr< ConnectionMesh > ConnectionMeshPtr;

#endif /* ASTER_HAVE_MPI */

#endif /* CONNECTIONMESH_H_ */
