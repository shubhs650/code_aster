#ifndef MESH_H_
#define MESH_H_

/**
 * @file Mesh.h
 * @brief Fichier entete de la classe Mesh
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

#include "Meshes/BaseMesh.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/ResultNaming.h"

#include <filesystem>

/**
 * @class Mesh
 * @brief Cette classe decrit un maillage Aster
 */
class Mesh : public BaseMesh {
  protected:
    /** @brief Pointeur de nom Jeveux '.NOMNOE' */
    NamesMapChar8 _nameOfNodes;
    /** @brief Pointeur de nom Jeveux '.NOMMAIL' */
    NamesMapChar8 _nameOfCells;

    /**
     * @brief Constructeur
     */
    Mesh( const std::string name, const std::string type )
        : BaseMesh( name, type ),
          _nameOfNodes( getName() + ".NOMNOE" ),
          _nameOfCells( getName() + ".NOMMAI" ) {};

  public:
    /**
     * @typedef MeshPtr
     * @brief Pointeur intelligent vers un Mesh
     */
    typedef std::shared_ptr< Mesh > MeshPtr;

    /**
     * @brief Constructeur
     */
    Mesh()
        : BaseMesh( ResultNaming::getNewResultName(), "MAILLAGE" ),
          _nameOfNodes( getName() + ".NOMNOE" ),
          _nameOfCells( getName() + ".NOMMAI" ) {};

    /**
     * @brief Constructeur
     */
    Mesh( const std::string name )
        : BaseMesh( name, "MAILLAGE" ),
          _nameOfNodes( getName() + ".NOMNOE" ),
          _nameOfCells( getName() + ".NOMMAI" ) {};

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

    VectorLong getCells( const std::string name ) const;

    VectorLong getCells( const VectorString &names = {} ) const;

    /**
     * @brief Return list of nodes
     * @param name name of group (if empty all the nodes)
     * @param local node id in local or global numbering
     * @param same_rank keep or not the nodes owned by the current domain
     * @return list of Nodes
     */
    VectorLong getNodes( const std::string name = std::string(), const bool localNumbering = false,
                         const ASTERINTEGER same_rank = PythonBool::None ) const;

    VectorLong getNodes( const VectorString &names, const bool localNumbering = true,
                         const ASTERINTEGER same_rank = PythonBool::None ) const;
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
     * @brief Get inner nodes
     * @return list of node ids
     */
    VectorLong getInnerNodes() const { return this->getNodes( std::string(), true ); };

    /**
     * @brief Get inner nodes
     * @return list of node ids
     */
    const JeveuxVectorLong getNodesOwner() const {
        return JeveuxVectorLong( getName() + ".NOEX",
                                 VectorLong( getNumberOfNodes(), getMPIRank() ) );
    };

    /**
     * @brief Get inner cells
     * @return list of cells ids
     */
    const JeveuxVectorLong getCellsOwner() const {
        return JeveuxVectorLong( getName() + ".MAEX",
                                 VectorLong( getNumberOfCells(), getMPIRank() ) );
    };

    bool isQuadratic( const bool local = false ) const;

    /** @brief Test if all cells in group are skin cells */
    bool isSkin( const std::string groupName ) const;

    /**
     * @brief Read a Aster Mesh file
     * @return retourne true si tout est ok
     */
    bool readAsterFile( const std::filesystem::path &fileName );

    /**
     * @brief Read a Gibi Mesh file
     * @return retourne true si tout est ok
     */
    bool readGibiFile( const std::filesystem::path &fileName );

    /**
     * @brief Read a Gmsh Mesh file
     * @return retourne true si tout est ok
     */
    bool readGmshFile( const std::filesystem::path &fileName );

    MeshPtr convertToLinear( const ASTERINTEGER info = 1 );

    MeshPtr convertToQuadratic( const ASTERINTEGER info = 1 );

    MeshPtr convertToBiQuadratic( const ASTERINTEGER info = 1 );

    MeshPtr convertToCubic( const ASTERINTEGER info = 1 );

    MeshPtr fix( const bool remove_orphan, const bool positive_measure, const bool outward_normal,
                 const bool double_nodes, const bool double_cells, const ASTERDOUBLE tole,
                 const ASTERINTEGER info = 1 );

    MeshPtr getOctreeMesh( const ASTERINTEGER nb_max_pt, const ASTERINTEGER nb_max_level );

    void addNodeLabels( const VectorString &labels );

    void addCellLabels( const VectorString &labels );

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
 * @typedef MeshPtr
 * @brief Pointeur intelligent vers un Mesh
 */
typedef std::shared_ptr< Mesh > MeshPtr;

#endif /* MESH_H_ */
