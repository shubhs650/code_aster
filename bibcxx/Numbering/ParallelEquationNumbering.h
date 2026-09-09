/**
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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/Joints.h"
#include "Numbering/EquationNumbering.h"

#ifdef ASTER_HAVE_MPI

/**
 * @class ParallelEquationNumbering
 * @brief Class definissant un NUME_EQUA
 */
class ParallelEquationNumbering : public EquationNumbering {
  protected:
    /** @brief Objet Jeveux '.NULG' */
    JeveuxVectorLong _localToGlobal;
    /** @brief Objet Jeveux '.PDDL' */
    JeveuxVectorLong _localToRank;
    /** @brief Objet Jeveux '.LILT' */
    JeveuxVectorLong _globalDofNumIds;
    /** @brief List of joints */
    JointsPtr _joints;

    std::unordered_map< ASTERINTEGER, ASTERINTEGER > _global2localMap;

    /**
     * @brief Build the mapping from global to local numbering of the dof
     */
    void _buildGlobal2LocalMap();

  public:
    /**
     * @typedef EquationNumberingPtr
     * @brief Pointeur intelligent vers un ParallelEquationNumbering
     */
    typedef std::shared_ptr< ParallelEquationNumbering > ParallelEquationNumberingPtr;

    ParallelEquationNumbering();

    ParallelEquationNumbering( const std::string &baseName );

    /**
     * @brief Returns the vector of local to global numbering
     */
    const JeveuxVectorLong getLocalToGlobal() const { return _localToGlobal; };

    /**
     * @brief Returns the vector of the rank owning the local dof number
     */
    const JeveuxVectorLong getLocalToRank() const { return _localToRank; };

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    VectorPairLong getNodeAndComponentIdFromDOF( const bool local = true ) const;

    /**
     * @brief Returns a vector with node index and component name for each DOFs
     */
    std::vector< std::pair< ASTERINTEGER, std::string > >
    getNodeAndComponentFromDOF( const bool local = true ) const;

    /**
     * @brief Returns a pair with node index and component name for given DOF
     */
    std::pair< ASTERINTEGER, std::string >
    getNodeAndComponentFromDOF( const ASTERINTEGER dof, const bool local = true ) const;

    /**
     * @brief Returns a pair with node index and component Id for given DOF
     */

    PairLong getNodeAndComponentIdFromDOF( const ASTERINTEGER dof, const bool local = true ) const;
    /**
     * @brief Get The Component Associated To A Given Row
     */
    std::string getComponentFromDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    ASTERINTEGER getNodeFromDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Get (Nodes And Components) and dofs corresponding to component name and list of node
     * groups
     */
    std::pair< std::pair< VectorLong, VectorString >, VectorLong >
    getDOFsWithDescription( const VectorString &, const VectorString &, const bool local = false,
                            const ASTERINTEGER same_rank = PythonBool::None ) const;

    std::pair< std::pair< VectorLong, VectorString >, VectorLong >
    getDOFsWithDescription( const VectorString &, const VectorLong &, const bool local = false,
                            const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Return true if a physical dof is Associated To A Given Row
     */
    bool isPhysicalDOF( const ASTERINTEGER dof, const bool local = false ) const;

    /**
     * @brief Get The total number of Dofs
     */
    ASTERINTEGER getNumberOfDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    VectorLong getPhysicalDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to all Ghost Dof
     */
    VectorLong getGhostDOFs( const bool local = false, const bool lastLayerOnly = false ) const;

    /**
     * @brief Get Rows owned locally (aka not Ghost)
     */
    VectorLong getNoGhostDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    VectorLong getLagrangeDOFs( const bool local = false ) const;

    /**
     * @brief Get Rows Associated to first and second Lagrange Multipliers Dof
     */

    std::map< ASTERINTEGER, VectorLong > getDictOfLagrangeDOFs( const bool local = false ) const;

    /**
     * @brief Get Assigned Components
     */
    VectorString getComponents() const;
    SetLong getComponentsId() const;

    /**
     * @brief Maps between name of components and the number
     */
    std::map< std::string, ASTERINTEGER > getComponentsNameToId() const;

    /**
     * @brief Get the mapping between local ang global numbering of the Dof
     */
    const JeveuxVectorLong getLocalToGlobalMapping() const { return getLocalToGlobal(); };

    /**
     * @brief Return the local number of a global Dof
     * @return Return the local number if the dof if present on the subdomain ; otherwise
     * raise an exception
     */
    ASTERINTEGER globalToLocalDOF( const ASTERINTEGER ) const;

    /**
     * @brief Return the global number of a local Dof
     * @return Return the global number if the dof if present on the subdomain ; otherwise
     * raise an exception
     */
    ASTERINTEGER localToGlobalDOF( const ASTERINTEGER ) const;

    /**
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    bool useLagrangeDOF() const;

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    bool useSingleLagrangeDOF() const;

    bool isParallel() const { return true; };

    bool build();
};

/**
 * @typedef ParallelEquationNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un ParallelEquationNumbering
 */
typedef std::shared_ptr< ParallelEquationNumbering > ParallelEquationNumberingPtr;

#endif
