#ifndef GENERALIZEDLOAD_H_
#define GENERALIZEDLOAD_H_

/**
 * @file GeneralizedLoad.h
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

#include "aster_fort_superv.h"

#include "DataStructures/DataStructure.h"
#include "LinearAlgebra/GeneralizedAssemblyVector.h"
#include "Numbering/GeneralizedDOFNumbering.h"

/**
 * @class GeneralizedLoad
 * @brief Define a generalized load
 */
class GeneralizedLoad : public DataStructure {
  private:
    /** @brief generalized dof_numbering */
    GeneralizedDOFNumberingPtr _nume;
    /** @brief attributes for multi-point constraints */
    std::vector< std::tuple< VectorInt, VectorReal, double > > _mpcs;

  public:
    /**
     * @typedef GeneralizedLoadPtr
     * @brief Pointeur intelligent vers un GeneralizedLoad
     */
    typedef std::shared_ptr< GeneralizedLoad > GeneralizedLoadPtr;

    /**
     * @brief Constructor
     */
    GeneralizedLoad() : GeneralizedLoad( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructor
     */
    GeneralizedLoad( const std::string name )
        : DataStructure( name, 8, "CHAR_GENE" ), _nume( nullptr ) {};

    /**
     * @brief Setter for GeneralizedDOFNumbering
     */
    void setDOFNumbering( const GeneralizedDOFNumberingPtr nume ) { _nume = nume; }

    /**
     * @brief Getter for GeneralizedDOFNumbering
     */
    const GeneralizedDOFNumberingPtr getDOFNumbering() const { return _nume; }

    /**
     * @brief Getter for multi-point constraints
     */
    const auto &getMPCs() const { return _mpcs; }

    /**
     * @brief Setter for multi-point constraints
     */
    void setMPCs( std::vector< std::tuple< VectorInt, VectorReal, double > > &mpcs );

    /**
     * @brief return empty RHS
     */
    const GeneralizedAssemblyVectorRealPtr getAssemblyVector() const;
};

#endif /* GENERALIZEDLOAD_H_ */
