#ifndef DATAFIELD_H_
#define DATAFIELD_H_

/**
 * @file DataField.h
 * @brief Fichier entete de la classe DataField
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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"

/**
 * @class DataField
 * @brief class which describe a field of data
 * @author Nicolas Sellenet
 */
class DataField : public DSWithCppPickling {
  private:
  public:
    /**
     * @typedef DataFieldPtr
     * @brief Pointeur intelligent vers un DataField
     */
    typedef std::shared_ptr< DataField > DataFieldPtr;

    /**
     * @brief Constructor
     * @param name Jeveux name
     */
    DataField( const std::string name, const std::string type )
        : DSWithCppPickling( name, 19, type ) {};

    /**
     * @brief Constructor
     */
    DataField( const std::string type ) : DSWithCppPickling( 19, type ) {};

    /**
     * @brief Copy Constructor
     * @param other DataField to copy
     */
    DataField( const DataField &other )
        : DSWithCppPickling( other.getName().size(), other.getType() ) {};

    /**
     * @brief Move Constructor
     */
    DataField( DataField &&other ) : DSWithCppPickling( std::move( other ) ) {};

    /**
     * @brief Constructor
     */
    DataField() : DSWithCppPickling( 19, "CHAM_GD" ) {};

    std::string getFieldType() const;

    std::string getFieldScalar() const;

    virtual bool exists() const {
        AS_ABORT( "Not implemented" );
        return false;
    };
};

/**
 * @typedef DataFieldPtrReal
 */
typedef std::shared_ptr< DataField > DataFieldPtr;

#endif /* DATAFIELD_H_ */
