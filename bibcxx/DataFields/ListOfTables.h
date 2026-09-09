#ifndef LISTOFTABLES_H_
#define LISTOFTABLES_H_

/**
 * @file ListOfTables.h
 * @brief Fichier entete de la classe ListOfTables
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

#include "DataFields/Table.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ListOfTables
 * @brief Cette classe template permet de definir une table Aster
 * @author Nicolas Sellenet
 */
class ListOfTables {
  private:
    /** @brief Nom de la sd */
    std::string _name;

  protected:
    /** @brief Vecteur Jeveux '.LTNT' */
    JeveuxVectorChar16 _dsId;
    /** @brief Vecteur Jeveux '.LTNS' */
    JeveuxVectorChar24 _dsName;

    /** @typedef std::map d'une chaine et des pointers vers les objets 'Table' */
    typedef std::map< std::string, TablePtr > mapStrTable;
    /** @brief Liste des objets 'Table' */
    mapStrTable _mapTables;

  public:
    /**
     * @typedef ListOfTablesPtr
     * @brief Definition of a smart pointer to a ListOfTables
     */
    typedef std::shared_ptr< ListOfTables > ListOfTablesPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux
     */
    ListOfTables( const std::string &name );

    /**
     * @brief Constructeur
     */
    ListOfTables() : ListOfTables( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Construire une sd_l_tabke à partir d'objet produit dans le Fortran
     * @return true si l'allocation s'est bien passée, false sinon
     */
    bool update_tables();

    /**
     * @brief Return a Table
     * @param id Table identifier
     */
    TablePtr getTable( const std::string id );

    /**
     * @brief Set a Table
     */
    void setTable( const std::string, const TablePtr );
};

#endif /* LISTOFTABLES_H_ */
