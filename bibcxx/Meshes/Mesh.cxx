/**
 * @file Mesh.cxx
 * @brief Implementation de Mesh
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

#include "Meshes/Mesh.h"

#include "aster_fort_mesh.h"

#include "IOManager/AsterToMedWriter.h"
#include "Utilities/Tools.h"

bool Mesh::readAsterFile( const std::filesystem::path &fileName ) {
    readMeshFile( fileName, "ASTER", 0 );
    return true;
}

bool Mesh::readGibiFile( const std::filesystem::path &fileName ) {
    readMeshFile( fileName, "GIBI", 0 );
    return true;
}

bool Mesh::readGmshFile( const std::filesystem::path &fileName ) {
    readMeshFile( fileName, "GMSH", 0 );
    return true;
}

bool Mesh::hasGroupOfCells( const std::string &name, const bool ) const {
    return _groupsOfCells->contains( name );
}

bool Mesh::hasGroupOfNodes( const std::string &name, const bool ) const {
    return _groupsOfNodes->contains( name );
}

VectorString Mesh::getGroupsOfCells( const bool ) const {
    ASTERINTEGER size = _nameOfGrpCells->size();
    VectorString names;
    for ( auto i = 0; i < size; i++ ) {
        names.push_back( strip( _nameOfGrpCells->getStringFromIndex( i + 1 ) ) );
    }
    return names;
}

VectorString Mesh::getGroupsOfNodes( const bool ) const {
    ASTERINTEGER size = _nameOfGrpNodes->size();
    VectorString names;
    for ( auto i = 0; i < size; i++ ) {
        names.push_back( strip( _nameOfGrpNodes->getStringFromIndex( i + 1 ) ) );
    }
    return names;
}

void Mesh::setGroupOfCells( const std::string &name, const VectorLong &cell_ids ) {
    if ( !name.empty() && !cell_ids.empty() ) {
        ASTERLOGICAL isAdded = false;
        const auto name_s = ljust( strip( name ), 24, ' ' );
        auto cell_ids_u = unique( cell_ids );
        std::for_each( cell_ids_u.begin(), cell_ids_u.end(), []( ASTERINTEGER &d ) { d += 1; } );
        ASTERINTEGER size = cell_ids_u.size(), un = 1;
        CALLO_ADDGROUPELEM( getName(), &un );
        CALLO_ADDGRPMA( getName(), name_s, cell_ids_u.data(), &size, (ASTERLOGICAL *)&isAdded );
        if ( isAdded ) {
            _groupsOfCells->build( true );
        }
    }
};

void Mesh::setGroupOfNodes( const std::string &name, const VectorLong &node_ids,
                            const bool localNumbering ) {
    if ( !name.empty() && !node_ids.empty() ) {
        ASTERLOGICAL isAdded = false;
        const auto name_s = ljust( strip( name ), 24, ' ' );
        auto node_ids_u = unique( node_ids );
        std::for_each( node_ids_u.begin(), node_ids_u.end(), []( ASTERINTEGER &d ) { d += 1; } );
        ASTERINTEGER size = node_ids_u.size(), un = 1;
        CALLO_ADDGROUPNODE( getName(), &un );
        CALLO_ADDGRPNO( getName(), name_s, node_ids_u.data(), &size, (ASTERLOGICAL *)&isAdded );
        if ( isAdded ) {
            _groupsOfNodes->build( true );
        }
    }
};

VectorLong Mesh::getCells( const std::string name ) const {
    return getCells( VectorString( { name } ) );
}

VectorLong Mesh::getCells( const VectorString &names ) const {

    if ( names.empty() ) {
        return irange( (ASTERINTEGER)0, (ASTERINTEGER)( getNumberOfCells() - 1 ) );
    }

    std::vector< VectorLong > cells;
    cells.reserve( names.size() );

    for ( auto &name : names ) {
        if ( hasGroupOfCells( name ) ) {
            cells.push_back( ( *_groupsOfCells )[name]->toVector() );
        }
    }

    auto all_cells = unique( concatenate( cells ) );
    for ( auto &cell : all_cells ) {
        cell -= 1;
    }

    return all_cells;
}

VectorLong Mesh::getNodes( const VectorString &names, const bool, const ASTERINTEGER ) const {

    if ( names.empty() ) {
        return irange( (ASTERINTEGER)0, (ASTERINTEGER)( getNumberOfNodes() - 1 ) );
    }

    std::vector< VectorLong > nodes;
    nodes.reserve( names.size() );

    for ( auto &name : names ) {
        if ( hasGroupOfNodes( name ) ) {
            nodes.push_back( ( *_groupsOfNodes )[name]->toVector() );
        }
    }

    if ( nodes.empty() ) {
        return VectorLong();
    }

    VectorLong all_nodes;

    if ( nodes.size() == 1 ) {
        all_nodes = nodes[0];
    } else {
        all_nodes = unique( concatenate( nodes ) );
    }

    for ( auto &node : all_nodes ) {
        node -= 1;
    }

    return all_nodes;
}

VectorLong Mesh::getNodes( const std::string name, const bool, const ASTERINTEGER ) const {
    if ( name.empty() ) {
        return getNodes( VectorString() );
    }

    return getNodes( VectorString( { name } ) );
}

VectorLong Mesh::getNodesFromCells( const VectorLong &cells, const bool,
                                    const ASTERINTEGER ) const {

    if ( cells.empty() ) {
        return VectorLong();
    }

    CALL_JEMARQ();

    const auto &connecExp = getConnectivityExplorer();

    SetLong nodes;

    for ( auto &cellId : cells ) {
        const auto cell = connecExp[cellId];
        for ( auto &node : cell ) {
            auto ret = nodes.insert( node );
        }
    }

    CALL_JEDEMA();
    return VectorLong( nodes.begin(), nodes.end() );
}

VectorLong Mesh::getNodesFromCells( const VectorString &names, const bool,
                                    const ASTERINTEGER ) const {
    return getNodesFromCells( this->getCells( names ) );
};

VectorLong Mesh::getNodesFromCells( const std::string name, const bool, const ASTERINTEGER ) const {
    return getNodesFromCells( this->getCells( name ) );
};

bool Mesh::isQuadratic( const bool local ) const {
    CALL_JEMARQ();

    auto cellsType = getMedCellsTypes();
    cellsType->updateValuePointer();

    for ( auto &cellType : cellsType ) {
        if ( cellType == 103 || cellType == 104 || cellType == 206 || cellType == 207 ||
             cellType == 208 || cellType == 209 || cellType == 310 || cellType == 315 ||
             cellType == 318 || cellType == 313 || cellType == 320 || cellType == 327 ) {
            CALL_JEDEMA();
            return true;
        }
    }

    CALL_JEDEMA();

    return false;
}

MeshPtr Mesh::fix( const bool remove_orphan, const bool positive_measure, const bool outward_normal,
                   const bool double_nodes, const bool double_cells, const ASTERDOUBLE tole,
                   const ASTERINTEGER info ) {
    auto mesh_out = std::make_shared< Mesh >();
    ASTERINTEGER inf = info, fro, fpv, fon, fdn, fdc;
    fro = static_cast< int >( remove_orphan );
    fpv = static_cast< int >( positive_measure );
    fon = static_cast< int >( outward_normal );
    fdn = static_cast< int >( double_nodes );
    fdc = static_cast< int >( double_cells );
    CALL_FIX_MESH( getName(), mesh_out->getName(), &fro, &fpv, &fon, &fdn, &fdc, &tole, &inf );
    mesh_out->build();
    return mesh_out;
}

MeshPtr Mesh::getOctreeMesh( const ASTERINTEGER nb_max_pt, const ASTERINTEGER nb_max_level ) {
    auto mesh_out = std::make_shared< Mesh >();
    ASTERINTEGER max_pt = nb_max_pt, max_level = nb_max_level;
    CALL_EXPORT_OCTREE( getName(), mesh_out->getName(), &max_pt, &max_level );
    mesh_out->build();
    return mesh_out;
}

MeshPtr Mesh::convertToLinear( const ASTERINTEGER info ) {
    auto mesh_out = std::make_shared< Mesh >();
    ASTERINTEGER un = 1, inf = info;
    CALL_CMBQBQ( getName(), mesh_out->getName(), &un, &inf );
    mesh_out->build();
    return mesh_out;
};

MeshPtr Mesh::convertToQuadratic( const ASTERINTEGER info ) {
    auto mesh_out = std::make_shared< Mesh >();
    ASTERINTEGER deux = 2, inf = info;
    CALL_CMBQBQ( getName(), mesh_out->getName(), &deux, &inf );
    mesh_out->build();
    return mesh_out;
};

MeshPtr Mesh::convertToBiQuadratic( const ASTERINTEGER info ) {
    auto mesh_out = std::make_shared< Mesh >();
    ASTERINTEGER trois = 3, inf = info;
    CALL_CMBQBQ( getName(), mesh_out->getName(), &trois, &inf );
    mesh_out->build();
    return mesh_out;
};

MeshPtr Mesh::convertToCubic( const ASTERINTEGER info ) {
    auto mesh_out = std::make_shared< Mesh >();
    ASTERINTEGER quatre = 4, inf = info;
    CALL_CMBQBQ( getName(), mesh_out->getName(), &quatre, &inf );
    mesh_out->build();
    return mesh_out;
};

void Mesh::addNodeLabels( const VectorString &labels ) {
    const auto &size = labels.size();
    if ( getNumberOfNodes() != size ) {
        throw std::runtime_error( "Bad label size" );
    }
    if ( _nameOfNodes->exists() ) {
        _nameOfNodes->deallocate();
    }
    _nameOfNodes->allocate( size );
    for ( int i = 0; i < size; ++i ) {
        _nameOfNodes->add( i + 1, labels[i] );
    }
}

void Mesh::addCellLabels( const VectorString &labels ) {
    const auto &size = labels.size();
    if ( getNumberOfCells() != size ) {
        throw std::runtime_error( "Bad label size" );
    }
    if ( _nameOfCells->exists() ) {
        _nameOfCells->deallocate();
    }
    _nameOfCells->allocate( size );
    for ( int i = 0; i < size; ++i ) {
        _nameOfCells->add( i + 1, labels[i] );
    }
}

bool Mesh::printMedFile( const std::filesystem::path &fileName, bool local,
                         std::array< int, 3 > version, const std::string &mode ) const {
#ifdef ASTER_HAVE_MED
    auto aTMWriter = AsterToMedWriter();
    return aTMWriter.printMesh( *this, fileName, local, "", mode );
#else
    return false;
#endif
};
