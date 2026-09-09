#ifndef GENERALIZEDASSEMBLYMATRIX_H_
#define GENERALIZEDASSEMBLYMATRIX_H_

/**
 * @file GeneralizedAssemblyMatrix.h
 * @brief Fichier entete de la classe GeneralizedAssemblyMatrix
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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Numbering/ForwardGeneralizedDOFNumbering.h"
#include "Results/ForwardGeneralizedModeResult.h"
#include "Results/ForwardModeResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class GenericGeneralizedAssemblyMatrix
 * @brief Cette classe correspond a un matr_asse_gene
 * @author Nicolas Sellenet
 */
class GenericGeneralizedAssemblyMatrix : public DataStructure {
  private:
    /** @brief Objet Jeveux '.REFE' */
    JeveuxVectorChar24 _refe;
    /** @brief GeneralizedModeResult */
    ForwardGeneralizedModeResultPtr _geneModeC;

  protected:
    /** @brief Objet Jeveux '.DESC' */
    JeveuxVectorLong _desc;
    /** @brief ModeResult */
    ForwardModeResultPtr _mecaModeC;
    /** @brief GeneralizedDOFNumbering */
    ForwardGeneralizedDOFNumberingPtr _dofNum;
    /** @brief V K24 '.REFA' */
    JeveuxVectorChar24 _refa;

  public:
    /**
     * @typedef GeneralizedAssemblyMatrixPtr
     * @brief Pointeur intelligent vers un GeneralizedAssemblyMatrix
     */
    typedef std::shared_ptr< GenericGeneralizedAssemblyMatrix > GenericGeneralizedAssemblyMatrixPtr;

    /**
     * @brief Constructeur
     */
    GenericGeneralizedAssemblyMatrix( const std::string name )
        : DataStructure( name, 19, "MATR_ASSE_GENE" ),
          _desc( JeveuxVectorLong( getName() + ".DESC" ) ),
          _refe( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _refa( JeveuxVectorChar24( getName() + ".REFA" ) ),
          _dofNum( nullptr ),
          _mecaModeC( nullptr ),
          _geneModeC( nullptr ) {};

    /**
     * @brief Get GeneralizedDOFNumbering
     */
    GeneralizedDOFNumberingPtr getGeneralizedDOFNumbering() {
        if ( _dofNum.isSet() )
            return _dofNum.getPointer();
        return GeneralizedDOFNumberingPtr( nullptr );
    };

    /**
     * @brief Get GeneralizedModeResult
     */
    GeneralizedModeResultPtr getModalBasisFromGeneralizedModeResult() {
        if ( _geneModeC.isSet() )
            return _geneModeC.getPointer();
        return GeneralizedModeResultPtr( nullptr );
    };

    /**
     * @brief Get ModeResult
     */
    ModeResultPtr getModalBasisFromModeResult() {
        if ( _mecaModeC.isSet() )
            return _mecaModeC.getPointer();
        return ModeResultPtr( nullptr );
    };

    /**
     * @brief Set GeneralizedDOFNumbering
     */
    bool setGeneralizedDOFNumbering( const GeneralizedDOFNumberingPtr &dofNum ) {
        if ( dofNum != nullptr ) {
            _dofNum = dofNum;
            return true;
        }
        return false;
    };

    /**
     * @brief Set GeneralizedModeResult
     */
    bool setModalBasis( const GeneralizedModeResultPtr &mecaModeC ) {
        if ( mecaModeC != nullptr ) {
            _geneModeC = mecaModeC;
            _mecaModeC = nullptr;
            return true;
        }
        return false;
    };

    /**
     * @brief Set ModeResult
     */
    bool setModalBasis( const ModeResultPtr &mecaModeC ) {
        if ( mecaModeC != nullptr ) {
            _mecaModeC = mecaModeC;
            _geneModeC = nullptr;
            return true;
        }
        return false;
    };

    bool isDiagonal() const {
        _desc->updateValuePointer();
        return ( *_desc )[2] == 1;
    }

    bool isDense() const {
        _desc->updateValuePointer();
        return ( *_desc )[2] == 2;
    }

    ASTERINTEGER size() const {
        _desc->updateValuePointer();
        return ( *_desc )[1];
    }

    bool exists() const { return _desc.exists(); }
};

/**
 * @class GeneralizedAssemblyMatrix
 * @brief Cette classe correspond a un matr_asse_gene
 * @author Nicolas Sellenet
 */
template < class ValueType >
class GeneralizedAssemblyMatrix : public GenericGeneralizedAssemblyMatrix {
  private:
    /** @brief Objet Jeveux '.VALM' */
    JeveuxCollection< ValueType > _matrixValues;
    /** @brief V Jeveux C or I '.CONL' */
    JeveuxVector< ValueType > _conl;

    /**
     * @brief definir le type
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERDOUBLE >::value, void >::type
    setMatrixType() {
        setType( "MATR_ASSE_GENE_R" );
    };

    /**
     * @brief definir le type
     */
    template < class type = ValueType >
    typename std::enable_if< std::is_same< type, ASTERCOMPLEX >::value, void >::type
    setMatrixType() {
        setType( "MATR_ASSE_GENE_C" );
    };

  public:
    /**
     * @typedef GeneralizedAssemblyMatrixPtr
     * @brief Pointeur intelligent vers un GeneralizedAssemblyMatrix
     */
    typedef std::shared_ptr< GeneralizedAssemblyMatrix< ValueType > > GeneralizedAssemblyMatrixPtr;

    /**
     * @brief Constructeur
     */
    GeneralizedAssemblyMatrix() : GeneralizedAssemblyMatrix( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     */
    GeneralizedAssemblyMatrix( const std::string name )
        : GenericGeneralizedAssemblyMatrix( name ),
          _matrixValues( JeveuxCollection< ValueType >( getName() + ".VALM" ) ),
          _conl( JeveuxVector< ValueType >( getName() + ".CONL" ) ) {
        GeneralizedAssemblyMatrix< ValueType >::setMatrixType();
    };

    JeveuxCollection< ValueType > getValues() const { return _matrixValues; };

    /**
     * @brief Allocate the matrix
     */
    void allocate( bool isSymmetric ) {
        if ( _mecaModeC.getPointer() == nullptr )
            throw std::runtime_error( "Unable to allocate, ModalBasis is not set" );
        if ( _dofNum.getPointer() == nullptr )
            throw std::runtime_error( "Unable to allocate, GeneralizedDOFNumbering is not set" );
        ASTERINTEGER size = _mecaModeC.getNumberOfIndexes();
        _desc->allocate( 3 );
        ( *_desc )[0] = 2;
        ( *_desc )[1] = size;
        ( *_desc )[2] = 2;
        //_desc->updateValuePointer();
        if ( isSymmetric )
            _matrixValues->allocate( 1 );
        else
            _matrixValues->allocate( 2 );
        _matrixValues->allocateObject( 1, size * ( size + 1 ) / 2 );
        if ( !isSymmetric )
            _matrixValues->allocateObject( 2, size * ( size + 1 ) / 2 );
        //_matrixValues->build();
        _refa->allocate( 20 );
        ( *_refa )[0] = _mecaModeC.getName();
        ( *_refa )[1] = _dofNum.getName();
        if ( isSymmetric )
            ( *_refa )[8] = "MS";
        else
            ( *_refa )[8] = "MR";
        ( *_refa )[10] = "MPI_COMPLET";
    };

    std::vector< ValueType > getUpperValues() const {
        _matrixValues->build();
        _matrixValues->updateValuePointer();
        return ( *_matrixValues )[1]->toVector();
    };

    std::vector< ValueType > getLowerValues() const {
        _matrixValues->build();
        _matrixValues->updateValuePointer();
        return ( *_matrixValues )[2]->toVector();
    };

    void setUpperValues( const std::vector< ValueType > &values ) {
        _matrixValues->build();
        _matrixValues->updateValuePointer();
        ( *_matrixValues )[1]->setValues( values );
    };

    void setLowerValues( const std::vector< ValueType > &values ) {
        _matrixValues->build();
        _matrixValues->updateValuePointer();
        ( *_matrixValues )[2]->setValues( values );
    };

    bool isSymmetric() const {
        _matrixValues->build();
        return _matrixValues->size() == 1;
    };

    bool build() { return _matrixValues->build( true ); }
};

/** @typedef Definition d'une matrice assemblee généralisée de double */
typedef GeneralizedAssemblyMatrix< ASTERDOUBLE > GeneralizedAssemblyMatrixReal;
/** @typedef Definition d'une matrice assemblee généralisée de complexe */
typedef GeneralizedAssemblyMatrix< ASTERCOMPLEX > GeneralizedAssemblyMatrixComplex;

/**
 * @typedef GenericGeneralizedAssemblyMatrixPtr
 * @brief Pointeur intelligent vers un GenericGeneralizedAssemblyMatrix
 */
typedef std::shared_ptr< GenericGeneralizedAssemblyMatrix > GenericGeneralizedAssemblyMatrixPtr;

/**
 * @typedef GeneralizedAssemblyMatrixRealPtr
 * @brief Pointeur intelligent vers un GeneralizedAssemblyMatrixReal
 */
typedef std::shared_ptr< GeneralizedAssemblyMatrixReal > GeneralizedAssemblyMatrixRealPtr;

/**
 * @typedef GeneralizedAssemblyMatrixComplexPtr
 * @brief Pointeur intelligent vers un GeneralizedAssemblyMatrixComplex
 */
typedef std::shared_ptr< GeneralizedAssemblyMatrixComplex > GeneralizedAssemblyMatrixComplexPtr;

#endif /* GENERALIZEDASSEMBLYMATRIX_H_ */
