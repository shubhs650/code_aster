#ifndef DISCRETEPROBLEM_H_
#define DISCRETEPROBLEM_H_

/**
 * @file DiscreteComputation.h
 * @brief Header of class DiscreteComputation
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

#include "Behaviours/BehaviourProperty.h"
#include "Discretization/Calcul.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/ElementaryVector.h"
#include "Studies/PhysicalProblem.h"

/**
 * @class DiscreteComputation
 * @brief Compute discrete operators (vectors and matrices)
 * All methods are in the same header but implementation is splitted in several files
 * - DiscreteComputation.cxx for generic methods
 * - DiscreteComputationVector.cxx to compute vector thar are independent of physics
 * - DiscreteComputationMechanicalVector.cxx to compute vector for mechanics
 * - DiscreteComputationMechanicalMatrix.cxx to compute matrix for mechanics
 * - DiscreteComputationThermalMatrix.cxx to compute matrix for thermic
 * - DiscreteComputationThermalVector.cxx to compute vector for thermic
 */
class DiscreteComputation {
  private:
    /** @brief Physical problem */
    PhysicalProblemPtr _phys_problem;

    /** @brief Create time field */
    ConstantFieldOnCellsRealPtr createTimeField( const ASTERDOUBLE time_curr,
                                                 const ASTERDOUBLE time_step = 0.0,
                                                 const ASTERDOUBLE theta = 0.0 ) const;

    /** @brief Create damping fluid field */
    ConstantFieldOnCellsLongPtr createDampingFluidField( const ASTERINTEGER damping,
                                                         const ASTERINTEGER onde_flui ) const;

    /** @brief Create wave type fluid field (for IMPE_MECA and AMOR_ACOU option) */
    ConstantFieldOnCellsLongPtr createWaveTypeFluidField( const ASTERINTEGER onde_flui ) const;

    /** @brief Preparation for non-linear computations */
    CalculPtr createCalculForNonLinear( const std::string option, const ASTERDOUBLE &time_prev,
                                        const ASTERDOUBLE &time_step,
                                        const FieldOnCellsRealPtr _externVarFieldPrev,
                                        const FieldOnCellsRealPtr _externVarFieldCurr,
                                        const VectorString &groupOfCells = VectorString() ) const;

    /** @brief Compute B elementary matrices fo dualized boundary conditions */
    void baseDualElasticStiffnessMatrix( CalculPtr &calcul,
                                         ElementaryMatrixDisplacementRealPtr &elemMatr ) const;

    /** @brief Compute B elementary matrices fo dualized thermal boundary conditions */
    void baseDualLinearConductivityMatrix( CalculPtr &calcul,
                                           ElementaryMatrixTemperatureRealPtr &elemMatr ) const;

    /** @brief Compute B elementary matrices fo dualized acoustic boundary conditions */
    void baseDualAcousticMatrix( CalculPtr &calcul,
                                 ElementaryMatrixPressureComplexPtr &elemMatr ) const;
    /**
     * @brief Compute Dirichlet reaction vector B^T * \lambda for mechanical case
     * @param lagr_curr Field on nodes for Lagrange multipliers
     * @return Nodal field for Dirichlet reaction vector
     */
    FieldOnNodesRealPtr dualMechanicalVector( FieldOnNodesRealPtr lagr_curr ) const;

    /**
     * @brief Compute Dirichlet reaction vector B^T * \lambda for thermal case
     * @param lagr_curr Field on nodes for Lagrange multipliers
     * @return Nodal field for Dirichlet reaction vector
     */
    FieldOnNodesRealPtr dualThermalVector( FieldOnNodesRealPtr lagr_curr ) const;

  public:
    /** @typedef DiscreteComputationPtr */
    typedef std::shared_ptr< DiscreteComputation > DiscreteComputationPtr;

    /** @brief Default constructor disabled */
    DiscreteComputation( void ) = delete;

    /**
     * @brief Constructor
     * @param PhysicalProblemPtr study
     */
    DiscreteComputation( const PhysicalProblemPtr &currPhysProblem )
        : _phys_problem( currPhysProblem ) {};

    /** @brief Destructor */
    ~DiscreteComputation() {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    DiscreteComputation( const py::tuple &tup )
        : DiscreteComputation( tup[0].cast< PhysicalProblemPtr >() ) {};
    py::tuple _getState() const { return py::make_tuple( _phys_problem ); };

    /** @brief Compute nodal field for external state variables RHS */
    std::variant< ElementaryVectorRealPtr, FieldOnNodesRealPtr > getExternalStateVariablesForces(
        const ASTERDOUBLE time_curr, const FieldOnCellsRealPtr varc_curr,
        const FieldOnCellsRealPtr varc_prev = nullptr,
        const FieldOnCellsRealPtr vari_curr = nullptr,
        const FieldOnCellsRealPtr stress_prev = nullptr, const ASTERINTEGER mode_fourier = 0,
        const bool assembly = true, const FieldOnCellsLongPtr maskField = nullptr ) const;

    /**
     * @brief Compute imposed displacement U_impo with Lagrange
     * @param time_curr Time
     * @return Nodal field for imposed displacement
     */
    std::variant< ElementaryVectorDisplacementRealPtr, FieldOnNodesRealPtr >
    getMechanicalImposedDualBC( const ASTERDOUBLE time_curr = 0.0,
                                const bool assembly = true ) const;

    std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
    getThermalImposedDualBC( const ASTERDOUBLE time_curr = 0.0, const bool assembly = true ) const;

    std::variant< ElementaryVectorPressureComplexPtr, FieldOnNodesComplexPtr >
    getAcousticImposedDualBC( const bool assembly = true ) const;

    /**
     * @brief Compute Dirichlet reaction vector B^T * \lambda
     * @param lagr_curr Field on nodes for Lagrange multipliers
     * @return Nodal field for Dirichlet reaction vector
     */
    FieldOnNodesRealPtr getDualForces( FieldOnNodesRealPtr lagr_curr ) const;

    /**
     * @brief Compute Dirichlet differential imposed dualized displacement B * U
     * @return Nodal field for dualized differential displacement
     */
    FieldOnNodesRealPtr getDifferentialDualDisplacement() const;

    /**
     * @brief Compute Dirichlet imposed dualized displacement B * U
     * @return Nodal field for dualized displacement
     */
    FieldOnNodesRealPtr getDualDisplacement( FieldOnNodesRealPtr disp_curr,
                                             ASTERDOUBLE scaling = 1.0 ) const;

    /**
     * @brief Compute Dirichlet imposed dualized displacement B * U
     * @return Nodal field for dualized displacement
     */
    FieldOnNodesRealPtr getDualTemperature( FieldOnNodesRealPtr temp_curr,
                                            ASTERDOUBLE scaling = 1.0 ) const;

    /**
     * @brief Compute Dirichlet imposed dualized primal B * U
     * @return Nodal field for dualized displacement
     */
    FieldOnNodesRealPtr getDualPrimal( FieldOnNodesRealPtr disp_curr,
                                       ASTERDOUBLE scaling = 1.0 ) const;
    /**
     * @brief Compute Neumann loads
     * @param TimeParameters Parameters for time
     * @return Nodal field for Neumann loads
     */
    std::variant< ElementaryVectorDisplacementRealPtr, FieldOnNodesRealPtr >
    getMechanicalNeumannForces( const ASTERDOUBLE time_curr = 0.0,
                                const ASTERDOUBLE time_step = 0.0, const ASTERDOUBLE theta = 1.0,
                                const ASTERINTEGER modeFourier = 0,
                                const FieldOnCellsRealPtr varc_curr = nullptr,
                                const bool assembly = true ) const;

    FieldOnNodesRealPtr getMechanicalForces( const ASTERDOUBLE time_curr = 0.0,
                                             const ASTERDOUBLE time_step = 0.0,
                                             const ASTERDOUBLE theta = 1.0,
                                             const ASTERINTEGER modeFourier = 0,
                                             const FieldOnCellsRealPtr varc_curr = nullptr ) const;

    std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
    getThermalNeumannForces( const ASTERDOUBLE time_curr = 0.0, const bool assembly = true ) const;

    std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
    getThermalNonLinearNeumannForces( const FieldOnNodesRealPtr temp_curr,
                                      const ASTERDOUBLE time_curr,
                                      const bool assembly = true ) const;

    std::variant< ElementaryVectorPressureComplexPtr, FieldOnNodesComplexPtr >
    getAcousticNeumannForces( const bool assembly = true ) const;

    /**
     * @brief Compute volumetric loads
     */
    std::variant< ElementaryVectorDisplacementRealPtr, FieldOnNodesRealPtr >
    getMechanicalVolumetricForces( const ASTERDOUBLE time_curr = 0.0,
                                   const ASTERDOUBLE time_step = 0.0, const ASTERDOUBLE theta = 1.0,
                                   const ASTERINTEGER modeFourier = 0,
                                   const FieldOnCellsRealPtr varc_curr = nullptr,
                                   const bool assembly = true ) const;

    std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
    getThermalVolumetricForces( const ASTERDOUBLE time_curr = 0.0,
                                const FieldOnCellsRealPtr varc_curr = nullptr,
                                const bool assembly = true ) const;

    std::variant< ElementaryVectorPressureComplexPtr, FieldOnNodesComplexPtr >
    getAcousticVolumetricForces( const bool assembly = true ) const;

    std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
    getThermalNonLinearVolumetricForces( const FieldOnNodesRealPtr temp_curr,
                                         const ASTERDOUBLE time_curr,
                                         const bool assembly = true ) const;

    /** @brief Compute nodal forces */
    std::variant< ElementaryVectorDisplacementRealPtr, FieldOnNodesRealPtr >
    getMechanicalNodalForces( const FieldOnCellsRealPtr stress,
                              const FieldOnNodesRealPtr disp = nullptr,
                              const ASTERINTEGER modeFourier = 0,
                              const FieldOnCellsRealPtr varc_curr = nullptr,
                              const ConstantFieldOnCellsChar16Ptr behaviourMap = nullptr,
                              const VectorString &groupOfCells = VectorString(),
                              const bool assembly = true ) const;

    /** @brief Compute reaction forces */
    FieldOnNodesRealPtr getMechanicalReactionForces(
        const FieldOnCellsRealPtr stress, const FieldOnNodesRealPtr disp = nullptr,
        const ASTERDOUBLE time_prev = 0.0, const ASTERDOUBLE time_curr = 0.0,
        const ASTERDOUBLE theta = 1.0, const ASTERINTEGER modeFourier = 0,
        const FieldOnCellsRealPtr varc_curr = nullptr,
        const ConstantFieldOnCellsChar16Ptr behaviourMap = nullptr ) const;

    /**
     * @brief Compute elementary matrices for mechanical stiffness (RIGI_MECA)
     */
    ElementaryMatrixDisplacementRealPtr getElasticStiffnessMatrix(
        const ASTERDOUBLE &time_curr = 0.0, const ASTERINTEGER &modeFourier = 0,
        const FieldOnCellsRealPtr varc_curr = nullptr,
        const VectorString &groupOfCells = VectorString(), const bool &with_dual = true ) const;

    /**
     * @brief Compute elementary matrices for mechanical stiffness (RIGI_FLUI_STRU)
     */
    ElementaryMatrixDisplacementRealPtr
    getFluidStructureStiffnessMatrix( const ASTERINTEGER &modeFourier = 0,
                                      const FieldOnCellsRealPtr varc_curr = nullptr,
                                      const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for rotational mechanical stiffness (RIGI_ROTA)
     */
    ElementaryMatrixDisplacementRealPtr getGeometricStiffnessMatrix(
        const FieldOnCellsRealPtr sief_elga, const FieldOnCellsRealPtr strx_elga = nullptr,
        const FieldOnNodesRealPtr displ = nullptr, const ASTERINTEGER &modeFourier = 0,
        const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute follower-load stiffness contributions for mechanical buckling.
     */
    ElementaryMatrixDisplacementRealPtr
    getMechanicalFollowerLoadStiffnessMatrix( const FieldOnNodesRealPtr displ = nullptr,
                                              const ASTERDOUBLE &time_curr = 0.0,
                                              const ASTERINTEGER &modeFourier = 0,
                                              const VectorString &groupOfCells =
                                                  VectorString() ) const;

    /**
     * @brief Compute elementary matrices for rotational mechanical stiffness (RIGI_ROTA)
     */
    ElementaryMatrixDisplacementRealPtr
    getRotationalStiffnessMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for rotational mechanical stiffness (RIGI_GYRO)
     */
    ElementaryMatrixDisplacementRealPtr
    getGyroscopicStiffnessMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for rotational mechanical damping (MECA_GYRO)
     */
    ElementaryMatrixDisplacementRealPtr
    getGyroscopicDampingMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for thermal model (RIGI_THER)
     */
    ElementaryMatrixTemperatureRealPtr
    getLinearConductivityMatrix( const ASTERDOUBLE time_curr, const ASTERINTEGER &modeFourier = -1,
                                 const FieldOnCellsRealPtr varc_curr = nullptr,
                                 const VectorString &groupOfCells = VectorString(),
                                 const bool &with_dual = true ) const;

    /**
     * @brief Compute elementary matrices for thermal model (RIGI_THER_TANG)
     */
    ElementaryMatrixTemperatureRealPtr getTangentConductivityMatrix(
        const FieldOnNodesRealPtr temp_prev, const FieldOnNodesRealPtr temp_step,
        const FieldOnCellsRealPtr varc_curr = nullptr,
        const VectorString &groupOfCells = VectorString(), const bool &with_dual = true ) const;

    /**
     * @brief Compute elementary matrices for acoustic model (RIGI_ACOU)
     */
    ElementaryMatrixPressureComplexPtr
    getLinearMobilityMatrix( const VectorString &groupOfCells = VectorString(),
                             const bool &with_dual = true ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    getMechanicalMassMatrix( const bool diagonal, const FieldOnCellsRealPtr varc_curr = nullptr,
                             const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for mechanical stiffness (MASS_FLUI_STRU)
     */
    ElementaryMatrixDisplacementRealPtr
    getFluidStructureMassMatrix( const FieldOnCellsRealPtr varc_curr = nullptr,
                                 const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_THER)
     */
    ElementaryMatrixTemperatureRealPtr
    getLinearCapacityMatrix( const ASTERDOUBLE time_curr,
                             const FieldOnCellsRealPtr varc_curr = nullptr,
                             const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_THER_TANG)
     */
    ElementaryMatrixTemperatureRealPtr
    getTangentCapacityMatrix( const FieldOnNodesRealPtr temp_prev,
                              const FieldOnNodesRealPtr temp_step,
                              const FieldOnCellsRealPtr varc_curr = nullptr,
                              const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for mass matrix (MASS_ACOU)
     */
    ElementaryMatrixPressureComplexPtr
    getCompressibilityMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for damping matrix (AMOR_MECA)
     */
    ElementaryMatrixDisplacementRealPtr getMechanicalDampingMatrix(
        const ElementaryMatrixDisplacementRealPtr &massMatrix = nullptr,
        const ElementaryMatrixDisplacementRealPtr &stiffnessMatrix = nullptr,
        const FieldOnCellsRealPtr varc_curr = nullptr,
        const VectorString &groupOfCells = VectorString(), const ASTERINTEGER &flui_int = 1,
        const ASTERINTEGER &onde_flui = 1 ) const;
    /**
     * @brief Compute elementary matrices for damping matrix (AMOR_ACOU)
     */
    ElementaryMatrixPressureComplexPtr
    getImpedanceMatrix( const ASTERINTEGER &onde_flui = 1 ) const;

    /**
     * @brief Compute third order elementary matrices for absorbing fluid elements (IMPE_MECA)
     */
    ElementaryMatrixDisplacementRealPtr
    getImpedanceBoundaryMatrix( const VectorString &groupOfCells = VectorString(),
                                const ASTERINTEGER &onde_flui = 1 ) const;

    ElementaryMatrixDisplacementRealPtr
    getImpedanceWaveMatrix( const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute elementary matrices for complex rigidity matrix (RIGI_MECA_HYST)
     */
    ElementaryMatrixDisplacementComplexPtr
    getHystereticStiffnessMatrix( const ElementaryMatrixDisplacementRealPtr &stiffnessMatrix,
                                  const FieldOnCellsRealPtr varc_curr = nullptr,
                                  const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute nodal field for kinematic boundary condition
     * @param time_curr Time
     * @return Nodal field for kinematic boundary condition
     */
    template < typename T >
    std::shared_ptr< FieldOnNodes< T > > _getDirichletBC( const ASTERDOUBLE time_curr = 0.0 ) const;
    FieldOnNodesRealPtr _getDirichletDifferentialBC() const;
    FieldOnNodesRealPtr getMechanicalDirichletBC( const ASTERDOUBLE time_curr = 0.0 ) const;
    FieldOnNodesRealPtr getThermalDirichletBC( const ASTERDOUBLE time_curr = 0.0 ) const;
    FieldOnNodesComplexPtr getAcousticDirichletBC( const ASTERDOUBLE time_curr = 0.0 ) const;

    /**
     * @brief Compute nodal field for incremental kinematic boundary condition
     * @param time_curr Time
     * @param disp_curr Current displacement
     * @return Nodal field for incremental kinematic boundary condition
     */
    FieldOnNodesRealPtr getIncrementalDirichletBC( const ASTERDOUBLE &time_curr,
                                                   const FieldOnNodesRealPtr disp_curr ) const;

    /**
     * @brief Compute B elementary matrices for dualized boundary conditions
     * @return Elementary matrices for dualized boundary conditions
     */
    ElementaryMatrixDisplacementRealPtr getDualElasticStiffnessMatrix() const;

    ElementaryMatrixPressureComplexPtr getDualLinearMobilityMatrix() const;

    ElementaryMatrixTemperatureRealPtr getDualLinearConductivityMatrix() const;

    ElementaryMatrixTemperatureRealPtr
    getThermalExchangeMatrix( const ASTERDOUBLE &time_curr ) const;

    /**
     * @brief Get physical problem
     * @return Physical problem
     */
    PhysicalProblemPtr getPhysicalProblem() const { return _phys_problem; };

    std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
    getThermalExchangeForces( const FieldOnNodesRealPtr temp_curr,
                              const ASTERDOUBLE time_curr = 0.0, const bool assembly = true ) const;

    FieldOnNodesRealPtr
    getTransientThermalForces( const ASTERDOUBLE time_curr, const ASTERDOUBLE time_step,
                               const ASTERDOUBLE theta,
                               const FieldOnNodesRealPtr previousPrimalField,
                               const FieldOnCellsRealPtr varc_curr = nullptr ) const;

    std::variant< ElementaryVectorTemperatureRealPtr, FieldOnNodesRealPtr >
    getTransientThermalLoadForces( const ASTERDOUBLE time_curr,
                                   const FieldOnNodesRealPtr temp_prev = nullptr,
                                   const bool assembly = true ) const;

    /**
     * @brief Compute internal forces, stress and internal state variables
     * @return Tuple with 5 objects:
     * field of exitcode
     * error code (integer)
     * internal state variables (VARI_ELGA)
     * Cauchy stress (SIEF_ELGA)
     * field of internal forces (`B^T \sigma`)
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
                FieldOnNodesRealPtr >
    getInternalMechanicalForces( const FieldOnNodesRealPtr displ_prev,
                                 const FieldOnNodesRealPtr displ_step,
                                 const FieldOnCellsRealPtr stress,
                                 const FieldOnCellsRealPtr internVar,
                                 const FieldOnCellsRealPtr internVarIter,
                                 const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                                 const FieldOnCellsRealPtr &varc_prev = nullptr,
                                 const FieldOnCellsRealPtr &varc_curr = nullptr,
                                 const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute internal forces, stress and internal state variables
     * @return Tuple with 5 objects:
     * field of internal forces (`B^T \sigma`)
     */
    std::tuple< ASTERINTEGER, FieldOnCellsRealPtr, FieldOnNodesRealPtr >
    getInternalThermalForces( const FieldOnNodesRealPtr temp_prev,
                              const FieldOnNodesRealPtr temp_step,
                              const FieldOnCellsRealPtr varc_curr = nullptr,
                              const VectorString &groupOfCells = VectorString() ) const;

    // MASS_THER_RESI
    FieldOnNodesRealPtr
    getNonLinearCapacityForces( const FieldOnNodesRealPtr temp_prev,
                                const FieldOnNodesRealPtr temp_step,
                                const FieldOnCellsRealPtr varc_curr = nullptr,
                                const VectorString &groupOfCells = VectorString() ) const;

    ElementaryMatrixTemperatureRealPtr
    getThermalTangentNonLinearVolumetricMatrix( const FieldOnNodesRealPtr temp_curr,
                                                const ASTERDOUBLE time_curr ) const;

    ElementaryMatrixTemperatureRealPtr
    getThermalTangentNonLinearNeumannMatrix( const FieldOnNodesRealPtr temp_curr,
                                             const ASTERDOUBLE time_curr,
                                             const FieldOnCellsRealPtr varc_curr = nullptr ) const;

    /**
     * @brief Compute tangent matrix (not assembled)
     * @return Tuple with 3 objects:
     * field of exitcode
     * error code (integer)
     * elementary tangent matrix
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
    getTangentStiffnessMatrix( const FieldOnNodesRealPtr displ_prev,
                               const FieldOnNodesRealPtr displ_step,
                               const FieldOnCellsRealPtr stress,
                               const FieldOnCellsRealPtr internVar,
                               const FieldOnCellsRealPtr internVarIter,
                               const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                               const FieldOnCellsRealPtr &varc_prev = nullptr,
                               const FieldOnCellsRealPtr &varc_curr = nullptr,
                               const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute tangent prediction matrix (not assembled)
     * @return Tuple with 3 objects:
     * field of exitcode
     * error code (integer)
     * elementary prediction matrix
     */
    std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, ElementaryMatrixDisplacementRealPtr >
    getPredictionTangentStiffnessMatrix( const FieldOnNodesRealPtr displ_prev,
                                         const FieldOnNodesRealPtr displ_step,
                                         const FieldOnCellsRealPtr stress,
                                         const FieldOnCellsRealPtr internVar,
                                         const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step,
                                         const FieldOnCellsRealPtr &varc_prev = nullptr,
                                         const FieldOnCellsRealPtr &varc_curr = nullptr,
                                         const VectorString &groupOfCells = VectorString() ) const;

    /**
     * @brief Compute contact forces
     */
    FieldOnNodesRealPtr
    getContactForces( const MeshCoordinatesFieldPtr geom, const FieldOnNodesRealPtr displ_prev,
                      const FieldOnNodesRealPtr displ_step, const ASTERDOUBLE &time_prev,
                      const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr data,
                      const FieldOnNodesRealPtr coef_cont,
                      const FieldOnNodesRealPtr coef_frot ) const;

    /**
     * @brief Compute contact matrix
     */
    ElementaryMatrixDisplacementRealPtr
    getContactMatrix( const MeshCoordinatesFieldPtr geom, const FieldOnNodesRealPtr displ_prev,
                      const FieldOnNodesRealPtr displ_step, const ASTERDOUBLE &time_prev,
                      const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr data,
                      const FieldOnNodesRealPtr coef_cont,
                      const FieldOnNodesRealPtr coef_frot ) const;

    /**
     * @brief Compute coupling forces
     */

    FieldOnNodesRealPtr getMechanicalCouplingForces( const FieldOnNodesRealPtr displ_prev,
                                                     const FieldOnNodesRealPtr displ_step,
                                                     const ASTERDOUBLE &time_prev,
                                                     const ASTERDOUBLE &time_step ) const;

    ElementaryMatrixDisplacementRealPtr
    getMechanicalCouplingMatrix( const FieldOnNodesRealPtr displ_prev,
                                 const FieldOnNodesRealPtr displ_step, const ASTERDOUBLE &time_prev,
                                 const ASTERDOUBLE &time_step ) const;

    ElementaryMatrixDisplacementRealPtr
    getMechanicalLinearCouplingMatrix( const FieldOnCellsRealPtr &externVar = nullptr ) const;

    /**
     * @brief Compute residual reference (for RESI_REFE_RELA)
     */
    FieldOnNodesRealPtr
    getResidualReference( const std::map< std::string, ASTERDOUBLE > &vale_by_name ) const;
};

using DiscreteComputationPtr = std::shared_ptr< DiscreteComputation >;

#endif /* DISCRETEPROBLEM_H_ */
