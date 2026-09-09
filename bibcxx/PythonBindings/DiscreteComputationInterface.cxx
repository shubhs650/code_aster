/**
 * @file DiscreteComputationInterface.cxx
 * @brief Interface python de DiscreteComputation
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

#include "PythonBindings/DiscreteComputationInterface.h"

#include "aster_pybind.h"

void exportDiscreteComputationToPython( py::module_ &mod ) {

    py::class_< DiscreteComputation, DiscreteComputation::DiscreteComputationPtr >(
        mod, "DiscreteComputation" )
        .def( py::init( &initFactoryPtr< DiscreteComputation, PhysicalProblemPtr > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( define_pickling< DiscreteComputation >() )
        .def( "getMechanicalImposedDualBC", &DiscreteComputation::getMechanicalImposedDualBC,
              R"(
      Return the mechanical imposed nodal BC elementary vector

      Arguments:
            time_curr (float): Current time (default: 0.0)
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorDisplacementReal: imposed dual vector
        )",
              py::arg( "time_curr" ) = 0.0, py::arg( "assembly" ) = true )

        .def( "getThermalImposedDualBC", &DiscreteComputation::getThermalImposedDualBC,
              R"(
      Return the thermal imposed nodal BC elementary vector

      Arguments:
            time_curr (float): Current time (default: 0.0)
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorThermalReal: imposed dual vector
        )",
              py::arg( "time_curr" ) = 0.0, py::arg( "assembly" ) = true )

        .def( "getAcousticImposedDualBC", &DiscreteComputation::getAcousticImposedDualBC,
              R"(
      Return the acoustic imposed nodal BC elementary vector

      Arguments:
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorPressureComplex: imposed dual vector
        )",
              py::arg( "assembly" ) = true )

        .def( "getDualForces", &DiscreteComputation::getDualForces,
              R"(
      Return the imposed displacement assembled vector

      Arguments:
            U (FieldOnNodes): current displacement vector

      Returns:
            FieldOnNodes: dual reaction vector (B^T*lambda)
        )",
              py::arg( "U" ) )
        .def( "getDualPrimal", &DiscreteComputation::getDualPrimal,
              R"(
      Return the Dirichlet load vector

      Arguments:
            U (FieldOnNodes): current displacement vector

      Returns:
            FieldOnNodes: Dirichlet load vector
              )",
              py::arg( "U" ), py::arg( "scaling" ) = 1.0 )
        .def( "getThermalExchangeForces", &DiscreteComputation::getThermalExchangeForces,
              R"(
      Return the elementary thermal Exchange forces vector

      Arguments:
            temp_curr (FieldOnNodesReal): thermal field at current time
            time_curr (float): Current time
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorThermalReal: elementary Exchange forces vector
        )",
              py::arg( "temp_curr" ), py::arg( "time_curr" ) = 0.0, py::arg( "assembly" ) = true )

        .def( "getMechanicalNeumannForces", &DiscreteComputation::getMechanicalNeumannForces,
              R"(
      Return the elementary mechanical Neumann forces vector

      Arguments:
            time_curr (float): Current time
            time_step (float): Time increment
            theta (float): Theta parameter for time-integration
            mode (int) : fourier mode
            varc_curr (FieldOnCellsReal): external state variables at current time
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorDisplacementReal: elementary Neumann forces vector
        )",
              py::arg( "time_curr" ) = 0.0, py::arg( "time_step" ) = 0.0, py::arg( "theta" ) = 1.0,
              py::arg( "mode" ) = 0, py::arg( "varc_curr" ) = nullptr,
              py::arg( "assembly" ) = true )

        .def( "getThermalNeumannForces", &DiscreteComputation::getThermalNeumannForces,
              R"(
      Return the elementary thermal Neumann forces vector

      Arguments:
            time_curr (float): Current time (default: 0.0)
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorThermalReal: elementary Neumann forces vector
        )",
              py::arg( "time_curr" ) = 0.0, py::arg( "assembly" ) = true )

        .def( "getAcousticNeumannForces", &DiscreteComputation::getAcousticNeumannForces,
              R"(
      Return the elementary acoustic Neumann forces vector

      Arguments:
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorPressureComplex: elementary Neumann forces vector
        )",
              py::arg( "assembly" ) = true )

        .def( "getMechanicalVolumetricForces", &DiscreteComputation::getMechanicalVolumetricForces,
              R"(
      Return the elementary mechanical Volumetric forces vector

      Arguments:
            time_curr (float): Current time
            time_step (float): Time increment
            theta (float): Theta parameter for time-integration
            mode (int) : fourier mode
            varc_curr (FieldOnCellsReal): external state variables at current time
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorDisplacementReal: elementary Volumetric forces vector
        )",
              py::arg( "time_curr" ) = 0.0, py::arg( "time_step" ) = 0.0, py::arg( "theta" ) = 1.0,
              py::arg( "mode" ) = 0, py::arg( "varc_curr" ) = nullptr,
              py::arg( "assembly" ) = true )

        .def( "getThermalVolumetricForces", &DiscreteComputation::getThermalVolumetricForces,
              R"(
      Return the elementary thermal Volumetric forces vector

      Arguments:
            time_curr (float): Current time
            varc_curr (FieldOnCellsReal): external state variables at current time
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorThermalReal: elementary Volumetric forces vector
        )",
              py::arg( "time_curr" ) = 0.0, py::arg( "varc_curr" ) = nullptr,
              py::arg( "assembly" ) = true )

        .def( "getAcousticVolumetricForces", &DiscreteComputation::getAcousticVolumetricForces,
              R"(
      Return the elementary acoustic volumetric forces vector

      Arguments:
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorPressureComplex: elementary volumetric forces vector
        )",
              py::arg( "assembly" ) = true )
        .def( "getExternalStateVariablesForces",
              &DiscreteComputation::getExternalStateVariablesForces, R"(
            Compute load from external state variables

            Arguments:
                  time_curr (float): Current time
                  varc_curr (FieldOnCellsReal): external state variables at current time
                  varc_prev (FieldOnCells): external state variables at begin of current time
                  vari_curr (FieldOnCellsReal): internal state variables at current time
                  stress_prev (FieldOnCellsReal): stress at begin of current time
                  mode (int): fourier mode
                  assembly (bool) : assemble or not
                  mask (FieldOnCellsLongPtr): mask to assemble

            Returns:
                  FieldOnNodes: load from external state variables
            )",
              py::arg( "time_curr" ), py::arg( "varc_curr" ), py::arg( "varc_prev" ) = nullptr,
              py::arg( "vari_curr" ) = nullptr, py::arg( "stress_prev" ) = nullptr,
              py::arg( "mode" ) = 0, py::arg( "assembly" ) = true, py::arg( "mask" ) = nullptr )

        .def( "getTransientThermalForces", &DiscreteComputation::getTransientThermalForces, R"(
            Compute Transient Thermal Load

            Arguments:
                  time_curr (float): Current time
                  time_step (float): Time increment
                  theta (float): Theta parameter for integration
                  varc_curr (FieldOnCellsReal): external state variables at current time
                  previousPrimalField (fieldOnNodesReal): solution field at previous time

            Returns:
                  FieldOnNodes: load
            )",
              py::arg( "time_curr" ), py::arg( "time_step" ), py::arg( "theta" ),
              py::arg( "varc_curr" ) = nullptr, py::arg( "previousPrimalField" ) = nullptr )

        .def( "getTransientThermalForces", &DiscreteComputation::getTransientThermalForces, R"(
            Compute Transient Thermal forces due to time scheme
            Option CHAR_THER_EVOL

            Arguments:
                  time_curr (float): Current time
                  time_step (float): Time increment
                  theta (float): Theta parameter for integration
                  previousPrimalField (fieldOnNodesReal): solution field at previous time
                  varc_curr (FieldOnCellsReal): external state variables at current time

            Returns:
                  FieldOnNodes: load
            )",
              py::arg( "time_curr" ), py::arg( "time_step" ), py::arg( "theta" ),
              py::arg( "previousPrimalField" ), py::arg( "varc_curr" ) = nullptr )

        .def( "getTransientThermalLoadForces", &DiscreteComputation::getTransientThermalLoadForces,
              R"(
            Compute Transient Thermal Load given by EVOL_CHAR.
            Option CHAR_THER.

            Arguments:
                  time_curr (float): Current time
                  temp_prev (FieldOnNodesReal): solution field at previous time
                  assembly (bool) : if True return assembled vector (default: True)

            Returns:
                  FieldOnNodes: load
            )",
              py::arg( "time_curr" ), py::arg( "temp_prev" ) = nullptr,
              py::arg( "assembly" ) = true )

        .def( "getMechanicalDirichletBC", &DiscreteComputation::getMechanicalDirichletBC,
              R"(
            Return the imposed displacement vector used to remove imposed DDL
            *for internal use - prefer to use getDirichletBC*

            Arguments:
                  time_curr (float): Current time (default 0.0)

            Returns:
                  FieldOnNodesReal: imposed displacement vector
        )",
              py::arg( "time_curr" ) = 0.0 )

        .def( "getThermalDirichletBC", &DiscreteComputation::getThermalDirichletBC,
              R"(
            Return the imposed thermal vector used to remove imposed DDL
            *for internal use - prefer to use getDirichletBC*

            Arguments:
                  time_curr (float): Current time (default 0.0)

            Returns:
                  FieldOnNodesReal: imposed thermal vector
        )",
              py::arg( "time_curr" ) = 0.0 )

        .def( "getAcousticDirichletBC", &DiscreteComputation::getAcousticDirichletBC,
              R"(
            Return the imposed acoustic vector used to remove imposed DDL
            *for internal use - prefer to use getDirichletBC*

            Arguments:
                  time_curr (float): Current time (default 0.0)

            Returns:
                  FieldOnNodesComplex: imposed accoustic vector
        )",
              py::arg( "time_curr" ) = 0.0 )

        .def( "getIncrementalDirichletBC", &DiscreteComputation::getIncrementalDirichletBC,
              R"(
            Return the incremental imposed displacement vector used to remove imposed DDL
            for incremental resolution.

            incr_disp = getDirichletBC(time_curr) - disp, with 0.0 for DDL not imposed

            Arguments:
                  time_curr (float): Current time
                  disp (FieldOnNodes): displacement field at current time

            Returns:
                  FieldOnNodes: incremental imposed displacement vector
        )",
              py::arg( "time_curr" ), py::arg( "disp" ) )

        .def( "getElasticStiffnessMatrix", &DiscreteComputation::getElasticStiffnessMatrix, R"(
            Return the elementary matrices for elastic Stiffness matrix.
            Option RIGI_MECA.

            Arguments:
                  time_curr (float): Current time for external state variable
                    evaluation (default: 0.0)
                  fourierMode (int): Fourier mode (default: -1)
                  varc_curr (FieldOnCellsReal): external state variables at current time
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
                  with_dual (bool): compute dual terms or not (default: True)
            Returns:
                  ElementaryMatrix: elementary elastic Stiffness matrix
            )",
              py::arg( "time_curr" ) = 0.0, py::arg( "fourierMode" ) = -1,
              py::arg( "varc_curr" ) = nullptr, py::arg( "groupOfCells" ) = VectorString(),
              py::arg( "with_dual" ) = true )

        .def( "getFluidStructureStiffnessMatrix",
              &DiscreteComputation::getFluidStructureStiffnessMatrix,
              R"(
            Return the elementary matrices for fluid-structure stiffness matrix.
            Option RIGI_FLUI_STRUC.

            Arguments:
                  fourierMode (int): Fourier mode (default: -1)
                  varc_curr (FieldOnCells): internal state variables at current time
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
            Returns:
                  ElementaryMatrixReal: elementary fluid-structure Stiffness matrix
            )",
              py::arg( "fourierMode" ) = -1, py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getFluidStructureMassMatrix", &DiscreteComputation::getFluidStructureMassMatrix,
              R"(
            Return the elementary matrices for fluid-structure mass matrix.
            Option MASS_FLUI_STRUC.

            Arguments:
                  varc_curr (FieldOnCellsReal): external state variables at current time
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
            Returns:
                  ElementaryMatrixReal: elementary fluid-structure mass matrix
            )",
              py::arg( "varc_curr" ) = nullptr, py::arg( "groupOfCells" ) = VectorString() )

        .def( "getPhysicalProblem", &DiscreteComputation::getPhysicalProblem, R"(
            Get physical probelm

            Returns:
                  PhysicalProblem: physical problem
            )" )

        .def( "getDualElasticStiffnessMatrix", &DiscreteComputation::getDualElasticStiffnessMatrix,
              R"(
            Return elementary matrices for dual mechanical BC

            Returns:
                ElementaryMatrix: elementary matrices
        )" )

        .def( "getDualLinearMobilityMatrix", &DiscreteComputation::getDualLinearMobilityMatrix,
              R"(
            Return elementary matrices for dual acoustic BC

            Returns:
                ElementaryMatrix: elementary matrices
        )" )

        .def( "getDualLinearConductivityMatrix",
              &DiscreteComputation::getDualLinearConductivityMatrix,
              R"(
            Return elementary matrices for dual thermal BC

            Returns:
                ElementaryMatrix: elementary matrices
        )" )

        .def( "getLinearConductivityMatrix", &DiscreteComputation::getLinearConductivityMatrix,
              R"(
            Return the elementary matices for linear thermal matrix.
            Option RIGI_THER.

            Arguments:
                  time_curr (float): Current time
                  fourierMode (int): Fourier mode (default: -1)
                  varc_curr (FieldOnCellsReal): external state variables at current time
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                  with_dual (bool): compute dual terms or not (default: True)
            Returns:
                  ElementaryMatrix: elementary linear thermal matrices
        )",
              py::arg( "time_curr" ), py::arg( "fourierMode" ) = 0,
              py::arg( "varc_curr" ) = nullptr, py::arg( "groupOfCells" ) = VectorString(),
              py::arg( "with_dual" ) = true )

        .def( "getTangentConductivityMatrix", &DiscreteComputation::getTangentConductivityMatrix,
              R"(
            Return the elementary matrices for tangent conductivity.
            Option MASS_THER_TANG.

            Arguments:
                temp_prev (FieldOnNodes): thermal field at begin of current time
                temp_step (FieldOnNodes): field of increment of temperature
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                with_dual (bool): compute dual terms or not (default: True)
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "temp_prev" ), py::arg( "temp_step" ), py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "with_dual" ) = true )

        .def( "getThermalNonLinearVolumetricForces",
              &DiscreteComputation::getThermalNonLinearVolumetricForces,
              R"(
            Return the elementary field for nonlinear volumetric forces.
            Option CHAR_THER_SOURNL.

            Arguments:
                temp_curr (FieldOnNodesReal): thermal field at end of current time
                time_curr (float): Current time
                assembly (bool): assemble or not the field
            Returns:
                ElementaryVector: elementary field
            )",
              py::arg( "temp_curr" ), py::arg( "time_curr" ), py::arg( "assembly" ) = true )

        .def( "getThermalNonLinearNeumannForces",
              &DiscreteComputation::getThermalNonLinearNeumannForces,
              R"(
            Return the elementary field for nonlinear neuamnn forces.
            Option CHAR_THER_FLUTNL, CHAR_THER_RAYO_F, CHAR_THER_RAYO_R.

            Arguments:
                temp_curr (FieldOnNodesReal): thermal field at end of current time
                time_curr (float): Current time
                assembly (bool): assemble or not the field
            Returns:
                ElementaryVector: elementary field
            )",
              py::arg( "temp_curr" ), py::arg( "time_curr" ), py::arg( "assembly" ) = true )

        .def( "getThermalTangentNonLinearVolumetricMatrix",
              &DiscreteComputation::getThermalTangentNonLinearVolumetricMatrix,
              R"(
            Return the elementary matrices for tangent nonlinear volumetric forces.
            Option MTAN_THER_SOURNL.

            Arguments:
                temp_curr (FieldOnNodesReal): thermal field at end of current time
                time_curr (float): Current time

            Returns:
                ElementaryMatrix: elementary matrix
            )",
              py::arg( "temp_curr" ), py::arg( "time_curr" ) )

        .def( "getThermalTangentNonLinearNeumannMatrix",
              &DiscreteComputation::getThermalTangentNonLinearNeumannMatrix,
              R"(
            Return the elementary matrices for tangent nonlinear neumann forces.
            Option MTAN_THER_FLUXNL, MTAN_THER_RAYO_R, MTAN_THER_RAYO_F.

            Arguments:
                temp_curr (FieldOnNodesReal): thermal field at end of current time
                time_curr (float): Current time
                varc_curr (FieldOnCellsReal): external state variables at current time

            Returns:
                ElementaryMatrix: elementary matrix
            )",
              py::arg( "temp_curr" ), py::arg( "time_curr" ), py::arg( "varc_curr" ) = nullptr )

        .def( "getThermalExchangeMatrix", &DiscreteComputation::getThermalExchangeMatrix,
              R"(
            Return the elementary matices for exhange thermal matrix.

            Arguments:
                time_curr (float): Current time
            Returns:
                ElementaryMatrix: elementary exchange thermal matrices
        )",
              py::arg( "time_curr" ) )

        .def( "getLinearMobilityMatrix", &DiscreteComputation::getLinearMobilityMatrix,
              R"(
            Return the elementary matices for linear mobility acoustic matrix
            Option RIGI_ACOU.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                with_dual (bool): compute dual terms or not (default: True)

            Returns:
                ElementaryMatrix: elementary linear acoustic matrices
        )",
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "with_dual" ) = true )

        .def( "getMechanicalMassMatrix", &DiscreteComputation::getMechanicalMassMatrix, R"(
            Return the elementary matrices for mechanical mass matrix
            Option MASS_MECA.

            Arguments:
                diagonal (bool) : True for diagonal mass matrix else False.
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "diagonal" ), py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getCompressibilityMatrix", &DiscreteComputation::getCompressibilityMatrix, R"(
            Return the elementary matrices for compressibility acoustic matrix.
            Option MASS_ACOU.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getLinearCapacityMatrix", &DiscreteComputation::getLinearCapacityMatrix, R"(
            Return the elementary matrices for linear Capacity matrix in thermal computation.
            Option MASS_THER.

            Arguments:
                time_curr (float): current time to evaluate rho_cp
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "time_curr" ), py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getTangentCapacityMatrix", &DiscreteComputation::getTangentCapacityMatrix, R"(
            Return the elementary matrices for nonlinear Capacity matrix in thermal computation.
            Option MASS_THER_TANG.

            Arguments:
                temp_prev (FieldOnNodes): thermal field at begin of current time
                temp_step (FieldOnNodes): field of increment of temperature
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "temp_prev" ), py::arg( "temp_step" ), py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getMechanicalDampingMatrix", &DiscreteComputation::getMechanicalDampingMatrix, R"(
            Return the elementary matrices for damping matrix.
            Option AMOR_MECA.

            Arguments:
                getMechanicalMassMatrix : elementary mass matrix
                stiffnessMatrix : elementary stiffness matrix
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                flui_int (int): integer to activate damping impedance fluid matrix
                onde_flui (int): integer to indicate if we have an outgoing or incoming wave
            Returns:
                ElementaryMatrixReal: elementary damping matrix
            )",
              py::arg( "getMechanicalMassMatrix" ) = nullptr,
              py::arg( "stiffnessMatrix" ) = nullptr, py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "flui_int" ) = 1,
              py::arg( "onde_flui" ) = 1 )

        .def( "getImpedanceMatrix", &DiscreteComputation::getImpedanceMatrix, R"(
            Return the elementary matrices for impedance (acoustic) damping matrix.
            Option AMOR_ACOU.

            Returns:
                ElementaryMatrixReal: elementary damping matrix
            )" )

        .def( "getImpedanceBoundaryMatrix", &DiscreteComputation::getImpedanceBoundaryMatrix, R"(
            Return the elementary matrices for impedance (mechanical) matrix.
            Option IMPE_MECA.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                onde_flui (int): integer to indicate if we have an outgoing or incoming wave
            Returns:
                ElementaryMatrixReal: impedance mechanical matrix
            )",
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "onde_flui" ) = 1 )

        .def( "getImpedanceWaveMatrix", &DiscreteComputation::getImpedanceWaveMatrix, R"(
            Return the elementary matrices for impedance (mechanical) matrix
            from an harmonic wave.
            Option ONDE_FLUI.

            Returns:
                ElementaryMatrixReal: impedance wave matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getHystereticStiffnessMatrix", &DiscreteComputation::getHystereticStiffnessMatrix,
              R"(
            Return the elementary matrices for viscoelastic Stiffness matrix.
            Option RIGI_MECA_HYST.

            Arguments:
                stiffnessMatrix : elementary stiffness matrix
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixComplex: elementary viscoelastic rigidity matrix
            )",
              py::arg( "stiffnessMatrix" ), py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getGeometricStiffnessMatrix", &DiscreteComputation::getGeometricStiffnessMatrix, R"(
            Return the elementary matrices for geometric Stiffness matrix.
            Option RIGI_MECA_HYST.

            Arguments:
                sief_elga (FieldOnCellsReal) : stress at Gauss points
                strx_elga (FieldOnCellsReal) : stress at Gauss points for structural element
                displ (FieldOnNodesReal) : displacement field
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixComplex: elementary geometric rigidity matrix
            )",
              py::arg( "sief_elga" ), py::arg( "strx_elga" ) = nullptr,
              py::arg( "displ" ) = nullptr, py::arg( "modeFourier" ) = -1,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getRotationalStiffnessMatrix", &DiscreteComputation::getRotationalStiffnessMatrix,
              R"(
            Return the elementary matrices for rotational Stiffness matrix.
            Option RIGI_ROTA.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixReal: elementary rotational rigidity matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getGyroscopicStiffnessMatrix", &DiscreteComputation::getGyroscopicStiffnessMatrix,
              R"(
            Return the elementary matrices for gyroscopic Stiffness matrix.
            Option RIGI_GYRO.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixReal: elementary gyroscopic rigidity matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getGyroscopicDampingMatrix", &DiscreteComputation::getGyroscopicDampingMatrix, R"(
            Return the elementary matrices for gyroscopic damping matrix.
            Option MECA_GYRO.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixReal: elementary gyroscopic damping matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getInternalMechanicalForces", &DiscreteComputation::getInternalMechanicalForces,
              R"(
            Compute internal forces (integration of behaviour)

            Arguments:
                displ_prev (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                stress (FieldOnCells): field of stress at begin of current time
                internVar (FieldOnCells): field of internal state variables at begin of current time
                internVarIter (FieldOnCells): field of internal state variables at begin of
                                              current newton iteration
                time_prev (float): time at begin of the step
                time_step (float): delta time between begin and end of the step
                varc_prev (FieldOnCells): external state variables at begin of current time
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.

            Returns:
                tuple (tuple): return code error (FieldOnCells),
                error code flag (int),
                internal state variables VARI_ELGA (FieldOnCells),
                Cauchy stress SIEF_ELGA (FieldOnCells),
                field of internal forces (FieldOnNodesReal),
        )",
              py::arg( "displ_prev" ), py::arg( "displ_step" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "internVarIter" ), py::arg( "time_prev" ),
              py::arg( "time_step" ), py::arg( "varc_prev" ) = nullptr,
              py::arg( "varc_curr" ) = nullptr, py::arg( "groupOfCells" ) = VectorString() )

        .def( "getInternalThermalForces", &DiscreteComputation::getInternalThermalForces,
              R"(
            Compute internal thermal forces (integration of behaviour)
            Option RAPH_THER.

            Arguments:
                temp_prev (FieldOnNodes): thermal field at begin of current time
                temp_step (FieldOnNodes): field of increment of temperature
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used

            Returns:
                tuple (tuple):
                error code flag (int),
                fluxes FLUX_ELGA (FieldOnCellsReal),
                internal forces (FieldOnNodesReal),
            )",
              py::arg( "temp_prev" ), py::arg( "temp_step" ), py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getNonLinearCapacityForces", &DiscreteComputation::getNonLinearCapacityForces,
              R"(
            Compute internal thermal forces (integration of behaviour)
            Option MASS_THER_RESI.

            Arguments:
                temp_prev (FieldOnNodes): thermal field at begin of current time
                temp_step (FieldOnNodes): field of increment of temperature
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "temp_prev" ), py::arg( "temp_step" ), py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getTangentStiffnessMatrix", &DiscreteComputation::getTangentStiffnessMatrix,
              R"(
            Compute jacobian matrix for Newton algorithm

            Arguments:
                displ_prev (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                stress (FieldOnCells): field of stress at begin of current time
                internVar (FieldOnCells): internal state variables at begin of current time
                internVarIter (FieldOnCells): field of internal state variables
                                              at begin of current newton iteration
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                varc_prev (FieldOnCellsReal): external state variables at begin of current time
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.

            Returns:
                tuple (tuple): return code error (FieldOnCellsLong),
                error code flag (int),
                elementary tangent matrix (ElementaryMatrixDisplacementReal)
        )",
              py::arg( "displ_prev" ), py::arg( "displ_step" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "internVarIter" ), py::arg( "time_prev" ),
              py::arg( "time_step" ), py::arg( "varc_prev" ) = nullptr,
              py::arg( "varc_curr" ) = nullptr, py::arg( "groupOfCells" ) = VectorString() )

        .def( "getPredictionTangentStiffnessMatrix",
              &DiscreteComputation::getPredictionTangentStiffnessMatrix, R"(
            Compute jacobian matrix for Newton algorithm, Euler prediction

            Arguments:
                displ_prev (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                stress (FieldOnCells): field of stress at begin of current time
                internVar (FieldOnCells): internal state variables at begin of current time
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                varc_prev (FieldOnCellsReal): external state variables at begin of current time
                varc_curr (FieldOnCellsReal): external state variables at current time
                groupOfCells (list[str]): compute matrices on given groups of cells.

            Returns:
                tuple (tuple): return code error (FieldOnCellsLong),
                error code flag (int),
                elementary tangent matrix (ElementaryMatrixDisplacementReal)
        )",
              py::arg( "displ_prev" ), py::arg( "displ_step" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "time_prev" ), py::arg( "time_step" ),
              py::arg( "varc_prev" ) = nullptr, py::arg( "varc_curr" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getContactForces", &DiscreteComputation::getContactForces, R"(
            Compute contact and friction forces

            Arguments:
                geom (MeshCoordinatesField): coordinates of mesh used to compute normal
                displ_prev (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                data (FieldOnCellsReal): contact data
                coef_cont (FieldOnNodesReal) : contact coefficient
                coef_frot (FieldOnNodesReal) : friction coefficient

            Returns:
                FieldOnNodesReal: contact and friction forces

        )",
              py::arg( "geom" ), py::arg( "displ_prev" ), py::arg( "displ_step" ),
              py::arg( "time_prev" ), py::arg( "time_step" ), py::arg( "data" ),
              py::arg( "coef_cont" ), py::arg( "coef_frot" ) )

        .def( "getContactMatrix", &DiscreteComputation::getContactMatrix, R"(
            Compute contact matrix

            Arguments:
                geom (MeshCoordinatesField): coordinates of mesh used to compute normal
                displ_prev (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                data (FieldOnCellsReal): contact data
                coef_cont (FieldOnNodesReal) : contact coefficient
                coef_frot (FieldOnNodesReal) : friction coefficient

            Returns:
                ElementaryMatrixDisplacementReal: contact and friction elementary matrix
        )",
              py::arg( "geom" ), py::arg( "displ_prev" ), py::arg( "displ_step" ),
              py::arg( "time_prev" ), py::arg( "time_step" ), py::arg( "data" ),
              py::arg( "coef_cont" ), py::arg( "coef_frot" ) )

        .def( "getMechanicalCouplingForces", &DiscreteComputation::getMechanicalCouplingForces, R"(
            Compute coupling for LIAISON_MASSIF.

            Arguments:
                displ_prev (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step

            Returns:
                FieldOnNodesReal: coupling forces

        )",
              py::arg( "displ_prev" ), py::arg( "displ_step" ), py::arg( "time_prev" ),
              py::arg( "time_step" ) )

        .def( "getMechanicalCouplingMatrix", &DiscreteComputation::getMechanicalCouplingMatrix, R"(
            Compute coupling for LIAISON_MASSIF.

            Arguments:
                displ_prev (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step

            Returns:
                ElementaryMatrixDisplacementReal: coupling elementary matrix.

        )",
              py::arg( "displ_prev" ), py::arg( "displ_step" ), py::arg( "time_prev" ),
              py::arg( "time_step" ) )

        .def( "getMechanicalLinearCouplingMatrix",
              &DiscreteComputation::getMechanicalLinearCouplingMatrix, R"(
            Compute coupling for LIAISON_MASSIF.

      Arguments:
            varc_curr (FieldOnCellsReal): external state variables for Nitsche method.

      Returns:
            ElementaryMatrixDisplacementReal: coupling elementary matrix.

        )",
              py::arg( "varc_curr" ) = nullptr )

        .def( "getMechanicalNodalForces", &DiscreteComputation::getMechanicalNodalForces,
              R"(
      Return the elementary mechanical nodal forces vector

      Arguments:
            stress (FieldOnCells): field of stresses
            disp (FieldOnNodes): displacement field (required for large strains hypothesis)
            modeFourier (int) : fourier mode
            varc_curr (FieldOnCellsReal): external state variables
            behaviourMap (FieldOnCellsReal): map for non-linear behaviour
            groupOfCells (list[str]): compute vector on given groups of cells.
            assembly (bool) : if True return assembled vector (default: True)

      Returns:
            ElementaryVectorDisplacementReal: elementary Neumann forces vector
        )",
              py::arg( "stress" ), py::arg( "disp" ) = nullptr, py::arg( "modeFourier" ) = 0,
              py::arg( "varc_curr" ) = nullptr, py::arg( "behaviourMap" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "assembly" ) = true )

        .def( "getMechanicalForces", &DiscreteComputation::getMechanicalForces,
              R"(
      Return the total mechanical Neumann forces vector

      Arguments:
            time_curr (float): Current time
            time_step (float): Time increment
            theta (float): Theta parameter for time-integration
            modeFourier (int) : fourier mode
            varc_curr (FieldOnCellsReal): external state variables at current time

      Returns:
            FieldOnNodesReal: forces vector
        )",
              py::arg( "time_curr" ) = 0.0, py::arg( "time_step" ) = 0.0, py::arg( "theta" ) = 1.0,
              py::arg( "modeFourier" ) = 0, py::arg( "varc_curr" ) = nullptr )

        .def( "getMechanicalReactionForces", &DiscreteComputation::getMechanicalReactionForces,
              R"(
      Return the reaction forces

      Arguments:
            stress (FieldOnCells): field of stresses
            disp (FieldOnNodes): displacement field (required for large strains hypothesis)
            time_prev (float): time at begin of the step
            time_curr (float): time at end of the step
            theta (float): Theta parameter for time-integration
            modeFourier (int) : fourier mode
            varc_curr (FieldOnCellsReal): external state variables at current time
            behaviourMap (FieldOnCellsReal): map for non-linear behaviour

      Returns:
            FieldOnNodesReal: forces vector
        )",
              py::arg( "disp" ), py::arg( "stress" ), py::arg( "time_prev" ) = 0.0,
              py::arg( "time_curr" ) = 0.0, py::arg( "theta" ) = 1.0, py::arg( "modeFourier" ) = 0,
              py::arg( "varc_curr" ) = nullptr, py::arg( "behaviourMap" ) = nullptr )
        .def( "getResidualReference", &DiscreteComputation::getResidualReference,
              R"(
      Return the residual reference (for RESI_REFE_RELA)

      Arguments:
            vale_by_name : dict :
                           keys are component names
                           values are the given reference value corresponding to component name

      Returns:
            FieldOnNodesReal: residual reference forces vector
        )" );
};
