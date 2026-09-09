# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------


from code_aster.Commands import *
from code_aster import CA


def coupled_mechanics(cpl):
    """Run mechanical coupling.

    Arguments:
        cpl (ExternalCoupling): Mechanical coupler
    """

    ################################################################################
    # setup the simulation
    ################################################################################
    # send signal 6 (abort) to produce a traceback

    test = CA.TestCase()

    CA.init("--test", comm=cpl.MPI.ASTER_COMM_WORLD, debug=False, ERREUR=_F(ERREUR_F="ABORT"))

    # Read the quadratic mesh - 2 cases : 1 or several procs
    if CA.MPI.ASTER_COMM_WORLD.size > 1:
        MQ = LIRE_MAILLAGE(FORMAT="MED", UNITE=20, PARTITIONNEUR="PTSCOTCH")
    else:
        MQ = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

    # Create the linear mesh
    ML = CREA_MAILLAGE(MAILLAGE=MQ, QUAD_LINE=_F(TOUT="OUI"), INFO=1)

    # Check the orientation of the boundary
    ML = MODI_MAILLAGE(reuse=ML, MAILLAGE=ML, ORIE_PEAU=_F(GROUP_MA_PEAU="L1"))

    # Check the orientation of the boundary
    MQ = MODI_MAILLAGE(reuse=MQ, MAILLAGE=MQ, ORIE_PEAU=_F(GROUP_MA_PEAU="L2"))

    # Assign Mechanical model
    MODE_MQ = AFFE_MODELE(
        MAILLAGE=ML, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D")
    )

    # Assign thermal model on linear model
    MODE_TL = AFFE_MODELE(
        MAILLAGE=ML, AFFE=_F(TOUT="OUI", MODELISATION="3D", PHENOMENE="THERMIQUE")
    )

    # Define the mechanical material
    ACIER_M = DEFI_MATERIAU(ELAS=_F(E=2.0e12, NU=0.3e00, RHO=1.0e03, ALPHA=1.0e-4))

    # Assign mechanical Dirichlet BC
    DIRI = AFFE_CHAR_CINE(MODELE=MODE_MQ, MECA_IMPO=(_F(GROUP_MA="L1", DX=0, DY=0.0, DZ=0.0),))

    # Assign mechanical loading
    CHAR = AFFE_CHAR_MECA(MODELE=MODE_MQ, PRES_REP=_F(GROUP_MA="L2", PRES=-100.0))

    ################################################################################
    # define one iteration
    ################################################################################

    class MechanicalSolver:
        """Class that allows to compute one iteration of the coupling.
           This class has to contains at least run() function.

        Attributes:
            result (NonLinearResult): result of the mechanical computation.
            evol_ther (ThermalResult) : evolution of the thermal field.
        """

        def __init__(self, cpl):
            """cpl (ExternalCoupling): coupler."""

            self._medcpl = cpl.medcpl

        def run_iteration(self, i_iter, current_time, delta_t, data):
            """Execute one iteration.

            Arguments:
                i_iter (int): Iteration number if the current time_step.
                current_time (float): Current time.
                delta_t (float): Time step.
                data (dict[*MEDCouplingField*]): dict of input fields.

            Returns:
                bool: True if solver has converged at the current time step, else False.
            """

            assert len(data) == 1, "expecting one field"

            # MEDC field => .med => code_aster field
            TEMP = self._medcpl.import_temperature("TEMP")

            self.evol_ther = CREA_RESU(
                TYPE_RESU="EVOL_THER",
                OPERATION="AFFE",
                AFFE=_F(NOM_CHAM="TEMP", CHAM_GD=TEMP, MODELE=MODE_TL, INST=current_time),
            )

            # Assign the mechanical material with the thermal field
            MATE_M = AFFE_MATERIAU(
                MAILLAGE=ML,
                AFFE=_F(TOUT="OUI", MATER=ACIER_M),
                AFFE_VARC=_F(
                    TOUT="OUI", EVOL=self.evol_ther, VALE_REF=0.0, NOM_VARC="TEMP", NOM_CHAM="TEMP"
                ),
            )

            # Solve the mechanical problem
            self.result = MECA_STATIQUE(
                MODELE=MODE_MQ,
                CHAM_MATER=MATE_M,
                EXCIT=(_F(CHARGE=CHAR), _F(CHARGE=DIRI)),
                SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_SP"),
            )

            displ = self.result.getField("DEPL", self.result.getLastIndex())
            self._medcpl.export_displacement("DEPL", displ)
            print("[Convert] Displacement field info:")
            print(self._medcpl.get_field("DEPL").getField().simpleRepr(), flush=True)

            return True

    ################################################################################
    # loop on time steps
    ################################################################################

    mech_solv = MechanicalSolver(cpl)

    cpl.setup(
        interface=(ML, "Volume"),
        input_fields=[("TEMP", ["TEMP"], "NODES")],
        output_fields=[("DEPL", ["DX", "DY", "DZ"], "NODES")],
    )

    cpl.run(mech_solv)

    # Extract the field from the result
    displ = mech_solv.result.getField("DEPL", 1)

    # Validate the result againt sequential run
    # norm = displ.norm("NORM_2")
    norm = 0.05732070604120847
    test.assertAlmostEqual(displ.norm("NORM_2"), norm)

    # test interpolation of thermics.
    temp = mech_solv.evol_ther.getField("TEMP", 1)
    test.assertAlmostEqual(temp.norm("NORM_2"), 218.6233131691481)

    test.printSummary()

    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    CA.close()
