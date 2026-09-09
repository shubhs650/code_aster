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


def coupled_thermics(cpl):
    """Run thermal coupling.

    Arguments:
        cpl (ExternalCoupling): Thermal coupler
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

    MQ = MODI_MAILLAGE(reuse=MQ, MAILLAGE=MQ, ORIE_PEAU=_F(GROUP_MA_PEAU="L2"))

    # Assign thermal model on linear model
    MODE_TL = AFFE_MODELE(
        MAILLAGE=ML, AFFE=_F(TOUT="OUI", MODELISATION="3D", PHENOMENE="THERMIQUE")
    )

    # Assign thermal loading
    CLIMT = AFFE_CHAR_THER(MODELE=MODE_TL, FLUX_REP=(_F(GROUP_MA="L1", FLUN=-400),))

    # Assign thermal Dirichlet BC
    BLOQT = AFFE_CHAR_CINE(MODELE=MODE_TL, THER_IMPO=_F(GROUP_MA="L2", TEMP=10))

    # Define the thermal material
    ACIER_T = DEFI_MATERIAU(THER=_F(LAMBDA=33.5, RHO_CP=526.0e4))

    # Assign the thermal material
    MATE_T = AFFE_MATERIAU(MAILLAGE=ML, AFFE=_F(TOUT="OUI", MATER=ACIER_T))

    L_INST = DEFI_LIST_REEL(VALE=(-1.0, 0.0))

    ################################################################################
    # define one iteration
    ################################################################################

    class ThermalSolver:
        """Class that allows to compute one iteration of the coupling.
           This class has to contains at least run() function.

        Attributes:
            result (ThermalResult): result of the thermal computation.
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

            mc_depl = data["DEPL"]
            if mc_depl:
                # MEDC field => .med => code_aster field
                depl = self._medcpl.import_displacement("DEPL")

            self.result = THER_LINEAIRE(
                MODELE=MODE_TL,
                CHAM_MATER=MATE_T,
                EXCIT=(_F(CHARGE=CLIMT), _F(CHARGE=BLOQT)),
                # SOLVEUR=_F(METHODE="PETSC", PRE_COND="HPDDM", RESI_RELA=1.0e-9),
                SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_DP", RESI_RELA=1.0e-9),
                TYPE_CALCUL="STAT",
                INCREMENT=_F(LIST_INST=DEFI_LIST_REEL(VALE=current_time)),
                INFO=1,
            )

            temp = self.result.getField("TEMP", self.result.getLastIndex())
            self._medcpl.export_temperature("TEMP", temp)
            print("[Convert] Temperature field info:")
            print(self._medcpl.get_field("TEMP").getField().simpleRepr(), flush=True)

            return True

    ################################################################################
    # loop on time steps
    ################################################################################

    ther_solv = ThermalSolver(cpl)

    cpl.setup(
        interface=(ML, ["Volume"]),
        input_fields=[("DEPL", ["DX", "DY", "DZ"], "NODES")],
        output_fields=[("TEMP", ["TEMP"], "NODES")],
    )

    cpl.run(ther_solv, time_list=L_INST.getValues())

    # Extract the field from the result
    temp = ther_solv.result.getField("TEMP", 1)
    test.assertAlmostEqual(temp.norm("NORM_2"), 218.6233131691481)

    # Assign thermal model on linear model
    MODE_TQ = AFFE_MODELE(
        MAILLAGE=MQ, AFFE=_F(TOUT="OUI", MODELISATION="3D", PHENOMENE="THERMIQUE")
    )

    # Project the field from the linear to the quadratic thermal model
    TEMP2 = PROJ_CHAMP(
        METHODE="COLLOCATION", RESULTAT=ther_solv.result, MODELE_1=MODE_TL, MODELE_2=MODE_TQ
    )

    # Validate the result
    temp2 = TEMP2.getField("TEMP", 1)
    test.assertAlmostEqual(temp2.norm("NORM_INFINITY"), temp.norm("NORM_INFINITY"))
    test.assertAlmostEqual(temp2.norm("NORM_INFINITY"), 21.7467518133881)

    test.assertAlmostEqual(temp2.norm("NORM_2"), 535.323717828298)

    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    test.printSummary()

    CA.close()
