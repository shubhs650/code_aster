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


def coupled_fluid(cpl):
    """Run Fluid coupling.

    Arguments:
        cpl (ExternalCoupling): Fluid coupler
    """
    ################################################################################
    # setup the simulation
    ################################################################################
    # send signal 6 (abort) to produce a traceback

    test = CA.TestCase()

    CA.init("--test", comm=cpl.MPI.ASTER_COMM_WORLD, debug=False, ERREUR=_F(ERREUR_F="ABORT"))

    # Read the  mesh - 2 cases : 1 or several procs
    if CA.MPI.ASTER_COMM_WORLD.size > 1:
        MAFLUIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=19, PARTITIONNEUR="PTSCOTCH")
    else:
        MAFLUIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=19)

    MAFLUIDE = MODI_MAILLAGE(
        reuse=MAFLUIDE,
        MAILLAGE=MAFLUIDE,
        ORIE_PEAU=_F(GROUP_MA_PEAU=("Face1", "Face2", "Face3", "Face4", "Face5", "Face6")),
    )

    # Assign model
    MOFLUIDE = AFFE_MODELE(
        MAILLAGE=MAFLUIDE,
        AFFE=_F(
            GROUP_MA=("Face1", "Face2", "Face3", "Face4", "Face5", "Face6"),
            PHENOMENE="MECANIQUE",
            MODELISATION="3D",
        ),
    )

    FORM = FORMULE(VALE="1.E-4*INST*(X+Y+Z)", NOM_PARA=["X", "Y", "Z", "INST"])

    PRES_F = CREA_CHAMP(
        TYPE_CHAM="NOEU_NEUT_F",
        OPERATION="AFFE",
        MODELE=MOFLUIDE,
        AFFE=_F(TOUT="OUI", NOM_CMP="X1", VALE_F=FORM),
        INFO=1,
    )

    L_INST0 = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=5))

    ################################################################################
    # define one iteration
    ################################################################################

    class FluidSolver:
        """Class that allows to compute one iteration of the coupling.
           This class has to contains at least run() function.

        Attributes:
            result (LoadResult): result of the thermal computation.
        """

        def __init__(self, cpl):
            """cpl (ExternalCoupling): coupler."""

            self._medcpl = cpl.medcpl
            self.result = []

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

            CHINST = CA.FieldOnNodesReal(MAFLUIDE, "INST_R", {"INST": current_time})

            CHXN = CREA_CHAMP(
                OPERATION="EXTR", TYPE_CHAM="NOEU_GEOM_R", NOM_CHAM="GEOMETRIE", MAILLAGE=MAFLUIDE
            )

            CHNEUT = CREA_CHAMP(
                OPERATION="EVAL", TYPE_CHAM="NOEU_NEUT_R", CHAM_F=PRES_F, CHAM_PARA=(CHXN, CHINST)
            )

            pres = CHNEUT.asPhysicalQuantity("PRES_R", {"X1": "PRES"})

            self.result.append(pres)

            self._medcpl.export_pressure("PRES", pres)

            return True

    ################################################################################
    # loop on time steps
    ################################################################################

    fluid_solv = FluidSolver(cpl)

    cpl.setup(
        interface=(MAFLUIDE, ["Face2", "Face3", "Face4", "Face5", "Face6"]),
        input_fields=[("DEPL", ["DX", "DY", "DZ"], "NODES")],
        output_fields=[("PRES", ["PRES"], "NODES")],
    )

    cpl.run(fluid_solv, time_list=L_INST0.getValues())

    # Extract the field from the result

    test.assertEqual(len(fluid_solv.result), 5)

    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    test.printSummary()

    CA.close()
