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


def coupled_fluid(cpl, UNITE_MA):
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
        MAFLUIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=UNITE_MA, PARTITIONNEUR="PTSCOTCH")
    else:
        MAFLUIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=UNITE_MA)

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
    FORM0 = FORMULE(VALE="0.", NOM_PARA=["X", "Y", "Z", "INST"])

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
            self.epsilon = 1e-7
            self.depl_prev = None
            self.result = []
            self._use_CFEMDEC = cpl._use_CFEMDEC

        def extent(self, depl):
            fns = depl.toSimpleFieldOnNodes()

            val, desc = fns.getValuesWithDescription()

            fns.setValues([0] * 3 * MAFLUIDE.getNumberOfNodes())
            fns.setValues(desc[0], desc[1], val)

            return fns.toFieldOnNodes()

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

            assert len(data) == 2, "expecting one field"

            mc_depl = data["mesh_displacement"]
            depl = None
            if mc_depl:
                # MEDC field => .med => code_aster field
                depl = self._medcpl.import_displacement("mesh_displacement")

            CHINST = CA.FieldOnNodesReal(MAFLUIDE, "INST_R", {"INST": current_time})

            CHXN = CREA_CHAMP(
                OPERATION="EXTR", TYPE_CHAM="NOEU_GEOM_R", NOM_CHAM="GEOMETRIE", MAILLAGE=MAFLUIDE
            )

            if depl:
                DEPL = self.extent(depl)
                CHXN += DEPL

            PRES_F = CREA_CHAMP(
                TYPE_CHAM="NOEU_NEUT_F",
                OPERATION="AFFE",
                MAILLAGE=MAFLUIDE,
                AFFE=_F(TOUT="OUI", NOM_CMP=("X1", "X2", "X3"), VALE_F=(FORM, FORM0, FORM0)),
                INFO=1,
            )

            CHNEUT = CREA_CHAMP(
                OPERATION="EVAL", TYPE_CHAM="NOEU_NEUT_R", CHAM_F=PRES_F, CHAM_PARA=(CHXN, CHINST)
            )

            force_noeu = CHNEUT.asPhysicalQuantity("FORC_R", {"X1": "FX", "X2": "FY", "X3": "FZ"})

            force_elem = force_noeu.toFieldOnCells(MOFLUIDE.getFiniteElementDescriptor(), "ELEM")

            if i_iter == 0:
                self.result.append(force_elem)
            else:
                self.result[-1] = force_elem

            # export
            if self._use_CFEMDEC:
                self._medcpl.set_field("boundary_stress", force_noeu)
            else:
                self._medcpl.set_field("boundary_stress", force_elem)

            # test convergence:
            has_cvg = False
            if i_iter > 0:
                if self.depl_prev:
                    depl_incr = depl - self.depl_prev
                    resi_rela = depl_incr.norm("NORM_INFINITY") / depl.norm("NORM_INFINITY")
                    has_cvg = resi_rela < self.epsilon
                    print(f"RESI_CPL: #iter {i_iter}, #resi: {resi_rela}")

            self.depl_prev = depl

            return has_cvg

    ################################################################################
    # loop on time steps
    ################################################################################

    fluid_solv = FluidSolver(cpl)

    cpl.setup(
        interface=(MAFLUIDE, ["Face2", "Face3", "Face4", "Face5", "Face6"]),
        init_time=0.0,
        final_time=1.0,
        nb_step=5,
        nb_iter=10,
        epsilon=fluid_solv.epsilon,
    )

    cpl.run(fluid_solv)

    # Extract the field from the result

    test.assertEqual(len(fluid_solv.result), 5)

    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    test.printSummary()

    CA.close()
