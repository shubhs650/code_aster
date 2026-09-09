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
    CA.init("--test", comm=cpl.MPI.ASTER_COMM_WORLD, debug=False, ERREUR=_F(ERREUR_F="ABORT"))

    MAIL = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

    cpl.setup(
        interface=(MAIL, ["M1"]),
        input_fields=[("DEPL", ["DX", "DY"], "NODES")],
        output_fields=[("TEMP", ["TEMP"], "NODES")],
    )

    model = AFFE_MODELE(
        AFFE=_F(MODELISATION="AXIS", PHENOMENE="THERMIQUE", TOUT="OUI"), MAILLAGE=MAIL
    )

    MAT = DEFI_MATERIAU(THER=_F(LAMBDA=1.0e-3, RHO_CP=0.0e-3))

    CM = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MAT))

    TIMPVAR = DEFI_FONCTION(NOM_PARA="INST", NOM_RESU="TEMP", VALE=(0.0e0, 0.0e0, 100.0e0, 100.0e0))

    CHTHER = AFFE_CHAR_THER_F(
        MODELE=model,
        TEMP_IMPO=(
            _F(GROUP_NO="GRNO1", TEMP=TIMPVAR),
            _F(GROUP_NO="GRNO2", TEMP=TIMPVAR),
            _F(GROUP_NO="GRNO3", TEMP=TIMPVAR),
            _F(GROUP_NO="GRNO4", TEMP=TIMPVAR),
        ),
    )

    t1 = 66.666

    t3 = 80.0

    t5 = 90.0

    L_INST = DEFI_LIST_REEL(
        DEBUT=0.0,
        INTERVALLE=(_F(JUSQU_A=t1, NOMBRE=1), _F(JUSQU_A=t3, NOMBRE=2), _F(JUSQU_A=t5, NOMBRE=2)),
    )

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

            T0 = CA.FieldOnNodesReal(model)
            T0.setValues(0.0)

            self.result = CREA_RESU(
                TYPE_RESU="EVOL_THER",
                OPERATION="AFFE",
                AFFE=_F(NOM_CHAM="TEMP", CHAM_GD=T0, INST=0.0, MODELE=model, CHAM_MATER=CM),
            )

            self._medcpl = cpl.medcpl

            self._medcpl.export_temperature("TEMP", T0)
            cpl.send_output_fields()

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

            previous_time = current_time - delta_t

            self.result = THER_LINEAIRE(
                reuse=self.result,
                RESULTAT=self.result,
                MODELE=model,
                CHAM_MATER=CM,
                EXCIT=_F(CHARGE=CHTHER),
                INCREMENT=_F(
                    LIST_INST=DEFI_LIST_REEL(VALE=(previous_time, current_time)),
                    INST_INIT=previous_time,
                ),
                ETAT_INIT=_F(EVOL_THER=self.result),
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

    cpl.run(ther_solv, time_list=L_INST.getValues())

    TEST_RESU(
        RESU=(
            _F(
                INST=90.0,
                RESULTAT=ther_solv.result,
                NOM_CHAM="TEMP",
                GROUP_NO="N1",
                NOM_CMP="TEMP",
                VALE_CALC=90.0,
            ),
        )
    )

    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    FIN()
