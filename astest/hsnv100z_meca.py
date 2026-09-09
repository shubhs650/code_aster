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
    CA.init("--test", comm=cpl.MPI.ASTER_COMM_WORLD, debug=False, ERREUR=_F(ERREUR_F="ABORT"))

    MAIL = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

    model = AFFE_MODELE(
        AFFE=_F(MODELISATION="AXIS", PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=MAIL
    )

    cpl.setup(
        interface=(MAIL, ["M1"]),
        input_fields=[("TEMP", ["TEMP"], "NODES")],
        output_fields=[("DEPL", ["DX", "DY"], "NODES")],
    )

    # DONNEES DE MODELISATION

    FCT1 = DEFI_FONCTION(
        NOM_PARA="EPSI",
        VALE=(0.200e-2, 400.0, 0.400e-2, 500.0),
        PROL_DROITE="LINEAIRE",
        PROL_GAUCHE="LINEAIRE",
    )

    #

    FCT2 = DEFI_FONCTION(
        NOM_PARA="EPSI",
        VALE=(0.100e-2, 200.0, 0.300e-2, 300.0),
        PROL_DROITE="LINEAIRE",
        PROL_GAUCHE="LINEAIRE",
    )

    #

    CTRACB = DEFI_NAPPE(
        NOM_PARA="TEMP",
        PARA=(0.0, 50.0),
        FONCTION=(FCT1, FCT2),
        PROL_DROITE="LINEAIRE",
        PROL_GAUCHE="LINEAIRE",
    )

    MAT = DEFI_MATERIAU(ELAS=_F(E=200.0e3, NU=0.3, ALPHA=10.0e-6), TRACTION=_F(SIGM=CTRACB))

    CHMECA = AFFE_CHAR_CINE(
        MODELE=model, MECA_IMPO=(_F(GROUP_NO="GRNO1", DY=0.0), _F(GROUP_NO="GRNO3", DY=0.0))
    )

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
            input_data = cpl.recv_input_fields()
            TEMPE = self._medcpl.import_temperature("TEMP")

            self.evol_ther = CREA_RESU(
                TYPE_RESU="EVOL_THER",
                OPERATION="AFFE",
                AFFE=_F(NOM_CHAM="TEMP", CHAM_GD=TEMPE, INST=0.0),
            )

            self.result = None

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
            TEMPE = self._medcpl.import_temperature("TEMP")

            self.evol_ther = CREA_RESU(
                reuse=self.evol_ther,
                RESULTAT=self.evol_ther,
                TYPE_RESU="EVOL_THER",
                OPERATION="AFFE",
                AFFE=_F(NOM_CHAM="TEMP", CHAM_GD=TEMPE, INST=current_time),
            )

            CTM = AFFE_MATERIAU(
                MAILLAGE=MAIL,
                AFFE=_F(TOUT="OUI", MATER=MAT),
                AFFE_VARC=_F(TOUT="OUI", NOM_VARC="TEMP", EVOL=self.evol_ther, VALE_REF=0.0),
            )

            opts = {}
            if self.result:
                opts["reuse"] = self.result
                opts["RESULTAT"] = self.result
                opts["ETAT_INIT"] = _F(EVOL_NOLI=self.result)

            previous_time = current_time - delta_t

            self.result = STAT_NON_LINE(
                MODELE=model,
                CHAM_MATER=CTM,
                COMPORTEMENT=_F(RELATION="VMIS_ISOT_TRAC"),
                EXCIT=_F(CHARGE=CHMECA),
                INCREMENT=_F(
                    LIST_INST=DEFI_LIST_REEL(VALE=(previous_time, current_time)),
                    INST_INIT=previous_time,
                ),
                **opts,
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

    cpl.run(mech_solv)

    mech_solv.result = CALC_CHAMP(
        reuse=mech_solv.result,
        RESULTAT=mech_solv.result,
        CONTRAINTE=("SIGM_ELNO",),
        DEFORMATION=("EPSI_ELNO",),
        VARI_INTERNE=("VARI_ELNO",),
    )

    TEST_RESU(
        RESU=(
            _F(
                INST=66.665999999999997,
                RESULTAT=mech_solv.result,
                NOM_CHAM="EPSI_ELNO",
                GROUP_NO="N1",
                NOM_CMP="EPXX",
                VALE_CALC=8.66658e-4,
                GROUP_MA="M1",
            ),
            _F(
                INST=66.665999999999997,
                RESULTAT=mech_solv.result,
                NOM_CHAM="SIGM_ELNO",
                GROUP_NO="N1",
                NOM_CMP="SIYY",
                VALE_CALC=-133.332,
                GROUP_MA="M1",
            ),
            _F(
                INST=80.0,
                RESULTAT=mech_solv.result,
                NOM_CHAM="EPSI_ELNO",
                GROUP_NO="N2",
                NOM_CMP="EPZZ",
                VALE_CALC=1.1000000000000001e-3,
                GROUP_MA="M1",
            ),
            _F(
                INST=80.0,
                RESULTAT=mech_solv.result,
                NOM_CHAM="VARI_ELNO",
                GROUP_NO="N2",
                NOM_CMP="V1",
                VALE_CALC=2.9999999999999997e-4,
                GROUP_MA="M1",
            ),
            _F(
                INST=80.0,
                RESULTAT=mech_solv.result,
                NOM_CHAM="SIGM_ELNO",
                GROUP_NO="N2",
                NOM_CMP="SIYY",
                VALE_CALC=-100.0,
                GROUP_MA="M1",
            ),
            _F(
                INST=90.0,
                RESULTAT=mech_solv.result,
                NOM_CHAM="EPSI_ELNO",
                GROUP_NO="N3",
                NOM_CMP="EPZZ",
                VALE_CALC=1.2750000000000001e-3,
                GROUP_MA="M1",
            ),
            _F(
                INST=90.0,
                RESULTAT=mech_solv.result,
                NOM_CHAM="VARI_ELNO",
                GROUP_NO="N3",
                NOM_CMP="V1",
                VALE_CALC=5.2499999999999997e-4,
                GROUP_MA="M1",
            ),
            _F(
                INST=90.0,
                RESULTAT=mech_solv.result,
                NOM_CHAM="SIGM_ELNO",
                GROUP_NO="N3",
                NOM_CMP="SIYY",
                VALE_CALC=-75.0,
                GROUP_MA="M1",
            ),
        )
    )

    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    FIN()
