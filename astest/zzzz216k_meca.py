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


def coupled_mechanics(cpl, UNITE_MA, test_vale):
    """Run mechanical coupling.

    Arguments:
        cpl (ExternalCoupling): Mechanical coupler
    """

    ################################################################################
    # setup the simulation
    ################################################################################
    # send signal 6 (abort) to produce a traceback

    CA.init("--test", comm=cpl.MPI.ASTER_COMM_WORLD, debug=False, ERREUR=_F(ERREUR_F="ABORT"))

    # Read the mesh - 2 cases : 1 or several procs
    if CA.MPI.ASTER_COMM_WORLD.size > 1:
        MASOLIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=UNITE_MA, PARTITIONNEUR="PTSCOTCH")
    else:
        MASOLIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=UNITE_MA)

    # Check the orientation of the boundary
    MASOLIDE = MODI_MAILLAGE(
        reuse=MASOLIDE,
        MAILLAGE=MASOLIDE,
        ORIE_PEAU=_F(GROUP_MA_PEAU=("Face1", "Face2", "Face3", "Face4", "Face5", "Face6")),
    )

    # Assign Mechanical model
    MOSOLIDE = AFFE_MODELE(
        MAILLAGE=MASOLIDE, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D")
    )

    ACIER = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=0.3, RHO=1.0))

    MATER = AFFE_MATERIAU(MAILLAGE=MASOLIDE, AFFE=_F(TOUT="OUI", MATER=ACIER))

    CHA_IMPO = AFFE_CHAR_CINE(
        MODELE=MOSOLIDE, MECA_IMPO=_F(GROUP_MA="Face1", DX=0.0, DY=0.0, DZ=0.0)
    )

    init_zero = CREA_CHAMP(
        TYPE_CHAM="NOEU_DEPL_R",
        OPERATION="AFFE",
        MODELE=MOSOLIDE,
        AFFE=_F(TOUT="OUI", NOM_CMP=("DX", "DY", "DZ"), VALE=(0.0, 0.0, 0.0)),
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

            self.cpl = cpl
            self.result = None

        def run_iteration(self, i_iter, current_time, delta_t):
            """Execute one iteration.

            Arguments:
                i_iter (int): Iteration number if the current time_step.
                current_time (float): Current time.
                delta_t (float): Time step.
            """

            FORCE = self.cpl.import_fluidLoad(MOSOLIDE, current_time)

            FORC = FORCE.getField("FSUR_3D", current_time, "INST")

            PRES = FORC.asPhysicalQuantity("PRES_R", {"FX": "PRES"})

            evol_char = CREA_RESU(
                TYPE_RESU="EVOL_CHAR",
                OPERATION="AFFE",
                AFFE=_F(NOM_CHAM="PRES", CHAM_GD=PRES, INST=current_time),
            )

            CHA_PROJ = AFFE_CHAR_MECA(MODELE=MOSOLIDE, EVOL_CHAR=evol_char)

            previous_time = current_time - delta_t

            opts = {}
            if self.result:
                opts["reuse"] = self.result
                opts["RESULTAT"] = self.result
                opts["ETAT_INIT"] = _F(EVOL_NOLI=self.result)
            else:
                opts["ETAT_INIT"] = _F(DEPL=init_zero, VITE=init_zero, ACCE=init_zero)

            # Solve the mechanical problem
            self.result = DYNA_NON_LINE(
                MODELE=MOSOLIDE,
                CHAM_MATER=MATER,
                EXCIT=(_F(CHARGE=CHA_IMPO), _F(CHARGE=CHA_PROJ, TYPE_CHARGE="SUIV")),
                COMPORTEMENT=_F(RELATION="ELAS"),
                INCREMENT=_F(
                    LIST_INST=DEFI_LIST_REEL(VALE=(previous_time, current_time)),
                    INST_INIT=previous_time,
                ),
                SCHEMA_TEMPS=_F(FORMULATION="DEPLACEMENT", SCHEMA="NEWMARK"),
                **opts,
            )

            displ = self.result.getField("DEPL", self.result.getLastIndex())
            self.cpl.export_displacement(displ)

            velo = self.result.getField("VITE", self.result.getLastIndex())
            self.cpl.export_velocity(velo)

    ################################################################################
    # loop on time steps
    ################################################################################

    mech_solv = MechanicalSolver(cpl)

    cpl.setup(interface=(MASOLIDE, ["Face2", "Face3", "Face4", "Face5", "Face6"]))

    assert cpl.run(mech_solv)

    RESU = mech_solv.result

    TEST_RESU(
        RESU=(
            _F(
                INST=0.6,
                RESULTAT=RESU,
                NOM_CHAM="DEPL",
                GROUP_NO="N134",
                NOM_CMP="DX",
                VALE_CALC=test_vale[0],
                REFERENCE="AUTRE_ASTER",
                VALE_REFE=-2.146353385213317,
                PRECISION=0.01,
            ),
            _F(
                INST=1.0,
                RESULTAT=RESU,
                NOM_CHAM="DEPL",
                GROUP_NO="N134",
                NOM_CMP="DX",
                VALE_CALC=test_vale[1],
                REFERENCE="AUTRE_ASTER",
                VALE_REFE=-9.795579740520052,
                PRECISION=0.01,
            ),
            _F(
                INST=1.0,
                RESULTAT=RESU,
                NOM_CHAM="DEPL",
                GROUP_NO="N134",
                NOM_CMP="DY",
                VALE_CALC=test_vale[2],
                REFERENCE="AUTRE_ASTER",
                VALE_REFE=-0.6839627581316297,
                PRECISION=0.01,
            ),
            _F(
                INST=1.0,
                RESULTAT=RESU,
                NOM_CHAM="DEPL",
                GROUP_NO="N134",
                NOM_CMP="DZ",
                VALE_CALC=test_vale[3],
                REFERENCE="AUTRE_ASTER",
                VALE_REFE=-2.2849657446018616,
                PRECISION=0.025,
            ),
        )
    )
    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    FIN()
