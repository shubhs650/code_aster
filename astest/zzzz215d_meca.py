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

    # Read the mesh - 2 cases : 1 or several procs
    if CA.MPI.ASTER_COMM_WORLD.size > 1:
        MASOLIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=20, PARTITIONNEUR="PTSCOTCH")
    else:
        MASOLIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

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

    ACIER = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=0.3))

    MATER = AFFE_MATERIAU(MAILLAGE=MASOLIDE, AFFE=_F(TOUT="OUI", MATER=ACIER))

    CHA_IMPO = AFFE_CHAR_MECA(
        MODELE=MOSOLIDE, FACE_IMPO=_F(GROUP_MA="Face1", DX=0.0, DY=0.0, DZ=0.0)
    )

    RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1.0, 1.0))

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
            PRES = self._medcpl.import_pressure("PRES")

            fed = MOSOLIDE.getFiniteElementDescriptor().restrict(
                self._medcpl.mesh_interf.getGroupsOfCells()
            )
            PRES_elno = PRES.toFieldOnCells(fed, "ELNO")

            RES_PROJ = CREA_RESU(
                OPERATION="AFFE",
                TYPE_RESU="EVOL_CHAR",
                AFFE=_F(NOM_CHAM="PRES", CHAM_GD=PRES_elno, INST=current_time),
            )

            CHA_PROJ = AFFE_CHAR_MECA(MODELE=MOSOLIDE, EVOL_CHAR=RES_PROJ)

            opts = {}
            if self.result:
                opts["reuse"] = self.result
                opts["RESULTAT"] = self.result
                opts["ETAT_INIT"] = _F(EVOL_NOLI=self.result)

            previous_time = current_time - delta_t

            # Solve the mechanical problem
            self.result = STAT_NON_LINE(
                MODELE=MOSOLIDE,
                CHAM_MATER=MATER,
                EXCIT=(_F(CHARGE=CHA_IMPO), _F(CHARGE=CHA_PROJ)),
                COMPORTEMENT=_F(RELATION="ELAS"),
                INCREMENT=_F(
                    LIST_INST=DEFI_LIST_REEL(VALE=(previous_time, current_time)),
                    INST_INIT=previous_time,
                ),
                **opts,
            )

            displ = self.result.getField("DEPL", self.result.getLastIndex())
            self._medcpl.export_displacement("DEPL", displ)

            return True

    ################################################################################
    # loop on time steps
    ################################################################################

    mech_solv = MechanicalSolver(cpl)

    cpl.setup(
        interface=(MASOLIDE, ["Face2", "Face3", "Face4", "Face5", "Face6"]),
        input_fields=[("PRES", ["PRES"], "NODES")],
        output_fields=[("DEPL", ["DX", "DY", "DZ"], "NODES")],
    )

    cpl.run(mech_solv)

    RESU = mech_solv.result

    TEST_RESU(
        RESU=(
            _F(
                INST=1.0,
                RESULTAT=RESU,
                NOM_CHAM="DEPL",
                GROUP_NO="N134",
                NOM_CMP="DX",
                VALE_CALC=-9.323277465885031,
            ),
            _F(
                INST=1.0,
                RESULTAT=RESU,
                NOM_CHAM="DEPL",
                GROUP_NO="N134",
                NOM_CMP="DY",
                VALE_CALC=-2.3117208947423626,
            ),
            _F(
                INST=1.0,
                RESULTAT=RESU,
                NOM_CHAM="DEPL",
                GROUP_NO="N134",
                NOM_CMP="DZ",
                VALE_CALC=-6.949019476525089,
            ),
        )
    )
    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    FIN()
