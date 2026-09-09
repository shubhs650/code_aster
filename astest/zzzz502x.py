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

# Verification test of LIAISON_GROUP in HPC mode

from code_aster.Commands import *
from code_aster import CA

test = CA.TestCase()

CA.init("--test")

nProc = CA.MPI.ASTER_COMM_WORLD.Get_size()
rank = CA.MPI.ASTER_COMM_WORLD.Get_rank()

idx = 1

pMesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH", INFO_MED=1)

pMesh = MODI_MAILLAGE(
    reuse=pMesh, MAILLAGE=pMesh, ORIE_PEAU=_F(GROUP_MA_PEAU=("HAUT", "BAS", "DROITE", "GAUCHE"))
)
pMesh = DEFI_GROUP(reuse=pMesh, MAILLAGE=pMesh, CREA_GROUP_NO=_F(TOUT_GROUP_MA="OUI"))

pmodel = AFFE_MODELE(
    MAILLAGE=pMesh, AFFE=(_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI"),)
)

MA = DEFI_MATERIAU(ELAS=_F(E=210000e6, NU=0.3))

pMATE = AFFE_MATERIAU(MAILLAGE=pMesh, AFFE=_F(TOUT="OUI", MATER=MA))

pload0 = AFFE_CHAR_MECA(MODELE=pmodel, DDL_IMPO=(_F(GROUP_MA="GAUCHE", DX=-1.0),))

pload3 = AFFE_CHAR_CINE(MODELE=pmodel, MECA_IMPO=(_F(GROUP_MA="BAS", DY=0.0),))

solverList = [_F(METHODE="MUMPS"), _F(METHODE="PETSC", RESI_RELA=1e-12, PRE_COND="SANS")]

dblLagList = ["OUI", "NON"]

LINST = DEFI_LIST_REEL(VALE=(0.0, 1.0))

LINST2 = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=LINST), ECHEC=_F(SUBD_NIVEAU=5, SUBD_PAS=10))

LINEDEPL = DEFI_FONCTION(NOM_PARA="INST", ABSCISSE=(0, 1), ORDONNEE=(1, 1))


for dblLag in dblLagList:
    for solver in solverList:

        pload1 = AFFE_CHAR_MECA(
            MODELE=pmodel,
            DOUBLE_LAGRANGE=dblLag,
            LIAISON_GROUP=(
                _F(
                    GROUP_MA_1="GAUCHE",
                    DDL_1="DX",
                    GROUP_MA_2="DROITE",
                    DDL_2="DX",
                    COEF_MULT_1=-1.0,
                    COEF_MULT_2=1.0,
                    COEF_IMPO=1,
                ),
            ),
        )

        pload2 = AFFE_CHAR_MECA(
            MODELE=pmodel, DOUBLE_LAGRANGE=dblLag, PRES_REP=(_F(GROUP_MA="HAUT", PRES=1.0),)
        )

        pRESU = MECA_STATIQUE(
            MODELE=pmodel,
            CHAM_MATER=pMATE,
            EXCIT=(_F(CHARGE=pload0), _F(CHARGE=pload1), _F(CHARGE=pload2), _F(CHARGE=pload3)),
            INST=1.0,
            SOLVEUR=solver,
        )

        g = pRESU.getField("DEPL", idx).restrict(["DX"], ["GAUCHE"]).getValues()

        d = pRESU.getField("DEPL", idx).restrict(["DX"], ["DROITE"]).getValues()

        if nProc == 1:
            test.assertAlmostEqual(d[-1], 0)
            test.assertAlmostEqual(g[-1], -1)
        else:
            if rank == 0:
                test.assertAlmostEqual(d[-1], 0)
            else:
                test.assertAlmostEqual(g[-1], -1)

        if dblLag == "OUI":
            pRESUSNL = STAT_NON_LINE(
                MODELE=pmodel,
                CHAM_MATER=pMATE,
                EXCIT=(_F(CHARGE=pload0), _F(CHARGE=pload1), _F(CHARGE=pload2), _F(CHARGE=pload3)),
                INCREMENT=_F(LIST_INST=LINST2),
                COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", TOUT="OUI"),),
                SOLVEUR=solver,
            )

            g = pRESUSNL.getField("DEPL", 1).restrict(["DX"], ["GAUCHE"]).getValues()

            d = pRESUSNL.getField("DEPL", 1).restrict(["DX"], ["DROITE"]).getValues()

            if nProc == 1:
                test.assertAlmostEqual(d[-1], 0)
                test.assertAlmostEqual(g[-1], -1)
            else:
                if rank == 0:
                    test.assertAlmostEqual(d[-1], 0)
                else:
                    test.assertAlmostEqual(g[-1], -1)


# --------------------------------------------------------------
#  We are now testing a LIAISON_GROUP that spans two subdomains
# --------------------------------------------------------------

pload0 = AFFE_CHAR_MECA(MODELE=pmodel, DDL_IMPO=(_F(GROUP_MA="BAS", DY=-1.0),))

pload3 = AFFE_CHAR_CINE(MODELE=pmodel, MECA_IMPO=(_F(GROUP_MA="BAS", DX=0.0),))

for dblLag in dblLagList:
    for solver in solverList:

        pload1 = AFFE_CHAR_MECA(
            MODELE=pmodel,
            DOUBLE_LAGRANGE=dblLag,
            LIAISON_GROUP=(
                _F(
                    GROUP_MA_1="BAS",
                    DDL_1="DY",
                    GROUP_MA_2="HAUT",
                    DDL_2="DY",
                    COEF_MULT_1=-1.0,
                    COEF_MULT_2=1.0,
                    COEF_IMPO=1,
                ),
            ),
        )

        pRESU = MECA_STATIQUE(
            MODELE=pmodel,
            CHAM_MATER=pMATE,
            EXCIT=(_F(CHARGE=pload0), _F(CHARGE=pload1), _F(CHARGE=pload3)),
            INST=1.0,
            SOLVEUR=solver,
        )

        b = pRESU.getField("DEPL", idx).restrict(["DY"], ["BAS"]).getValues()

        h = pRESU.getField("DEPL", idx).restrict(["DY"], ["HAUT"]).getValues()

        if nProc == 1:
            test.assertAlmostEqual(b[-1], -1)
            test.assertAlmostEqual(h[-1], 0)
        else:
            if rank == 0:
                test.assertAlmostEqual(b[-1], -1)
            else:
                test.assertAlmostEqual(h[-1], 0)

        if dblLag == "OUI":
            pRESUSNL = STAT_NON_LINE(
                MODELE=pmodel,
                CHAM_MATER=pMATE,
                EXCIT=(_F(CHARGE=pload0), _F(CHARGE=pload1), _F(CHARGE=pload2), _F(CHARGE=pload3)),
                INCREMENT=_F(LIST_INST=LINST2),
                COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", TOUT="OUI"),),
                SOLVEUR=solver,
            )

            b = pRESUSNL.getField("DEPL", idx).restrict(["DY"], ["BAS"]).getValues()

            h = pRESUSNL.getField("DEPL", idx).restrict(["DY"], ["HAUT"]).getValues()

            if nProc == 1:
                test.assertAlmostEqual(b[-1], -1)
                test.assertAlmostEqual(h[-1], 0)
            else:
                if rank == 0:
                    test.assertAlmostEqual(b[-1], -1)
                else:
                    test.assertAlmostEqual(h[-1], 0)

test.printSummary()

FIN()
