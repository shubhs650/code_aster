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
from code_aster.CA import MPI

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()

filename = "zzzz155e.med"
filename2 = "zzzz155e.23"

mesh = CA.ParallelMesh()
mesh.readMedFile(filename, deterministic=True)
mesh.printMedFile("test.med", mode="w")

mesh2 = CA.ParallelMesh()
mesh2.readMedFile(filename2, deterministic=True)
# overwrite file
mesh2.printMedFile("test.med", mode="w")

mesh2 = DEFI_GROUP(reuse=mesh2, MAILLAGE=mesh2, CREA_GROUP_NO=_F(TOUT_GROUP_MA="OUI"))

vol1 = set(mesh2.getNodes("VOL1"))
vol2 = set(mesh2.getNodes("VOL2"))
vol3 = vol2 - vol1
mesh2.setGroupOfNodes("VOL3", list(vol3))

modRed = AFFE_MODELE(
    MAILLAGE=mesh2, AFFE=(_F(GROUP_MA=("VOL1", "VOL2"), PHENOMENE="MECANIQUE", MODELISATION="3D"),)
)

f1 = CA.FieldOnNodesReal(modRed)
f1.setValues(0.0)

f1Red = f1.restrict(["DX"], ["VOL3"])

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))

cl = AFFE_CHAR_CINE(
    MODELE=model,
    MECA_IMPO=(_F(GROUP_MA="BAS", DX=0.0, DY=0.0, DZ=0.0), _F(GROUP_MA="HAUT", DZ=1.0)),
)

resu = MECA_STATIQUE(MODELE=model, CHAM_MATER=mater, EXCIT=_F(CHARGE=cl), INST=0.0)

LREEL = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1, NOMBRE=1))

resu2 = STAT_NON_LINE(
    MODELE=model, CHAM_MATER=mater, INCREMENT=_F(LIST_INST=LREEL), EXCIT=_F(CHARGE=cl)
)

test.assertTrue(True)

FIN()
