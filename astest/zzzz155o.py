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

CA.init("--test", IGNORE_ALARM="MODELISA8_15")

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()

test = CA.TestCase()

filename = "zzzz155o.med"

refMesh = CA.Mesh()
refMesh.readMedFile(filename)
refCoords = refMesh.getCoordinates()
refG2L = [i for i in range(refMesh.getNumberOfNodes())]

newMesh = CA.ParallelMesh()
newMesh.readMedFile(filename)


listeGroupmaModele = [
    "vis",
    "cloison",
    "renfort",
    "VisAppuiTete",
    "CloisonSymetrieX",
    "CloisonSymetrieZ",
    "RenfortSymetrieX",
    "RenfortBloquage",
    "DiscretsContact",
    "VisBasR",
]

reducedMesh = CREA_MAILLAGE(
    MAILLAGE=newMesh,
    RESTREINT=_F(GROUP_MA=listeGroupmaModele, TOUT_GROUP_MA="OUI", TOUT_GROUP_NO="OUI"),
)

cleanCoords = reducedMesh.getCoordinates()
cleanL2G = reducedMesh.getLocalToGlobalNodeIds()
assert len(cleanL2G) == len(cleanCoords.getValues()) / 3
for locNodeId in range(len(cleanL2G)):
    globNodeId = cleanL2G[locNodeId]
    refLocNodeId = globNodeId
    refNode = refCoords[refLocNodeId]
    node = cleanCoords[locNodeId]
    for iCoord in range(3):
        test.assertAlmostEqual(refNode[iCoord], node[iCoord], 12)

FIN()
