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

from code_aster.CA import MPI
from code_aster.Commands import *
from code_aster import CA

# Convergence verification in STAT_NON_LINE and MECA_STATIQUE with LIAISON_MAIL
# Sequential/Parallel comparison

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()

import numpy as np

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()


def buildCompleteFieldOnNodes(field):
    """
    Build complete (over processes) field on nodes
    Arguments:
        field (SimpleFieldOnNodes): field to complete

    Returns:
        list: list of lists containing all field values over processes on global numbering
    """
    field.updateValuePointers()
    mesh = field.getMesh()
    lTGN = mesh.getLocalToGlobalNodeIds()
    maxNodes = max(lTGN)
    maxNodes = MPI.ASTER_COMM_WORLD.allreduce(maxNodes, MPI.MAX)

    innerNodesSet = set(mesh.getInnerNodes())
    values, mask = field.getValues()
    nbNode = field.getNumberOfNodes()
    nbCmp = field.getNumberOfComponents()

    completeField = np.zeros((maxNodes + 1, nbCmp))
    cmpt = 0
    for idNode in range(nbNode):
        if idNode in innerNodesSet:
            globNodeId = lTGN[idNode]
            toAdd = []
            for iCmp in range(nbCmp):
                toAdd.append(values[idNode, iCmp])
            completeField[globNodeId] = np.array(toAdd)
        cmpt += 1

    return MPI.ASTER_COMM_WORLD.allreduce(completeField, MPI.SUM)


# First: Parallel
filename = "zzzz155m.med"
parallelMesh = CA.ParallelMesh()
parallelMesh.readMedFile(filename)

materialDef = DEFI_MATERIAU(ELAS=_F(E=9.84e4, NU=0.3))

parallelModel = AFFE_MODELE(
    MAILLAGE=parallelMesh,
    AFFE=(_F(GROUP_MA=("vol1", "vol2"), PHENOMENE="MECANIQUE", MODELISATION="3D"),),
)

parallelMaterial = AFFE_MATERIAU(
    MAILLAGE=parallelMesh, AFFE=(_F(GROUP_MA=("vol1", "vol2"), MATER=materialDef),)
)

parallelBloc = AFFE_CHAR_CINE(
    MODELE=parallelModel, MECA_IMPO=(_F(GROUP_MA=("gauche"), DX=0, DY=0, DZ=0),)
)

parallelLoad = AFFE_CHAR_CINE(
    MODELE=parallelModel, MECA_IMPO=(_F(GROUP_MA=("droite"), DX=1, DY=1),)
)

parallelGlue = AFFE_CHAR_MECA(
    MODELE=parallelModel, LIAISON_MAIL=(_F(GROUP_NO_ESCL="collage", GROUP_MA_MAIT="vol2"),)
)

LINST = DEFI_LIST_REEL(DEBUT=0, INTERVALLE=(_F(JUSQU_A=2, NOMBRE=3),))

LINST2 = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=LINST), ECHEC=_F(SUBD_NIVEAU=5, SUBD_PAS=10))

LINEDEPL = DEFI_FONCTION(NOM_PARA="INST", ABSCISSE=(0, 1, 2), ORDONNEE=(0, 8, 0))

parallelResuMS = MECA_STATIQUE(
    MODELE=parallelModel,
    CHAM_MATER=parallelMaterial,
    EXCIT=(
        _F(CHARGE=parallelBloc),
        _F(CHARGE=parallelGlue),
        _F(CHARGE=parallelLoad, FONC_MULT=LINEDEPL),
    ),
    LIST_INST=LINST,
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_DP"),
)

parallelResuSNL = STAT_NON_LINE(
    MODELE=parallelModel,
    CHAM_MATER=parallelMaterial,
    EXCIT=(
        _F(CHARGE=parallelBloc),
        _F(CHARGE=parallelGlue),
        _F(CHARGE=parallelLoad, FONC_MULT=LINEDEPL),
    ),
    INCREMENT=_F(LIST_INST=LINST2),
    COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", GROUP_MA=("vol1", "vol2")),),
    CONVERGENCE=_F(ITER_GLOB_MAXI=100, RESI_GLOB_RELA=1e-6),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_DP"),
)

# With MUMPS
parallelResuMSMumps = MECA_STATIQUE(
    MODELE=parallelModel,
    CHAM_MATER=parallelMaterial,
    EXCIT=(
        _F(CHARGE=parallelBloc),
        _F(CHARGE=parallelGlue),
        _F(CHARGE=parallelLoad, FONC_MULT=LINEDEPL),
    ),
    LIST_INST=LINST,
    SOLVEUR=_F(METHODE="MUMPS"),
)

# With MUMPS
parallelResuSNLMumps = STAT_NON_LINE(
    MODELE=parallelModel,
    CHAM_MATER=parallelMaterial,
    EXCIT=(
        _F(CHARGE=parallelBloc),
        _F(CHARGE=parallelGlue),
        _F(CHARGE=parallelLoad, FONC_MULT=LINEDEPL),
    ),
    INCREMENT=_F(LIST_INST=LINST2),
    COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", GROUP_MA=("vol1", "vol2")),),
    CONVERGENCE=_F(ITER_GLOB_MAXI=100, RESI_GLOB_RELA=1e-6),
    SOLVEUR=_F(METHODE="MUMPS"),
)

# Then: Sequential
sequentialMesh = CA.Mesh()
sequentialMesh.readMedFile(filename)

sequentialModel = AFFE_MODELE(
    MAILLAGE=sequentialMesh,
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
    AFFE=(_F(GROUP_MA=("vol1", "vol2"), PHENOMENE="MECANIQUE", MODELISATION="3D"),),
)

sequentialMaterial = AFFE_MATERIAU(
    MAILLAGE=sequentialMesh, AFFE=(_F(GROUP_MA=("vol1", "vol2"), MATER=materialDef),)
)

sequentialBloc = AFFE_CHAR_CINE(
    MODELE=sequentialModel, MECA_IMPO=(_F(GROUP_MA=("gauche"), DX=0, DY=0, DZ=0),)
)

sequentialLoad = AFFE_CHAR_CINE(
    MODELE=sequentialModel, MECA_IMPO=(_F(GROUP_MA=("droite"), DX=1, DY=1),)
)

sequentialGlue = AFFE_CHAR_MECA(
    MODELE=sequentialModel, LIAISON_MAIL=(_F(GROUP_NO_ESCL="collage", GROUP_MA_MAIT="vol2"),)
)

sequentialResuMS = MECA_STATIQUE(
    MODELE=sequentialModel,
    CHAM_MATER=sequentialMaterial,
    EXCIT=(
        _F(CHARGE=sequentialBloc),
        _F(CHARGE=sequentialGlue),
        _F(CHARGE=sequentialLoad, FONC_MULT=LINEDEPL),
    ),
    LIST_INST=LINST,
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_DP"),
)

sequentialResuSNL = STAT_NON_LINE(
    MODELE=sequentialModel,
    CHAM_MATER=sequentialMaterial,
    EXCIT=(
        _F(CHARGE=sequentialBloc),
        _F(CHARGE=sequentialGlue),
        _F(CHARGE=sequentialLoad, FONC_MULT=LINEDEPL),
    ),
    INCREMENT=_F(LIST_INST=LINST2),
    COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", GROUP_MA=("vol1", "vol2")),),
    CONVERGENCE=_F(ITER_GLOB_MAXI=100, RESI_GLOB_RELA=1e-6),
    SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_DP"),
)

parallelDispSNL = parallelResuSNL.getField("DEPL", 2).toSimpleFieldOnNodes()
# Parallel field completion
sParallelDispSNL = buildCompleteFieldOnNodes(parallelDispSNL)

sequentialDispSNL = sequentialResuSNL.getField("DEPL", 2).toSimpleFieldOnNodes()

parallelDispMS = parallelResuMS.getField("DEPL", 2).toSimpleFieldOnNodes()
# Parallel field completion
sParallelDispMS = buildCompleteFieldOnNodes(parallelDispMS)

parallelDispSNLM = parallelResuSNLMumps.getField("DEPL", 2).toSimpleFieldOnNodes()
# Parallel field completion
sParallelDispSNLM = buildCompleteFieldOnNodes(parallelDispSNLM)

parallelDispMSM = parallelResuMSMumps.getField("DEPL", 2).toSimpleFieldOnNodes()
# Parallel field completion
sParallelDispMSM = buildCompleteFieldOnNodes(parallelDispMSM)

sequentialDispMS = sequentialResuMS.getField("DEPL", 2).toSimpleFieldOnNodes()

pValuesSNL = sParallelDispSNL
sValuesSNL = sequentialDispSNL.getValues()[0]
pValuesMS = sParallelDispMS
sValuesMS = sequentialDispMS.getValues()[0]
pValuesSNLM = sParallelDispSNLM
pValuesMSM = sParallelDispMSM

# 1e-12 sequential/parallel SNL comparison after 2 iterations
for initArray, newArray in zip(pValuesSNL, sValuesSNL):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, 12)

# 1e-12 sequential/parallel MS comparison after 2 iterations
for initArray, newArray in zip(pValuesMS, sValuesMS):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, 12)

# 1e-12 SNL/MS comparison after 2 iterations
for initArray, newArray in zip(pValuesSNL, pValuesMS):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, 12)

# 1e-12 sequential/parallel SNL comparison after 2 iterations
for initArray, newArray in zip(pValuesSNL, pValuesSNLM):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, 12)

# 1e-12 sequential/parallel MS comparison after 2 iterations
for initArray, newArray in zip(pValuesMS, pValuesMSM):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, 12)

FIN()
