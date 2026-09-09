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

from code_aster.Utilities import logger
import numpy as N

from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import PETSc

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

from code_aster.LinearAlgebra import MatrixScaler

test = CA.TestCase()

rank = CA.MPI.ASTER_COMM_WORLD.Get_rank()

pMesh = CA.ParallelMesh()
pMesh.readMedFile("mesh004b/%d.med" % rank, partitioned=True)
# MAIL.debugPrint()

MATER = DEFI_MATERIAU(THER=_F(LAMBDA=6.0e9, RHO_CP=1.0))

affectMat = CA.MaterialField(pMesh)
affectMat.addMaterialOnMesh(MATER)
affectMat.build()

MODT = AFFE_MODELE(
    MAILLAGE=pMesh,
    AFFE=_F(TOUT="OUI", PHENOMENE="THERMIQUE", MODELISATION="PLAN"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)


charCine = CA.ThermalDirichletBC(MODT)
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Temp, 0.0, "EncastN")
charCine.build()
charCine.debugPrint()

CHT1 = AFFE_CHAR_THER(MODELE=MODT, FLUX_REP=_F(GROUP_MA="Press", FLUN=10.0), INFO=1)

# study = CA.PhysicalProblem(MODT, affectMat)
# study.addDirichletBC(charCine)
# study.addLoad(CHT1)
# dProblem = CA.DiscreteComputation(study)

vect_elem = CALC_VECT_ELEM(OPTION="CHAR_THER", CHARGE=CHT1)
matr_elem = CALC_MATR_ELEM(OPTION="RIGI_THER", MODELE=MODT, CHAM_MATER=affectMat)

monSolver = CA.PetscSolver(RENUM="SANS", PRE_COND="SANS")

numeDDL = CA.ParallelDOFNumbering()
numeDDL.computeNumbering([matr_elem])
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")

matrAsse = CA.AssemblyMatrixTemperatureReal()
matrAsse.setDOFNumbering(numeDDL)
matrAsse.assemble(matr_elem, charCine)
test.assertEqual(matrAsse.getType(), "MATR_ASSE_TEMP_R")

# retour = vect_elem.assembleWithLoadFunctions( numeDDL )
vecass = ASSE_VECTEUR(VECT_ELEM=vect_elem, NUME_DDL=numeDDL)
vcine = CALC_CHAR_CINE(NUME_DDL=numeDDL, CHAR_CINE=charCine)

matrAsse = FACTORISER(reuse=matrAsse, MATR_ASSE=matrAsse, METHODE="PETSC", PRE_COND="JACOBI")

resu = RESOUDRE(MATR=matrAsse, CHAM_NO=vecass, CHAM_CINE=vcine, ALGORITHME="CR", RESI_RELA=1e-9)

TEST_RESU(
    CHAM_NO=_F(
        CRITERE="ABSOLU",
        GROUP_NO="Noeud1",
        NOM_CMP="TEMP",
        REFERENCE="AUTRE_ASTER",
        CHAM_GD=resu,
        VALE_CALC=8.33333333333333e-10,
        VALE_REFE=8.33333333333333e-10,
    )
)

TEST_RESU(
    CHAM_NO=_F(
        CRITERE="ABSOLU",
        GROUP_NO="Noeud2",
        NOM_CMP="TEMP",
        REFERENCE="AUTRE_ASTER",
        CHAM_GD=resu,
        VALE_CALC=8.33333333333333e-10,
        VALE_REFE=8.33333333333333e-10,
    )
)


# tests in local numbering
physicalRows = numeDDL.getPhysicalDOFs(local=True)
test.assertListEqual(physicalRows, list(range(len(pMesh.getNodes(localNumbering=True)))))
multipliersRows = numeDDL.getLagrangeDOFs(local=True)
test.assertListEqual(multipliersRows, [])
test.assertFalse(numeDDL.useLagrangeDOF())
test.assertFalse(numeDDL.useSingleLagrangeDOF())
test.assertEqual(numeDDL.getComponents(), ["TEMP"])
test.assertEqual(numeDDL.getComponentFromNode(0, local=True), ["TEMP"])
test.assertEqual(numeDDL.getNodeFromDOF(0, local=True), 0)
test.assertTrue(numeDDL.isPhysicalDOF(0, local=True))
test.assertEqual(numeDDL.getNumberOfDOFs(local=True), len(pMesh.getNodes(localNumbering=True)))
test.assertEqual(numeDDL.getNumberOfDOFs(local=False), 8)
test.assertEqual(numeDDL.getPhysicalQuantity(), "TEMP_R")
ghostRows = numeDDL.getGhostDOFs(local=True)
test.assertListEqual(ghostRows, [[3, 5], [2, 4]][rank])


# tests in global numbering
physicalRows = numeDDL.getPhysicalDOFs(local=False)
test.assertListEqual(
    physicalRows, [numeDDL.localToGlobalDOF(i) for i in range(pMesh.getNumberOfNodes())]
)
test.assertListEqual(physicalRows, [[0, 1, 2, 6, 3, 7], [4, 5, 2, 6, 3, 7]][rank])

ghostRows = numeDDL.getGhostDOFs(local=False)
test.assertListEqual(ghostRows, [numeDDL.localToGlobalDOF(i) for i in [[3, 5], [2, 4]][rank]])
test.assertListEqual(ghostRows, [[6, 7], [2, 3]][rank])


# Some petsc4py manipulations
pA = matrAsse.toPetsc()

v = PETSc.Viewer()
pA.view(v)
v = PETSc.Viewer().createASCII("test.txt")
v.pushFormat(PETSc.Viewer.Format.ASCII_MATLAB)

rank = pA.getComm().getRank()
rs, re = pA.getOwnershipRange()
ce, _ = pA.getSize()
rows = N.array(list(range(rs, re)), dtype=PETSc.IntType)
cols = N.array(list(range(0, ce)), dtype=PETSc.IntType)
rows = PETSc.IS().createGeneral(rows, comm=pA.getComm())
cols = PETSc.IS().createGeneral(cols, comm=pA.getComm())
(S,) = pA.createSubMatrices(rows, cols)
v = PETSc.Viewer().createASCII("mesh004c_rank" + str(rank) + ".out", comm=S.getComm())
S.view(v)

# Scaling validation
pA_unscaled = matrAsse.toPetsc()

S = MatrixScaler.MatrixScaler()
logger.setLevel(2)
S.computeScaling(matrAsse, merge_dof=[["DX", "DY", "DZ"], ["DRX", "DRY", "DRZ"]])
S.scaleMatrix(matrAsse)

pA_scaled = matrAsse.toPetsc()
pA_scaled.view()
nt = PETSc.NormType.NORM_INFINITY
test.assertAlmostEqual(pA_unscaled.norm(nt), 24666666666.69744, delta=1e-2)
test.assertAlmostEqual(pA_scaled.norm(nt), 1.0)

logger.setLevel(0)

test.printSummary()

FIN()
