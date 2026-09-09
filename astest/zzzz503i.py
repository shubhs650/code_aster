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
from code_aster.Utilities import PETSc

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()

mesh = LIRE_MAILLAGE(FORMAT="MED")
model = AFFE_MODELE(
    MAILLAGE=mesh,
    AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"),
    DISTRIBUTION=_F(METHODE="SOUS_DOMAINE", NB_SOUS_DOMAINE=2),
)

material = DEFI_MATERIAU(ELAS=_F(E=100.0, NU=0.3))
mat_field = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=material))

dirichlet = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=_F(GROUP_MA="Encast", DX=1.0, DY=2.0))

phys_pb = CA.PhysicalProblem(model, mat_field)
phys_pb.addLoad(dirichlet)

# Numbering of equations
phys_pb.computeDOFNumbering()

# Create behaviour
phys_pb.computeBehaviourProperty(COMPORTEMENT=(_F(RELATION="ELAS", TOUT="OUI"),))

# Create discrete problem
disc_comp = CA.DiscreteComputation(phys_pb)

# Build matrices
matr_elem = disc_comp.getElasticStiffnessMatrix()
matr_elem_dual = disc_comp.getDualStiffnessMatrix()

matrAsse = CA.AssemblyMatrixDisplacementReal()
matrAsse.setDOFNumbering(phys_pb.getDOFNumbering())
matrAsse.assemble([matr_elem, matr_elem_dual])


# ------------------------------------
# tests in local numbering
numeDDL = phys_pb.getDOFNumbering()
physicalRows = numeDDL.getPhysicalDOFs(local=True)
test.assertListEqual(physicalRows, [2, 3, 8, 9, 12, 13, 14, 15, 18, 19, 24, 25, 28, 29, 30, 31])
multipliersRows = numeDDL.getLagrangeDOFs(local=True)
test.assertListEqual(multipliersRows, [0, 1, 4, 5, 6, 7, 10, 11, 16, 17, 20, 21, 22, 23, 26, 27])
test.assertTrue(numeDDL.useLagrangeDOF())
test.assertFalse(numeDDL.useSingleLagrangeDOF())
test.assertEqual(numeDDL.getComponents(), ["DX", "DY", "LAGR:DX", "LAGR:DY"])
test.assertEqual(numeDDL.getComponentFromNode(0, local=True), ["DX", "DY"])
test.assertEqual(numeDDL.getDOFsAssociatedToComponent("LAGR:DX"), [0, 4, 6, 10, 16, 20, 22, 26])
test.assertEqual(numeDDL.getDOFsAssociatedToComponent("LAGR:DY"), [1, 5, 7, 11, 17, 21, 23, 27])
test.assertEqual(numeDDL.getNodeFromDOF(0, local=True), 0)
test.assertFalse(numeDDL.isPhysicalDOF(0, local=True))
test.assertEqual(numeDDL.getNumberOfDOFs(local=True), 32)
test.assertEqual(numeDDL.getNumberOfDOFs(local=False), 32)
test.assertEqual(numeDDL.getPhysicalQuantity(), "DEPL_R")

# TODO A compléter après correction de issue32247
# ------------------------------------
# tests in global numbering
# physicalRows = numeDDL.getPhysicalDOFs(local=False)
# test.assertListEqual(physicalRows,  [numeDDL.localToGlobalDOF(d)
#                                      for d in [i*2+j for i in range(mesh.getNumberOfNodes())
#                                      for j in range(2)]])
# test.assertListEqual(physicalRows,  [[0, 1, 2, 3, 4, 5, 12, 13, 6, 7, 14, 15],
#                                      [8, 9, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15]][rank])

# -------------------------------------------
# MatrixScaler verification section - No context manager

from code_aster.LinearAlgebra import MatrixScaler
import numpy as np
from scipy.linalg import norm

# Compute reference solution
rhs = CA.FieldOnNodesReal(numeDDL)
rhs.setValues(1)

mySolver = CA.MumpsSolver()
mySolver.factorize(matrAsse)

ref_sol = mySolver.solve(rhs)

# Compute scaled solution
S = MatrixScaler.MatrixScaler()

nt = PETSc.NormType.NORM_INFINITY
test.assertAlmostEqual(matrAsse.norm("NORM_INFINITY"), 1527.7777777794063)
test.assertAlmostEqual(matrAsse.toPetsc().norm(nt), norm(matrAsse.toNumpy(), np.inf))

S.computeScaling(matrAsse, merge_dof=[["DX", "DY", "DZ"], ["DRX", "DRY", "DRZ"]], verbose=True)
S.scaleMatrix(matrAsse)

test.assertAlmostEqual(matrAsse.toPetsc().norm(nt), 1.0977480819609067)
test.assertAlmostEqual(matrAsse.toPetsc().norm(nt), norm(matrAsse.toNumpy(), np.inf))

S.scaleRHS(rhs)

mySolver.factorize(matrAsse)
scaled_sol = mySolver.solve(rhs)

# The solution *must* be unscaled !
S.unscaleSolution(scaled_sol)

test.assertAlmostEqual(scaled_sol.toPetsc().norm(nt), ref_sol.toPetsc().norm(nt))

# MatrixScaler verification section - with context manager
rhs = CA.FieldOnNodesReal(numeDDL)
rhs.setValues(1)

mySolver = CA.MumpsSolver()
mySolver.factorize(matrAsse)

ref_sol = mySolver.solve(rhs)

with MatrixScaler.matrixScaler(
    matrAsse,
    rhs,
    merge_dof=[["DX", "DY", "DZ"], ["DRX", "DRY", "DRZ"]],
    scaling_type=True,
    verbose=True,
) as scaling:
    mySolver.factorize(matrAsse)
    scaled_sol = mySolver.solve(rhs)
    # The solution *must* be unscaled !
    scaling.unscaleSolution(scaled_sol)

test.assertAlmostEqual(scaled_sol.toPetsc().norm(nt), ref_sol.toPetsc().norm(nt))

test.printSummary()


FIN()
