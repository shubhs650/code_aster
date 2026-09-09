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

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

###################################################################
#
#   Solve Poisson problem with HHO (Axi-Symmetry)
#
#   -div(A*grad(u))) = f
#   Continuous:
#   (A * grad u, grad v) = (f, v)
#   with f given and A > 0
#
#   Solution is a polynomial of order k
#   So method of order k should have a null error
#
#   The script to compute solution is given in zzzz512p.datg
#
#   HHO:
#   sum_{T \in Th} (A * GkT(huT), GkT(hvT))_T = (f, v_T)_T
#
###################################################################


# define material

u = {
    "CONSTANTE": FORMULE(VALE="(X-5)", NOM_PARA=("X", "Y")),
    "LINEAIRE": FORMULE(VALE="X*X+Y*Y", NOM_PARA=("X", "Y")),
    "QUADRATIQUE": FORMULE(VALE="X*X*X+Y*Y*Y-X*X*Y", NOM_PARA=("X", "Y")),
    "CUBIQUE": FORMULE(VALE="X*X*X*X+Y*Y*Y*Y-X*X*Y", NOM_PARA=("X", "Y")),
    "QUARTIQUE": FORMULE(VALE="X*X*X*X*X+Y*Y*Y*Y*Y-X*X*Y", NOM_PARA=("X", "Y")),
}

f = {
    "CONSTANTE": FORMULE(VALE="-1./X", NOM_PARA=("X", "Y")),
    "LINEAIRE": FORMULE(VALE="-6", NOM_PARA=("X", "Y")),
    "QUADRATIQUE": FORMULE(VALE="-9*X-2*Y", NOM_PARA=("X", "Y")),
    "CUBIQUE": FORMULE(VALE="-16*X*X-12*Y*Y+4*Y", NOM_PARA=("X", "Y")),
    "QUARTIQUE": FORMULE(VALE="-25*X*X*X-20*Y*Y*Y+4*Y", NOM_PARA=("X", "Y")),
}

mesh0 = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

mesh = DEFI_GROUP(
    reuse=mesh, MAILLAGE=mesh, CREA_GROUP_MA=_F(NOM="2D", TOUT="OUI", TYPE_MAILLE="2D")
)

coeff = DEFI_MATERIAU(THER=_F(LAMBDA=1.0, RHO_CP=1.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

formu = ["CONSTANTE", "LINEAIRE", "QUADRATIQUE", "CUBIQUE", "QUARTIQUE"]

delta_l2 = {
    "CONSTANTE": 8e-4,
    "LINEAIRE": 5e-5,
    "QUADRATIQUE": 6e-7,
    "CUBIQUE": 3e-8,
    "QUARTIQUE": 8e-10,
}
delta_ml2 = {
    "CONSTANTE": 5e-7,
    "LINEAIRE": 1e-9,
    "QUADRATIQUE": 1e-9,
    "CUBIQUE": 1e-10,
    "QUARTIQUE": 1e-10,
}

for form in formu:
    # define model
    model = AFFE_MODELE(
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", MODELISATION="AXIS_HHO", FORMULATION=form, PHENOMENE="THERMIQUE"),
    )

    bc = AFFE_CHAR_CINE_F(MODELE=model, THER_IMPO=_F(TEMP=u[form], GROUP_MA="BOUNDARIES"), INFO=1)
    load = AFFE_CHAR_THER_F(MODELE=model, SOURCE=_F(GROUP_MA="2D", SOUR=f[form]))

    # define discrete object
    phys_pb = CA.PhysicalProblem(model, mater)
    phys_pb.addDirichletBC(bc)
    phys_pb.addLoad(load)

    disc_comp = CA.DiscreteComputation(phys_pb)

    # compute DOF numbering
    phys_pb.computeDOFNumbering()

    # compute (GkT(huT), GkT(hvT))_T
    matEK = disc_comp.getLinearStiffnessMatrix()
    matK = CA.AssemblyMatrixTemperatureReal(phys_pb)
    matK.assemble(matEK, phys_pb.getListOfLoads())

    # compute (u_T, v_T) _T
    matEM = disc_comp.getMassMatrix()
    matM = CA.AssemblyMatrixTemperatureReal(phys_pb)
    matM.assemble(matEM, phys_pb.getListOfLoads())

    # compute (f, v_T)_T
    rhs = disc_comp.getVolumetricForces()

    # lhs matrix
    lhs = matK

    # BC
    diriBC = disc_comp.getDirichletBC()

    # solve linear system
    mySolver = CA.MumpsSolver()
    mySolver.factorize(lhs)
    u_sol = mySolver.solve(rhs, diriBC)

    del mySolver

    hho = CA.HHO(phys_pb)

    # project function
    u_hho = hho.projectOnHHOSpace(u[form])

    u_diff = u_hho - u_sol

    if form == "CONSTANTE":
        # not enougth quadrature point
        delta = 1e-3
    test.assertAlmostEqual(u_diff.norm("NORM_2") / u_hho.norm("NORM_2"), 0.0, delta=delta_l2[form])

    l2_diff = (matM * u_diff).dot(u_diff)
    l2_ref = (matM * u_hho).dot(u_hho)
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=delta_ml2[form])

    # project HHO solution
    h1_field = hho.projectOnLagrangeSpace(u_sol)
    hho_field = hho.projectOnHHOSpace(h1_field)

    # to test TEMP_ELGA
    f_elga = hho.evaluateAtQuadraturePoints(u_sol)
    f_hho = hho.projectOnHHOCellSpace(f_elga)
    u_diff = f_hho - u_sol
    l2_diff = (matM * u_diff).dot(u_diff)
    delta = 1e-7
    if form == "CONSTANTE":
        # not enougth quadrature point
        delta = 2.0e-5
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=delta)

test.printSummary()

FIN()
