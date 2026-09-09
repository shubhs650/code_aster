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

CA.init("--test")

test = CA.TestCase()

###################################################################################
#
#   Patch test with analytical solution - linear elasticity
#   Solution is a polynomial of order k
#   So method of order k should have a null error
#
#   The script to compute solution is given in zzzz512f.datg
#
####################################################################################

E = 200000.0
Nu = 0.3

lamb = E * Nu / (1 + Nu) / (1 - 2 * Nu)
mu = E / 2 / (1 + Nu)

uX = {
    "CONSTANTE": FORMULE(VALE="0.8*X+1", NOM_PARA=("X", "Y")),
    "LINEAIRE": FORMULE(VALE="0.8*X*Y-X+1", NOM_PARA=("X", "Y")),
    "QUADRATIQUE": FORMULE(VALE="0.8*X*X*Y-X+1", NOM_PARA=("X", "Y")),
    "CUBIQUE": FORMULE(VALE="0.8*X*X*Y*Y-X+1", NOM_PARA=("X", "Y")),
    "QUARTIQUE": FORMULE(VALE="0.8*X*X*Y*Y*X-X+1", NOM_PARA=("X", "Y")),
}
uY = {
    "CONSTANTE": FORMULE(VALE="-Y+0.1*X", NOM_PARA=("X", "Y")),
    "LINEAIRE": FORMULE(VALE="-Y*Y+0.1*X", NOM_PARA=("X", "Y")),
    "QUADRATIQUE": FORMULE(VALE="-Y*Y*Y+0.1*X-0.5*X*Y", NOM_PARA=("X", "Y")),
    "CUBIQUE": FORMULE(VALE="-Y*Y*Y*Y+0.1*X-0.5*X*Y", NOM_PARA=("X", "Y")),
    "QUARTIQUE": FORMULE(VALE="-Y*Y*Y*Y*X+0.1*X-0.5*X*Y", NOM_PARA=("X", "Y")),
}

zero = FORMULE(VALE="0", NOM_PARA=("X", "Y"))


fX = {
    "CONSTANTE": zero,
    "LINEAIRE": zero,
    "QUADRATIQUE": FORMULE(
        VALE="-lamb*(1.6*Y - 0.5) - 3.2*mu*Y + 0.5*mu", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu
    ),
    "CUBIQUE": FORMULE(
        VALE="-lamb*(1.6*Y*Y - 0.5) - 3.2*mu*Y*Y - mu*(1.6*X*X-0.5)",
        NOM_PARA=("X", "Y"),
        lamb=lamb,
        mu=mu,
    ),
    "QUARTIQUE": FORMULE(
        VALE="lamb*(- 4.8*X*Y*Y  + 4.0*Y*Y*Y  + 0.5) - 9.6*mu*X*Y*Y  + mu*(- 1.6*X*X*X  + 4.0*Y*Y*Y  + 0.5)",
        NOM_PARA=("X", "Y"),
        lamb=lamb,
        mu=mu,
    ),
}
fY = {
    "CONSTANTE": zero,
    "LINEAIRE": FORMULE(VALE="1.2*lamb+3.2*mu", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu),
    "QUADRATIQUE": FORMULE(
        VALE="-lamb*(1.6*X - 6.0*Y) - 1.6*mu*X + 12.0*mu*Y", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu
    ),
    "CUBIQUE": FORMULE(
        VALE="-Y*(lamb*(3.2*X - 12*Y) + 3.2*mu*X - 24.0*mu*Y)",
        NOM_PARA=("X", "Y"),
        lamb=lamb,
        mu=mu,
    ),
    "QUARTIQUE": FORMULE(
        VALE="-X*Y*(lamb*(4.8*X - 12.0*Y) + 4.8*mu*X - 24.0*mu*Y)",
        NOM_PARA=("X", "Y"),
        lamb=lamb,
        mu=mu,
    ),
}

mesh0 = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

mesh = DEFI_GROUP(
    reuse=mesh, MAILLAGE=mesh, CREA_GROUP_MA=_F(NOM="2D", TOUT="OUI", TYPE_MAILLE="2D")
)
# define material
coeff = DEFI_MATERIAU(ELAS=_F(E=E, NU=Nu, RHO=1.0), HHO=_F(COEF_STAB=2 * mu))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

for form in ["LINEAIRE", "QUADRATIQUE", "CUBIQUE", "QUARTIQUE"]:
    model = AFFE_MODELE(
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", MODELISATION="D_PLAN_HHO", FORMULATION=form, PHENOMENE="MECANIQUE"),
    )

    bc = AFFE_CHAR_CINE_F(
        MODELE=model, MECA_IMPO=_F(GROUP_MA="BOUNDARIES", DX=uX[form], DY=uY[form]), INFO=1
    )

    load = AFFE_CHAR_MECA_F(MODELE=model, FORCE_INTERNE=_F(GROUP_MA="2D", FX=fX[form], FY=fY[form]))

    # solve linear system
    LREEL = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1, NOMBRE=1))

    resu = MECA_NON_LINE(
        MODELE=model,
        CHAM_MATER=mater,
        INCREMENT=_F(LIST_INST=LREEL),
        EXCIT=(_F(CHARGE=bc), _F(CHARGE=load)),
        SOLVEUR=_F(METHODE="PETSC", ALGORITHME="FGMRES", PRE_COND="LDLT_SP", PCENT_PIVOT=20),
    )

    u_sol = resu.getField("DEPL", para="INST", value=1.0)

    # define discrete object
    phys_pb = CA.PhysicalProblem(model, mater)
    phys_pb.addDirichletBC(bc)
    phys_pb.computeDOFNumbering()

    hho = CA.HHO(phys_pb)

    # project function
    u_hho = hho.projectOnHHOSpace([uX[form], uY[form]])

    u_diff = u_hho - u_sol

    test.assertAlmostEqual(u_diff.norm("NORM_2") / u_hho.norm("NORM_2"), 0.0, delta=1e-7)

    dc = CA.DiscreteComputation(phys_pb)
    mass = dc.getMassMatrix(assembly=True)
    l2_diff = (mass * u_diff).dot(u_diff)
    l2_ref = (mass * u_hho).dot(u_hho)
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=1e-10)

    # print("u=", u_hho.getValues())
    # print("uh=", u_sol.getValues())
    # print("diif=", u_diff.getValues())

FIN()
