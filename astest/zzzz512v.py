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

# --------------------------------------------------------------------
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modifZ
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
#   Patch test with analytical solution
#   Solution is a polynomial of order k
#   So method of order k should have a null error
#
#   The script to compute solution is given in zzzz512v.datg
#
####################################################################################

E = 200000.0
Nu = 0.3

lamb = E * Nu / (1 + Nu) / (1 - 2 * Nu)
mu = E / 2 / (1 + Nu)

uR = {
    "CONSTANTE": FORMULE(VALE="X", NOM_PARA=("X", "Y")),
    "LINEAIRE": FORMULE(VALE="X*(1+X+Y)", NOM_PARA=("X", "Y")),
    "QUADRATIQUE": FORMULE(VALE="X*(1+X*X+Y*Y+X*Y)", NOM_PARA=("X", "Y")),
    "CUBIQUE": FORMULE(VALE="X*(1+X*X*X+Y*Y*Y+X*Y)", NOM_PARA=("X", "Y")),
    "QUARTIQUE": FORMULE(VALE="X*(1+X*X*X*X+Y*Y*Y*Y+X*Y)", NOM_PARA=("X", "Y")),
}
uZ = {
    "CONSTANTE": FORMULE(VALE="Y", NOM_PARA=("X", "Y")),
    "LINEAIRE": FORMULE(VALE="Y*(Y+1)", NOM_PARA=("X", "Y")),
    "QUADRATIQUE": FORMULE(VALE="Y*(Y*Y - X*X + 1)", NOM_PARA=("X", "Y")),
    "CUBIQUE": FORMULE(VALE="Y*(Y*Y*Y-X*X*X + 1)", NOM_PARA=("X", "Y")),
    "QUARTIQUE": FORMULE(VALE="Y*(Y*Y*Y*Y-X*X*X*X + 1)", NOM_PARA=("X", "Y")),
}

zero = FORMULE(VALE="0", NOM_PARA=("X", "Y"))

fR = {
    "CONSTANTE": zero,
    "LINEAIRE": FORMULE(VALE="-3*lamb-6*mu", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu),
    "QUADRATIQUE": FORMULE(
        VALE="-6*lamb*X-3*lamb*Y-16*mu*X-6*mu*Y", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu
    ),
    "CUBIQUE": FORMULE(
        VALE="-12*lamb*X*X-3*lamb*Y-27*mu*X*X -6*mu*X*Y -6*mu*Y",
        NOM_PARA=("X", "Y"),
        lamb=lamb,
        mu=mu,
    ),
    "QUARTIQUE": FORMULE(
        VALE="-20*lamb*X*X*X-3*lamb*Y-44*mu*X*X*X -12*mu*X*Y*Y -6*mu*Y",
        NOM_PARA=("X", "Y"),
        lamb=lamb,
        mu=mu,
    ),
}
fZ = {
    "CONSTANTE": zero,
    "LINEAIRE": FORMULE(VALE="-4*lamb-6*mu", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu),
    "QUADRATIQUE": FORMULE(
        VALE="-3*lamb*X-10*lamb*Y-3*mu*X-12*mu*Y", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu
    ),
    "CUBIQUE": FORMULE(
        VALE="-3*lamb*X - 18*lamb*Y*Y+9*mu*X*Y-3*mu*X-30*mu*Y*Y",
        NOM_PARA=("X", "Y"),
        lamb=lamb,
        mu=mu,
    ),
    "QUARTIQUE": FORMULE(
        VALE="-3*lamb*X-28*lamb*Y*Y*Y+16*mu*X*X*Y -3*mu*X -48*mu*Y*Y*Y",
        NOM_PARA=("X", "Y"),
        lamb=lamb,
        mu=mu,
    ),
}

mesh0 = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

mesh = MODI_MAILLAGE(reuse=mesh, MAILLAGE=mesh, ORIE_PEAU=_F(GROUP_MA_PEAU=("BOUNDARIES")))

mesh = DEFI_GROUP(
    reuse=mesh, MAILLAGE=mesh, CREA_GROUP_MA=_F(NOM="2D", TOUT="OUI", TYPE_MAILLE="2D")
)

# define material
coeff = DEFI_MATERIAU(ELAS=_F(E=E, NU=Nu, RHO=1.0), HHO=_F(COEF_STAB=2 * mu))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

for form in ["LINEAIRE", "QUADRATIQUE", "CUBIQUE", "QUARTIQUE"]:
    model = AFFE_MODELE(
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", MODELISATION="AXIS_HHO", FORMULATION=form, PHENOMENE="MECANIQUE"),
    )

    bc = AFFE_CHAR_CINE_F(
        MODELE=model, MECA_IMPO=_F(GROUP_MA="BOUNDARIES", DX=uR[form], DY=uZ[form])
    )

    load = AFFE_CHAR_MECA_F(
        MODELE=model,
        FORCE_INTERNE=_F(GROUP_MA="2D", FX=fR[form], FY=fZ[form]),
        PRES_REP=_F(GROUP_MA="RIGHT", PRES=zero),
        FORCE_CONTOUR=_F(GROUP_MA="RIGHT", FX=zero, FY=zero),
    )

    # fake load - for coverage
    load0 = AFFE_CHAR_MECA(
        MODELE=model,
        FORCE_INTERNE=_F(GROUP_MA="2D", FX=0.0, FY=0.0),
        FORCE_CONTOUR=_F(GROUP_MA="RIGHT", FX=0.0, FY=0.0),
    )

    # solve linear system
    LREEL = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1, NOMBRE=1))

    resu = STAT_NON_LINE(
        MODELE=model,
        CHAM_MATER=mater,
        INCREMENT=_F(LIST_INST=LREEL),
        EXCIT=(_F(CHARGE=bc), _F(CHARGE=load), _F(CHARGE=load0)),
    )

    u_sol = resu.getField("DEPL", para="INST", value=1.0)

    # define discrete object
    phys_pb = CA.PhysicalProblem(model, mater)
    phys_pb.addDirichletBC(bc)
    phys_pb.computeDOFNumbering()

    hho = CA.HHO(phys_pb)

    # project function
    u_hho = hho.projectOnHHOSpace([uR[form], uZ[form]])

    u_diff = u_hho - u_sol

    test.assertAlmostEqual(u_diff.norm("NORM_2") / u_hho.norm("NORM_2"), 0.0, delta=5e-6)

    dc = CA.DiscreteComputation(phys_pb)
    mass = dc.getMassMatrix(assembly=True)
    l2_diff = (mass * u_diff).dot(u_diff)
    l2_ref = (mass * u_hho).dot(u_hho)
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=1e-10)

    # compute (u_T, v_T) _T
    matEM = CALC_MATR_ELEM(MODELE=model, OPTION="MASS_MECA", CHAM_MATER=mater)

    matM = CA.AssemblyMatrixDisplacementReal(phys_pb)
    matM.assemble(matEM, phys_pb.getListOfLoads())

    l2_diff = (matM * u_diff).dot(u_diff)
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=1e-10)

    # for keyword coverage
    rigi_elem = CALC_MATR_ELEM(MODELE=model, OPTION="RIGI_MECA", CHAM_MATER=mater)
    amor_elem = CALC_MATR_ELEM(
        MODELE=model, OPTION="AMOR_MECA", CHAM_MATER=mater, RIGI_MECA=rigi_elem, MASS_MECA=matEM
    )

FIN()
