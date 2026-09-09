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

import os
import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True

except ImportError:
    HAS_MATPLOTLIB = False

from math import sqrt, log

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

###################################################################################
#
#   Analytical solution
#   -laplacian(u) = f
#   u = sin(Pi*x)*sin(Pi*y)
#   f = 2.0*Pi*Pi*sin(Pi*x)*sin(Pi*y)
#
#   Weak form: (grad u, grad v) = (f,v)
#   HHO unknowns : huT = (uT, udT)
#
####################################################################################


def diameter_cell(coor, nodes):
    nb_nodes = len(nodes)
    pts = np.array([coor.getNode(node).getValues() for node in nodes])

    return max(
        np.linalg.norm(pts[i] - pts[j]) for i in range(nb_nodes) for j in range(i + 1, nb_nodes)
    )


def diameter(mesh):
    mesh_lin = mesh.convertToLinear()
    nbCells = mesh_lin.getNumberOfCells()
    coor = mesh_lin.getCoordinates()
    connec = mesh_lin.getConnectivity()

    diam = -1.0
    for c_id in range(nbCells):
        diam = max(diam, diameter_cell(coor, connec[c_id]))

    return diam


# number of refinement
nb_reff = 6

E = 200000.0
Nu = 0.3

lamb = E * Nu / (1 + Nu) / (1 - 2 * Nu)
mu = E / 2 / (1 + Nu)

# define analytical solution
u_X = FORMULE(VALE="0.8*X*X*Y*Y-X+1", NOM_PARA=("X", "Y"))
u_Y = FORMULE(VALE="-Y*Y*Y*Y+0.1*X-0.5*X*Y", NOM_PARA=("X", "Y"))

# define load function
f_X = FORMULE(
    VALE="-lamb*(1.6*Y*Y - 0.5) - 3.2*mu*Y*Y - mu*(1.6*X*X-0.5)",
    NOM_PARA=("X", "Y"),
    lamb=lamb,
    mu=mu,
)
f_Y = FORMULE(
    VALE="-Y*(lamb*(3.2*X - 12*Y) + 3.2*mu*X - 24.0*mu*Y)", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu
)

# error save
error = {}
conv_order = {}

# initial_mesh - domain [0,1]^2 - 4 quads
mesh0_quad = CA.Mesh.buildSquare(refine=1)
# split in triangular mesh
mesh0_tri = CREA_MAILLAGE(MAILLAGE=mesh0_quad, MODI_MAILLE=_F(TOUT="OUI", OPTION="QUAD_TRIA3"))
# convert for hho-cells
mesh0_hho = CREA_MAILLAGE(MAILLAGE=mesh0_tri, MODI_HHO=_F(TOUT="OUI"))

for order in ("LINEAIRE", "QUADRATIQUE"):
    error[order] = {"h": [], "L2": [], "H1": []}
    mesh = mesh0_hho
    for i_reff in range(nb_reff):
        ## DEFINE PROBLEM
        # create mesh - refine previous mesh
        mesh = mesh.refine(1)

        # define material
        coeff = DEFI_MATERIAU(ELAS=_F(E=E, NU=Nu, RHO=1.0), HHO=_F(COEF_STAB=2 * mu))

        # apply material on mesh
        mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

        # define finite element model
        model = AFFE_MODELE(
            MAILLAGE=mesh,
            AFFE=_F(
                TOUT="OUI", MODELISATION="D_PLAN_HHO", FORMULATION=order, PHENOMENE="MECANIQUE"
            ),
        )

        # define Dirichlet BC
        bc = AFFE_CHAR_CINE_F(
            MODELE=model, MECA_IMPO=_F(GROUP_MA=("RIGHT", "LEFT", "TOP", "BOTTOM"), DX=u_X, DY=u_Y)
        )

        # define external load
        load = AFFE_CHAR_MECA_F(
            MODELE=model, FORCE_INTERNE=_F(GROUP_MA=("SURFACE"), FX=f_X, FY=f_Y)
        )

        ## COMPUTE DISCRETE SOLUTION
        # define physical problem
        phys_pb = CA.PhysicalProblem(model, mater)
        phys_pb.addLoad(load)
        phys_pb.addDirichletBC(bc)

        # compute DOF numbering
        phys_pb.computeDOFNumbering()

        # create discrete computation
        disc_comp = CA.DiscreteComputation(phys_pb)

        # compute rigidity matrix: (lambda * GkT(huT), GkT(hvT))_T + lambda * stab(huT, hvT)
        rigidity = disc_comp.getLinearStiffnessMatrix(assembly=True)

        # compute external load: (f, vT)_T
        rhs = disc_comp.getVolumetricForces()

        # compute Dirichlet BC to apply
        diriBCs = disc_comp.getDirichletBC()

        # define linear solver - MUMPS
        mySolver = CA.PetscSolver(PRE_COND="LDLT_SP", RESI_RELA=1e-14)

        # factorize and solve
        mySolver.factorize(rigidity)
        u_hho = mySolver.solve(rhs, diriBCs)

        ## COMPUTE ERROR
        # create hho handler
        hho = CA.HHO(phys_pb)

        # Project analytical solution on HHO space
        u_proj = hho.projectOnHHOSpace([u_X, u_Y])

        # compute difference
        u_diff = u_hho - u_proj

        # Compute mass matrix M = (vT, wT) - RHO_CP == 1 in DEFI_MATERIAU
        mass = disc_comp.getMassMatrix(assembly=True)

        # compute L2 and H1-errors
        error[order]["h"].append(diameter(mesh))
        error[order]["L2"].append(sqrt((mass * u_diff).dot(u_diff)))
        error[order]["H1"].append(sqrt((rigidity * u_diff).dot(u_diff)))

# compute convergence order
for order in ("LINEAIRE", "QUADRATIQUE"):
    conv_order[order] = {}
    for norm in ("L2", "H1"):
        xlog = [log(h) for h in error[order]["h"]]
        ylog = [log(y) for y in error[order][norm]]
        A = np.vstack([xlog, np.ones(len(xlog))]).T
        m, c = np.linalg.lstsq(A, ylog, rcond=None)[0]
        conv_order[order][norm] = [m, c]

    print(
        "Convergence order %s -> L2-norm: %f and H1-norm: %f"
        % (order, conv_order[order]["L2"][0], conv_order[order]["H1"][0])
    )

    # test convergence order
    test.assertAlmostEqual(
        conv_order[order]["L2"][0],
        {"LINEAIRE": 2.7188706153832713, "QUADRATIQUE": 3.9300875165659875}[order],
        delta=1e-4,
    )
    test.assertAlmostEqual(
        conv_order[order]["H1"][0],
        {"LINEAIRE": 1.8901992494377016, "QUADRATIQUE": 2.9764630230205915}[order],
        delta=1e-4,
    )

# plot figure with matplotlib
if HAS_MATPLOTLIB and os.getenv("DISPLAY"):

    # disable floating point exceptions from matplotlib
    with CA.disable_fpe():

        ylim = {"L2": [1e-10, 0.1], "H1": [1e-7, 1.0]}
        the_conv = {
            "LINEAIRE": {"L2": {"o": 3, "c": 2e-2}, "H1": {"o": 2, "c": 0.2}},
            "QUADRATIQUE": {"L2": {"o": 4, "c": 1e-3}, "H1": {"o": 3, "c": 0.02}},
        }
        # plot the data
        for norm in ("L2", "H1"):
            plt.plot(1, 1, 1)

            plt.xscale("log")
            plt.yscale("log")

            # set the limits
            plt.xlim([0.02, 2])
            plt.ylim(ylim[norm])
            plt.xlabel("mesh-size")
            plt.ylabel("%s-error" % norm)

            plt.plot(
                error["LINEAIRE"]["h"],
                error["LINEAIRE"][norm],
                marker="o",
                color="tab:blue",
                label="k=1, computed",
            )
            m, c = conv_order["LINEAIRE"][norm]
            plt.plot(
                error["LINEAIRE"]["h"],
                [
                    the_conv["LINEAIRE"][norm]["c"] * h ** the_conv["LINEAIRE"][norm]["o"]
                    for h in error["LINEAIRE"]["h"]
                ],
                "k--",
                label="k=1, theorical",
                color="tab:blue",
            )

            plt.plot(
                error["QUADRATIQUE"]["h"],
                error["QUADRATIQUE"][norm],
                marker="o",
                color="tab:orange",
                label="k=2, computed",
            )
            m, c = conv_order["QUADRATIQUE"][norm]
            plt.plot(
                error["QUADRATIQUE"]["h"],
                [
                    the_conv["QUADRATIQUE"][norm]["c"] * h ** the_conv["QUADRATIQUE"][norm]["o"]
                    for h in error["LINEAIRE"]["h"]
                ],
                "k--",
                label="k=2, theorical",
                color="tab:orange",
            )

            plt.legend()
            plt.title("%s convergence error for HHO" % norm)
            plt.show()

            # save plot
            # savedir = "/tmp/" or os.getcwd()
            # plt.savefig(os.path.join(savedir, "%s_error.png" % norm))
            plt.clf()

# close
CA.close()
