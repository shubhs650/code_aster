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


def getConnectivityOfGroup(lgrma, mesh):
    """
    Get connectivity of a list of GROUP_MA
    """
    conn = mesh.getConnectivity()
    con = []
    for grma in lgrma:
        gr_cells = mesh.getCells(grma)
        con.extend([item + 1 for item in conn[cell]] for cell in gr_cells)
    return con


# ---------------------------------------------------------------------------- #
CA.init("--test", "--abort")

test = CA.TestCase()

# read mesh
mesh = CA.Mesh()
mesh.readAsterFile("zzzz254a.mail")

# get initial connectivity
lgrma = ["GR_DEUX", "GR_OPPX", "GR_NOORI"]
con1 = getConnectivityOfGroup(lgrma, mesh)

# apply orientation
mesh = MODI_MAILLAGE(
    reuse=mesh,
    MAILLAGE=mesh,
    ORIE_INTERF_POU=_F(GROUP_MA=("GR_DEUX", "GR_OPPX", "GR_NOORI"), VECT_ORIE=[0.0, 0.0, 1.0]),
)

# get final connectivity
con2 = getConnectivityOfGroup(lgrma, mesh)


# ---------------------------------------------------------------------------- #
# check the reorientation of x-vector in parametric space
#    ckeck node number of nodes in position 23-1 and 25-1 (x-vector) of element
#    table of connectivity
# GR_DEUX
test.assertEqual(con2[0][23 - 1], 38)
test.assertEqual(con2[0][25 - 1], 33)
test.assertEqual(con2[1][23 - 1], 39)
test.assertEqual(con2[1][25 - 1], 38)
# GR_OPPX
#  local x-vector orientation should be inverted
test.assertEqual(con2[2][23 - 1], con1[2][25 - 1])
test.assertEqual(con2[2][25 - 1], con1[2][23 - 1])
# GR_NOORI
#  local x-vector orientation should be the same
test.assertEqual(con2[3][23 - 1], con1[3][23 - 1])
test.assertEqual(con2[3][25 - 1], con1[3][25 - 1])

# ---------------------------------------------------------------------------- #
CA.close()
