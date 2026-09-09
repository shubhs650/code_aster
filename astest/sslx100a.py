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

import numpy as np


CA.init("--test", "--continue")

test = CA.TestCase()

# unittests for ComponentOnCells objects with masked values
sief = DEP2.getField("SIEF_ELGA", 1)
sief_s = sief.toSimpleFieldOnCells()
sixx = sief_s.SIXX
siyy = sief_s.SIYY
vy = sief_s.VY
mfz = sief_s.MFZ

for norme in ("NORM_1", "NORM_2", "NORM_INFINITY"):
    test.assertAlmostEqual(sief.norm(norme), sief_s.norm(norme), 8, msg="norm")

cmps = sief_s.getComponents()
test.assertEqual(len(cmps), 12, msg="cmps")

test.assertAlmostEqual(vy.min(), -1.0, 8, msg="min")
test.assertAlmostEqual(vy.max(), -1.0, 8, msg="max")
test.assertAlmostEqual(vy.mean(), -1.0, 8, msg="mean")

no_error = mfz / vy
zero = abs(no_error + mfz).max()
test.assertAlmostEqual(zero, 0.0, 8, msg="no error division")

values = vy.values
test.assertAlmostEqual(vy.sum(), -6.0, 8, msg="sum: -1.0 * 6")

m_one = np.ones(values.shape, dtype=float) * -1.0
test.assertAlmostEqual(m_one.sum(), -69.0, 8, msg="-1.0 * 69")

vy.values = m_one
test.assertAlmostEqual(vy._values.sum(), -6.0, 8, msg="sum: only 6 values assigned, other 0.0")
test.assertAlmostEqual(vy.sum(), -6.0, 8, msg="sum")

coor = CALC_CHAM_ELEM(MODELE=MO, OPTION="COOR_ELGA", CARA_ELEM=CARA)
coor_s = coor.toSimpleFieldOnCells()
weight1 = coor_s.W

post = CALC_CHAMP(RESULTAT=DEP2, GROUP_MA="GROUP_MA_GRMA1", CRITERES="SIEQ_ELGA")
sieq = post.getField("SIEQ_ELGA", 1).toSimpleFieldOnCells()
vmis1 = sieq.VMIS

coor = CALC_CHAM_ELEM(MODELE=MO, GROUP_MA="GROUP_MA_GRMA1", OPTION="COOR_ELGA", CARA_ELEM=CARA)
coor_s = coor.toSimpleFieldOnCells()
weight2 = coor_s.W

post = CALC_CHAMP(RESULTAT=DEP2, GROUP_MA="Volume1", CRITERES="SIEQ_ELGA")
sieq = post.getField("SIEQ_ELGA", 1).toSimpleFieldOnCells()
vmis2 = sieq.VMIS

# weight1 has more points than vmis1
wr1 = weight1.onSupportOf(vmis1)
test.assertAlmostEqual(wr1.sum(), 30.0, 6, msg="wr1 sum")

# weight2 has less points than vmis2 (and not defined there)
wr2 = weight2.onSupportOf(vmis2)
test.assertAlmostEqual(wr2.sum(), 0.0, 6, msg="wr2 sum")

wr1p = weight1.onSupportOf(vmis1, strict=True)
test.assertAlmostEqual((wr1 - wr1p).sum(), 0.0, 6, msg="wr1 == wr1p")

with test.assertRaises(IndexError):
    weight2.onSupportOf(vmis2, strict=True)

test.printSummary()

CA.close()
