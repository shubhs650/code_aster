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


import sys
import unittest
from glob import glob
from pathlib import Path

# Note: "--hide-command" is used in the '.export' file.
from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities.logger import WARNING, logger
from run_aster.export import Export

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))
test = CA.TestCase()

logger.setLevel(WARNING)


def _test_module(module):
    print(f"\n\n+++ testing {module}...\n", flush=True)
    result = unittest.main(argv=["comm"], module=module, exit=False, verbosity=2).result
    # to flush printings from unittest
    sys.stdout.flush()
    sys.stderr.flush()
    test.assertTrue(result.wasSuccessful())


# build Export object
lfexp = glob("*.export")
assert len(lfexp) == 1, lfexp
export = Export(lfexp[0])

for input in export.datafiles:
    if input.filetype != "nom":
        continue
    mod = Path(input.path).stem
    _test_module(mod)


test.printSummary()

CA.close()
