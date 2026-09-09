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

from cataelem.Tools.base_objects import SetOfNodes, ElrefeLoc
from cataelem.Elements.ther_plan_hho444 import THER2DQ9_HHO444
import cataelem.Commons.mesh_types as MT


# --------------------------------------------------------------------------------------------------
class THERAXQ9_HHO444(THER2DQ9_HHO444):
    """Thermics - HHO_QUAD - PLAN - QUAD"""

    meshType = MT.QUAD9
    nodes = (
        SetOfNodes("EN1", (5, 6, 7, 8)),
        SetOfNodes("EN2", (1, 2, 3, 4)),
        SetOfNodes("EN3", (9,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.QU9, gauss=("RIGI=FPG25", "FPG1=FPG1", "MASS=FPG25"), mater=("RIGI", "FPG1", "MASS")
        ),
    )


# --------------------------------------------------------------------------------------------------
class THERAXT7_HHO444(THERAXQ9_HHO444):
    """Thermics - HHO_QUAD - PLAN - TRIA"""

    meshType = MT.TRIA7
    nodes = (SetOfNodes("EN1", (4, 5, 6)), SetOfNodes("EN2", (1, 2, 3)), SetOfNodes("EN3", (7,)))
    elrefe = (
        ElrefeLoc(
            MT.TR7, gauss=("RIGI=FPG19", "FPG1=FPG1", "MASS=FPG19"), mater=("RIGI", "FPG1", "MASS")
        ),
    )
