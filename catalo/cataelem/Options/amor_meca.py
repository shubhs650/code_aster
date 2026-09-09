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


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT

PCOMPOR = InputParameter(phys=PHY.COMPOR)

PVARCPR = InputParameter(phys=PHY.VARI_R)

PAMORFL = InputParameter(phys=PHY.NEUT_I)


PCAORIE = InputParameter(
    phys=PHY.CAORIE_R,
    container="ORIE!.CARORIEN",
    comment="""  PCAORIE : ORIENTATION LOCALE D'UN ELEMENT DE POUTRE OU DE TUYAU  """,
)


AMOR_MECA = Option(
    para_in=(
        SP.PCADISA,
        PCAORIE,
        SP.PCINFDI,
        PCOMPOR,
        SP.PGEOMER,
        SP.PMASSEL,
        SP.PMATERC,
        SP.PRIGIEL,
        SP.PRIGINS,
        PVARCPR,
        SP.PVARCRR,
        SP.PNONLIN,
        SP.PVARIPG,
        PAMORFL,
    ),
    para_out=(SP.PMATUNS, SP.PMATUUR),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        CondCalcul("-", ((AT.FLUIDE, "OUI"), (AT.FORMULATION, "U_P"))),
        CondCalcul("-", ((AT.FLUIDE, "OUI"), (AT.FORMULATION, "U_P_PHI"))),
        CondCalcul("+", ((AT.FLUIDE, "OUI"), (AT.ABSO, "OUI"))),
        CondCalcul("+", ((AT.FSI, "OUI"), (AT.BORD, "-1"), (AT.FORMULATION, "U_PSI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.MODELI, "D2D"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.MODELI, "D3D"))),
        # CondCalcul("-", ((AT.PHENO, "ME"), (AT.TYPMOD2, "INTERFAC"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.INTERFACE, "OUI"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.TYPMOD2, "ELEMJOIN"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.TYPMOD2, "EJ_HYME"))),
    ),
)
