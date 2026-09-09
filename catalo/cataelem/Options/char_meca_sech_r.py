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

PVARCPR = InputParameter(phys=PHY.VARI_R, comment="""  PVARCPR : VARIABLES DE COMMANDE  """)

PCAORIE = InputParameter(
    phys=PHY.CAORIE_R,
    container="ORIE!.CARORIEN",
    comment="""  PCAORIE : ORIENTATION LOCALE D'UN ELEMENT DE POUTRE, TUYAU ...  """,
)

PCOMPOR = InputParameter(phys=PHY.COMPOR)

PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, container="CARA!.CANBSP", comment="""  PNBSP_I : NOMBRE DE SOUS_POINTS  """
)


# For HHO
PCHHOGT = InputParameter(phys=PHY.N1920R, comment=""" HHO - matrice du gradient local""")
PCHHOBS = InputParameter(phys=PHY.N3600R, comment=""" HHO - coefficient base locale""")


CHAR_MECA_SECH_R = Option(
    para_in=(
        SP.PCAGNBA,
        SP.PCAGNPO,
        SP.PCAMASS,
        PCAORIE,
        PCOMPOR,
        PNBSP_I,
        SP.PFIBRES,
        SP.PGEOMER,
        SP.PMATERC,
        SP.PINSTR,
        PVARCPR,
        SP.PVARCRR,
        PCHHOGT,
        PCHHOBS,
    ),
    para_out=(SP.PVECTUR,),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        #     Les elements d'interface ne sont pas concernes (issue24099) :
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.INTERFACE, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.FLUIDE, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.MODELI, "3FI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.MODELI, "AFI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.MODELI, "PFI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.MODELI, "D2D"))),
    ),
    comment=""" CHAR_MECA_SECH_R: (MOT-CLE : SECH_CALCULEE): CALCUL DU SECOND
           MEMBRE ELEMENTAIRE CORRESPONDANT AU CHAMP DE SECHAGE """,
)
