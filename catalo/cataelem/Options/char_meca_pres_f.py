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

PCHHOBS = InputParameter(phys=PHY.N3600R, comment=""" HHO - coefficient base locale""")


PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, container="CARA!.CANBSP", comment="""  PNBSP_I :  NOMBRE DE SOUS_POINTS  """
)


PCAORIE = InputParameter(
    phys=PHY.CAORIE_R,
    container="ORIE!.CARORIEN",
    comment="""  PCAORIE : ORIENTATION LOCALE D'UN ELEMENT DE POUTRE OU DE TUYAU  """,
)


PPINTTO = InputParameter(phys=PHY.N132_R)


PCNSETO = InputParameter(
    phys=PHY.N1280I,
    container="MODL!.TOPOSE.CNS",
    comment="""  XFEM - CONNECTIVITE DES SOUS-ELEMENTS  """,
)


PHEAVTO = InputParameter(phys=PHY.N512_I)


PLONCHA = InputParameter(
    phys=PHY.N120_I,
    container="MODL!.TOPOSE.LON",
    comment="""  XFEM - NBRE DE TETRAEDRES ET DE SOUS-ELEMENTS  """,
)

PBASLOR = InputParameter(phys=PHY.NEUT_R)

PLSN = InputParameter(phys=PHY.NEUT_R)


PLST = InputParameter(phys=PHY.NEUT_R)


PPINTER = InputParameter(phys=PHY.N816_R)


PAINTER = InputParameter(phys=PHY.N1360R)


PCFACE = InputParameter(phys=PHY.N720_I)


PLONGCO = InputParameter(phys=PHY.N120_I)


PSTANO = InputParameter(phys=PHY.N120_I)


PPMILTO = InputParameter(phys=PHY.N792_R)


PBASECO = InputParameter(phys=PHY.N2448R)


PFISNO = InputParameter(phys=PHY.NEUT_I)


PHEA_FA = InputParameter(phys=PHY.N240_I)


PHEA_NO = InputParameter(phys=PHY.N120_I)


PHEA_SE = InputParameter(phys=PHY.N512_I)


CHAR_MECA_PRES_F = Option(
    para_in=(
        SP.PABSCUR,
        PAINTER,
        PBASECO,
        SP.PCACOQU,
        SP.PCAGEPO,
        PCAORIE,
        PCFACE,
        PCNSETO,
        PFISNO,
        SP.PGEOMER,
        PHEAVTO,
        PHEA_FA,
        PHEA_NO,
        PHEA_SE,
        PLONCHA,
        PLONGCO,
        PLSN,
        PLST,
        PNBSP_I,
        PPINTER,
        PPINTTO,
        PPMILTO,
        SP.PPRESSF,
        PSTANO,
        SP.PINSTR,
        SP.PMATERC,
        PBASLOR,
        PCHHOBS,
    ),
    para_out=(SP.PVECTUR,),
    condition=(
        #     Pour le bord des elements massifs
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "-1"))),
        #     Pour les elements XFEM massifs
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"), (AT.LXFEM, "OUI"))),
        #     Pour les elements coque/plaque ...
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"), (AT.COQUE, "OUI"))),
        #     ... mais pas pour leurs aretes!
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.BORD, "-1"), (AT.COQUE, "OUI"))),
        #     Pour les elements tuyau
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"), (AT.TUYAU, "OUI"))),
        #     Pour les elements solid-shell, la contribution de la pression est volumique (pinching)
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.MODELI, "SSH"))),
        #     Pour les elements membrane
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.MODELI, "MMB"))),
        #     Pour les elements fluide, mais contribution nulle (te0099)
        CondCalcul("+", ((AT.FLUIDE, "OUI"), (AT.BORD, "-1"))),
        #     Pour les elements IFS
        CondCalcul("+", ((AT.FSI, "OUI"),)),
    ),
    comment=""" Second membre pour une pression fonction """,
)
