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

from cataelem.Tools.base_objects import LocatedComponents, ArrayOfComponents, SetOfNodes, ElrefeLoc
from cataelem.Tools.base_objects import Calcul, Element
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.located_components as LC
import cataelem.Commons.parameters as SP
import cataelem.Commons.mesh_types as MT
from cataelem.Options.options import OP

# ----------------------------------------------------------------------------------------------
# Located components
# ----------------------------------------------------------------------------------------------

DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        ("EN1", ("DX", "DY", "DZ")),
        ("EN2", ("DX", "DY", "DZ", "DRX", "DRY", "DRZ")),
        ("EN3", ()),
    ),
)

EDEPSPG = LocatedComponents(
    phys=PHY.DEPL_R, type="ELGA", location="RIGI", components=("DX", "DY", "DZ")
)

EDEPSNO = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))

CCAGEPO = LocatedComponents(phys=PHY.CAGEPO_R, type="ELEM", components=("HY1", "HZ1", "R1", "TSEC"))

CCAORIE = LocatedComponents(phys=PHY.CAORIE_R, type="ELEM", components=("ALPHA", "BETA", "GAMMA"))

EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)

ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))

NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))

EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)

CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))

EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))

EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))

ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R, type="ELGA", location="RIGI", components=("FX", "FY", "FZ")
)

ECONTNO = LocatedComponents(phys=PHY.SIEF_R, type="ELNO", components=("FX", "FY", "FZ"))

ECONTNC = LocatedComponents(phys=PHY.SIEF_C, type="ELNO", components=("FX", "FY", "FZ"))

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ----------------------------------------------------------------------------------------------
class MECA_S8P2(Element):
    """Interface element for soil-pile interaction based on 8-nodes support for the soil"""

    meshType = MT.HEXA27
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
        SetOfNodes("EN2", (23, 25)),
        SetOfNodes("EN3", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 26, 27)),
    )
    elrefe = (
        # ElrefeLoc(MT.H10, gauss=("RIGI=FIS3",), mater=("RIGI",)),
        # ElrefeLoc(MT.H10, gauss=("RIGI=FIS4",), mater=("RIGI",)),
        ElrefeLoc(MT.H10, gauss=("RIGI=FIS2",), mater=("RIGI",)),
    )
    calculs = (
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU3D),)
        ),
        OP.SAUT_ELGA(
            te=383,
            para_in=(
                (SP.PDEPLAR, DDL_MECA),
                (SP.PCAGNPO, LC.CCAGNP1),
                (SP.PCAORIE, CCAORIE),
                (OP.SAUT_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PDEPSGA, EDEPSPG),),
        ),
        OP.SAUT_ELNO(
            te=4,
            para_in=((OP.SAUT_ELNO.PDEPOPG, EDEPSPG),),
            para_out=((OP.SAUT_ELNO.PDEPSNO, EDEPSNO),),
        ),
        OP.SIEF_ELGA(
            te=363,
            para_in=(
                (SP.PDEPLAR, DDL_MECA),
                (SP.PCAGNPO, LC.CCAGNP1),
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.SIEF_ELGA.PCONTRR, ECONTPG),),
        ),
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, ECONTPG), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, ECONTNC), (OP.SIEF_ELNO.PSIEFNOR, ECONTNO)),
        ),
        OP.VARI_ELNO(
            te=4,
            para_in=((SP.PVARIGR, LC.ZVARIPG),),
            para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),),
        ),
        OP.FORC_NODA(
            te=376,
            para_in=(
                (SP.PCAGNPO, LC.CCAGNP1),
                (OP.FORC_NODA.PCAORIE, CCAORIE),
                (OP.FORC_NODA.PCOMPOR, LC.CCOMPOR),
                (OP.FORC_NODA.PSIEFR, ECONTPG),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=366,
            para_in=(
                (SP.PCAGNPO, LC.CCAGNP1),
                (SP.PCAGEPO, CCAGEPO),
                (OP.FULL_MECA.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.FULL_MECA.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, ECONTPG),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.RIGI_MECA_TANG(
            te=366,
            para_in=(
                (SP.PCAGNPO, LC.CCAGNP1),
                (SP.PCAGEPO, CCAGEPO),
                (OP.RIGI_MECA_TANG.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, ECONTPG),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.RIGI_MECA(
            te=367,
            para_in=(
                (SP.PCAGNPO, LC.CCAGNP1),
                (SP.PCAGEPO, CCAGEPO),
                (OP.RIGI_MECA.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PNEUT_I, LC.EGTINII),
            ),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
                (SP.PNEUT_F, LC.CNTINIF),
                (SP.PNEUT_R, LC.CNTINIR),
                (SP.PNEUT_I, LC.CNTINII),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, ENGEOM_R),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PSIEF_R, ECONTNO),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
    )
