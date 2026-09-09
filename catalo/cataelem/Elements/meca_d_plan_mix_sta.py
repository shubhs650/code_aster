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

NDEPLAR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))

DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "EPXX", "EPYY", "EPZZ", "EPXY")
)

ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIPXX", "SIPYY", "SIPZZ", "SIPXY"),
)


ECONTNO = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELNO",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIPXX", "SIPYY", "SIPZZ", "SIPXY"),
)

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MVECTDR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NDEPLAR)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MEMS_QU4(Element):
    """QUAD4 modélisation D_PLAN_MIX_STA formulation STA"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)),)
    elrefe = (
        ElrefeLoc(
            MT.QU4,
            gauss=("RIGI=FPG4", "FPG1=FPG1", "MASS=FPG4", "NOEU=NOEU"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )
    calculs = (
        OP.ADD_SIGM(
            te=581,
            para_in=((SP.PEPCON1, ECONTPG), (SP.PEPCON2, ECONTPG)),
            para_out=((SP.PEPCON3, ECONTPG),),
        ),
        OP.CHAR_MECA_EPSA_R(
            te=421,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_EPSA_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_EPSI_F(
            te=284,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PEPSINF, LC.CEPS2DF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_EPSI_F.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_EPSI_R(
            te=284,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PEPSINR, LC.EGPS2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_FF2D2D(
            te=94,
            para_in=((SP.PFF2D2D, LC.CFOR2DF), (SP.PGEOMER, LC.EGEOM2D), (SP.PINSTR, LC.MTEMPSR)),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_FR2D2D(
            te=93,
            para_in=((SP.PFR2D2D, LC.NFOR2DR), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_META_Z(
            te=353,
            para_in=(
                (OP.CHAR_MECA_META_Z.PCOMPOR, LC.CCOMPOR),
                (OP.CHAR_MECA_META_Z.PCONTMR, ECONTPG),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.CHAR_MECA_META_Z.PVARCPR, LC.ZVARCPG),
                (OP.CHAR_MECA_META_Z.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=85,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_ROTA_R(
            te=84,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PROTATR, LC.CROTATR),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_SECH_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_SECH_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.COOR_ELGA(
            te=479,
            para_in=((SP.PGEOMER, LC.EGEOM2D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU2D),),
        ),
        OP.EPEQ_ELGA(
            te=335,
            para_in=((OP.EPEQ_ELGA.PDEFORR, LC.EGPS2DR),),
            para_out=((OP.EPEQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPEQ_ELNO(
            te=335,
            para_in=((OP.EPEQ_ELNO.PDEFORR, LC.EEPS2DR),),
            para_out=((OP.EPEQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPGQ_ELGA(
            te=335,
            para_in=((OP.EPGQ_ELGA.PDEFORR, LC.EGPS2DR),),
            para_out=((OP.EPGQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPGQ_ELNO(
            te=335,
            para_in=((OP.EPGQ_ELNO.PDEFORR, LC.EEPS2DR),),
            para_out=((OP.EPGQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPFD_ELGA(
            te=528,
            para_in=(
                (OP.EPFD_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PVARIGR, LC.ZVARIPG),
            ),
            para_out=((OP.EPFD_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPFD_ELNO(
            te=4,
            para_in=((OP.EPFD_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPFP_ELGA(
            te=528,
            para_in=(
                (OP.EPFP_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPFP_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIGR, LC.ZVARIPG),
            ),
            para_out=((OP.EPFP_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPFP_ELNO(
            te=4,
            para_in=((OP.EPFP_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPME_ELGA(
            te=87,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (OP.EPME_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPME_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPME_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPME_ELNO(
            te=4,
            para_in=((OP.EPME_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPSG_ELGA(
            te=87,
            para_in=(
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPSG_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSG_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPSG_ELNO(
            te=4,
            para_in=((OP.EPSG_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPSI_ELGA(
            te=87,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPSI_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PDEFOPC, LC.EGPS2DC), (OP.EPSI_ELGA.PDEFOPG, LC.EGPS2DR)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONC, LC.EEPS2DC), (SP.PDEFONO, LC.EEPS2DR)),
        ),
        OP.EPVC_ELGA(
            te=529,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPVC_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPVC_ELGA.PDEFOPG, LC.EGVARC3D),),
        ),
        OP.EPVC_ELNO(
            te=4,
            para_in=((OP.EPVC_ELNO.PDEFOPG, LC.EGVARC3D),),
            para_out=((SP.PDEFONO, LC.NVARC3D),),
        ),
        OP.FORC_NODA(
            te=54,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=54,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
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
        OP.FULL_MECA_ELAS(
            te=54,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA_ELAS.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.FULL_MECA_ELAS.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA_ELAS.PCONTPR, ECONTPG),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA_ELAS.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.PAS_COURANT(
            te=404,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.PAS_COURANT.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.RAPH_MECA(
            te=54,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.RAPH_MECA.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, ECONTPG),
                (OP.RAPH_MECA.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=54,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PRESIREF, LC.CRESIREF),
                (SP.PRESICMP, LC.CRESICMP),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.REPERE_LOCAL(
            te=133,
            para_in=((SP.PCAMASS, LC.CCAMASS), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((SP.PREPLO1, LC.CGEOM2D), (SP.PREPLO2, LC.CGEOM2D), (SP.PREPLO3, LC.CGEOM2D)),
        ),
        OP.RIGI_MECA_ELAS(
            te=54,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_ELAS.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARIMR, LC.ZVARIPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_MECA_TANG(
            te=54,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
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
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, ECONTPG), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, LC.ESIG2DC), (OP.SIEF_ELNO.PSIEFNOR, ECONTNO)),
        ),
        OP.SIEQ_ELGA(
            te=335,
            para_in=((OP.SIEQ_ELGA.PCONTRR, LC.EGIG2DR),),
            para_out=((OP.SIEQ_ELGA.PCONTEQ, LC.ECOEQPG),),
        ),
        OP.SIEQ_ELNO(
            te=335,
            para_in=((OP.SIEQ_ELNO.PCONTRR, LC.ESIG2DR),),
            para_out=((OP.SIEQ_ELNO.PCONTEQ, LC.ECOEQNO),),
        ),
        OP.SIGM_ELGA(
            te=546,
            para_in=((SP.PSIEFR, LC.EGIG2DR),),
            para_out=((SP.PSIGMC, LC.EGIG2DC), (SP.PSIGMR, LC.EGIG2DR)),
        ),
        OP.SIGM_ELNO(
            te=4,
            para_in=((OP.SIGM_ELNO.PCONTRR, LC.EGIG2DR),),
            para_out=((SP.PSIEFNOC, LC.ESIG2DC), (OP.SIGM_ELNO.PSIEFNOR, LC.ESIG2DR)),
        ),
        OP.SIMY_ELGA(
            te=6,
            para_in=((OP.SIMY_ELGA.PCONTRR, LC.EGIG2DR), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((OP.SIMY_ELGA.PSIEFNOR, LC.EGIG2DR),),
        ),
        OP.SING_ELEM(te=99, para_out=((SP.PSING_R, LC.ESINGUL),)),
        OP.SING_ELNO(te=99, para_out=((SP.PSINGNO, LC.ESINGNO),)),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
                (SP.PNEUT_F, LC.CNTINIF),
                (SP.PNEUT_R, LC.CNTINIR),
                (SP.PNEUT_I, LC.CNTINII),
            ),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PDEPL_R, LC.EGDEP2D),
                (OP.TOU_INI_ELGA.PDOMMAG, LC.EDOMGGA),
                (OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO2D),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (SP.PVALO_R, LC.EGTINIV),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
                (OP.TOU_INI_ELGA.PNEUT_I, LC.EGTINII),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PDOMMAG, LC.EDOMGNO),
                (OP.TOU_INI_ELNO.PEPSI_R, LC.EEPS2DR),
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (SP.PTEMPN_R, LC.ETEMPNO),
                (OP.TOU_INI_ELNO.PSIEF_R, ECONTNO),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
        OP.VARC_ELGA(
            te=530,
            para_in=((OP.VARC_ELGA.PVARCPR, LC.ZVARCPG),),
            para_out=((SP.PVARC_R, LC.EVARC_R),),
        ),
        OP.VARC_ELNO(
            te=4, para_in=((SP.PVARCGR, LC.EVARC_R),), para_out=((SP.PVARCNR, LC.EVARCNR),)
        ),
        OP.VARI_ELNO(
            te=4,
            para_in=((SP.PVARIGR, LC.ZVARIPG),),
            para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------
class MEMS_TR3(MEMS_QU4):
    """TR3 modélisation D_PLAN_MIX_STA formulation STA"""

    meshType = MT.TRIA3
    nodes = SetOfNodes("EN1", (1, 2, 3))
    elrefe = (
        ElrefeLoc(
            MT.TR3,
            gauss=("RIGI=FPG3", "FPG1=FPG1", "MASS=FPG3", "NOEU=NOEU"),
            mater=("RIGI", "MASS", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )


# ------------------------------------------------------------
class MEMI_QU4(MEMS_QU4):
    """QUAD4 modélisation D_PLAN_MIX_STA formulation STA_INCO"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)),)
    elrefe = (
        ElrefeLoc(
            MT.QU4,
            gauss=("RIGI=FPG4", "FPG1=FPG1", "MASS=FPG4", "NOEU=NOEU"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )
