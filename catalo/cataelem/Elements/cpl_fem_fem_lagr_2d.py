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


# ELEMENTARY TREATMENT OF 2D COUPLING FEM-FEM
# AUGMENTED LAGRANGIAN FORMULATION

# Convention:
# EN1: Slave nodes - FEM
# EN2: Master nodes - FEM

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        # Slave nodes - FEM
        ("EN1", ("DX", "DY", "LAGX", "LAGY")),
        # Master nodes - FEM
        ("EN2", ("DX", "DY")),
    ),
)


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class CL_S2S3(Element):
    """
    THE CL_S2S3 CLASS ELEMENT : SEG2/SEG3
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 2D : elementary treatments
    """

    meshType = MT.SEG23
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5)))
    calculs = (
        OP.CHAR_MECA_CPL(
            te=523,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PPAIRR, LC.CPAIRR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_CPL(
            te=523,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PPAIRR, LC.CPAIRR)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_ELAS_CPL(
            te=523,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PPAIRR, LC.CPAIRR)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
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
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
    )


# ------------------------------------------------------------
class CL_S3S3(CL_S2S3):
    """
    CL_S3S3 DERIVED FROM THE CL_S2S3 CLASS ELEMENT : SEG3/SEG3
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 2D : elementary treatments
    """

    meshType = MT.SEG33
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6)))


# ------------------------------------------------------------
class CL_S3S2(CL_S2S3):
    """
    CL_S3S3 DERIVED FROM THE CL_S2S3 CLASS ELEMENT : SEG3/SEG2
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 2D : elementary treatments
    """

    meshType = MT.SEG32
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5)))


# ------------------------------------------------------------
class CL_S2S2(CL_S2S3):
    """
    CL_S3S3 DERIVED FROM THE CL_S2S2 CLASS ELEMENT : SEG2/SEG2
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 2D : elementary treatments
    """

    meshType = MT.SEG22
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4)))


# ------------------------------------------------------------
class CL_S4S4(CL_S2S3):
    """
    CL_S3S3 DERIVED FROM THE CL_S2S2 CLASS ELEMENT : SEG2/SEG2
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 2D : elementary treatments
    """

    meshType = MT.SEG44
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7, 8)))


# ------------------------------------------------------------
class CL_S4S3(CL_S2S3):
    """
    CL_S3S3 DERIVED FROM THE CL_S2S2 CLASS ELEMENT : SEG2/SEG2
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 2D : elementary treatments
    """

    meshType = MT.SEG43
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7)))


# ------------------------------------------------------------
class CL_POI2D(CL_S2S3):
    """
    CL_S3S3 DERIVED FROM THE CL_S2S2 CLASS ELEMENT : SEG2/SEG2
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 2D : elementary treatments
    """

    meshType = MT.POI1
    nodes = (SetOfNodes("EN1", (1,)), SetOfNodes("EN2", ()))
    calculs = (
        OP.CHAR_MECA_CPL(
            te=524,
            para_in=((SP.PDEPLMR, DDL_MECA), (SP.PDEPLPR, DDL_MECA), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_CPL(te=524, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PMATUUR, MMATUUR),)),
        OP.RIGI_ELAS_CPL(
            te=524, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PMATUUR, MMATUUR),)
        ),
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
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
    )
