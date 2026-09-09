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


PDEPOPG = InputParameter(
    phys=PHY.DEPL_R,
    container="RESU!SAUT_ELGA!N",
    comment="""  PDEPOPG : SAUTS DE DEPLACEMENT AUX POINTS DE GAUSS """,
)

PDEPSNO = OutputParameter(
    phys=PHY.DEPL_R, type="ELNO", comment="""Saut de déplacements des elems aux noeuds"""
)


SAUT_ELNO = Option(
    para_in=(PDEPOPG,),
    para_out=(PDEPSNO,),
    condition=(
        # CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.TYPMOD2, "INTSOLPI"))),
    ),
    comment="""  SAUT_ELNO : SAUT DE DEPLACEMENT PAR ELEM AUX NDS POUR MODELE INTERF_POU  """,
)
