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

import numpy as np

from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import config

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

mesh.setGroupOfCells("OBSERV", [149, 150])
mesh.setGroupOfNodes("OBSERV", [100, 101])

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

FTRACTUB = DEFI_FONCTION(
    NOM_PARA="EPSI",
    VALE=(1000.0, 2000.0, 2000.0, 5000.0),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3), ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=200.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))


encast = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="BAS", DX=0, DY=0.0, DZ=0.0),))

depl = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="HAUT", DZ=1.0),))

LIST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=2))

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1000.0, 1000.0))

NORME = FORMULE(VALE="sqrt(DX*DX+DY*DY+DZ*DZ)", NOM_PARA=["DX", "DY", "DZ"])
MULT = FORMULE(VALE="V1*V2", NOM_PARA=["V1", "V2"])

common_keywords = _F(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="GDEF_LOG"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-8),
    INCREMENT=_F(LIST_INST=LIST),
    OBSERVATION=(
        _F(NOM_CMP=("V1", "V2"), NOM_CHAM="VARI_ELGA", GROUP_MA="OBSERV", EVAL_ELGA="MAX"),
        _F(
            NOM_VARI=("EPSPEQ", "INDIPLAS"),
            NOM_CHAM="VARI_ELGA",
            GROUP_MA="OBSERV",
            EVAL_ELGA="MAX",
            INST=0.5,
        ),
        _F(NOM_CMP="V1", NOM_CHAM="VARI_ELGA", GROUP_MA="OBSERV", POINT=(1, 2)),
        _F(
            NOM_CMP=("V1", "V2"),
            NOM_CHAM="VARI_ELGA",
            GROUP_MA="OBSERV",
            POINT=(1, 2),
            EVAL_CMP="FORMULE",
            FORMULE=MULT,
        ),
        _F(NOM_CMP="V1", NOM_CHAM="VARI_ELGA", TOUT="OUI", POINT=1, EVAL_CHAM="MAX"),
        _F(NOM_CMP="V1", NOM_CHAM="VARI_ELGA", TOUT="OUI", POINT=(1, 2), EVAL_CHAM="MAX"),
        _F(NOM_CMP="SIZZ", NOM_CHAM="SIEF_ELGA", MAILLE=("150", "151"), POINT=2, PAS_OBSE=2),
        _F(
            NOM_CMP=("SIXX", "SIYY", "SIZZ"),
            NOM_CHAM="SIEF_ELGA",
            TOUT="OUI",
            EVAL_ELGA="MIN",
            EVAL_CHAM="MINI_ABS",
        ),
        _F(
            NOM_CMP=("DX", "DY"),
            NOM_CHAM="DEPL",
            GROUP_MA="OBSERV",
            LIST_INST=LIST,
            OBSE_ETAT_INIT="NON",
        ),
        _F(NOM_CMP="DZ", NOM_CHAM="DEPL", TOUT="OUI", EVAL_CHAM="MIN"),
        _F(
            NOM_CMP=("DX", "DY", "DZ"),
            EVAL_CMP="FORMULE",
            FORMULE=NORME,
            NOM_CHAM="DEPL",
            NOEUD=("101", "102"),
        ),
        _F(
            NOM_CMP=("DX", "DY", "DZ"),
            EVAL_CMP="FORMULE",
            FORMULE=NORME,
            NOM_CHAM="DEPL",
            GROUP_NO="OBSERV",
            EVAL_CHAM="MAX",
        ),
    ),
    INFO=1,
)

init_keywords = common_keywords.copy()
init_keywords["INCREMENT"] = _F(LIST_INST=LIST, INST_FIN=0.5)

# STAT_NON_LINE DE REFERENCE
ressnl = STAT_NON_LINE(**init_keywords)
ressnl = STAT_NON_LINE(reuse=ressnl, ETAT_INIT=_F(EVOL_NOLI=ressnl), **common_keywords)


resmnl = MECA_NON_LINE(**init_keywords)
resmnl = MECA_NON_LINE(reuse=resmnl, ETAT_INIT=_F(EVOL_NOLI=resmnl), **common_keywords)


def assert_same_table(table_1, table_2):
    for parameter in table_1.getParameters():
        column_1 = table_1.get_column(parameter)
        if all([x is None for x in column_1]):
            continue
        test.assertEqual(table_1.getColumnType(parameter), table_2.getColumnType(parameter))
        column_2 = table_2.get_column(parameter)
        if parameter == "VALE":
            test.assertArrayEqual(
                np.array(column_1), np.array(column_2), rtol=1.0e-10, atol=1.0e-10
            )
        else:
            test.assertEqual(column_1, column_2)


tabsnl = RECU_TABLE(CO=ressnl, NOM_TABLE="OBSERVATION")
tabmnl = RECU_TABLE(CO=resmnl, NOM_TABLE="OBSERVATION")
IMPR_TABLE(TABLE=tabsnl, UNITE=6)
IMPR_TABLE(TABLE=tabmnl, UNITE=6)
assert_same_table(tabsnl, tabmnl)


# =========================================================
#          DETERMINATION DE LA REFERENCE
# =========================================================

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant
SIGMA_REF = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_SIEF_R", NOM_CHAM="SIEF_ELGA", RESULTAT=ressnl, INST=1.0
)

VARI_REF = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_VARI_R", NOM_CHAM="VARI_ELGA", RESULTAT=ressnl, INST=1.0
)

SIGMA1 = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_SIEF_R", NOM_CHAM="SIEF_ELGA", RESULTAT=resmnl, INST=1.0
)

VARI1 = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_VARI_R", NOM_CHAM="VARI_ELGA", RESULTAT=resmnl, INST=1.0
)

# =========================================================
#            REALISATION DES TESTS
# =========================================================

DIF_SIG1 = SIGMA_REF - SIGMA1
DIF_VAR1 = VARI_REF - VARI1

TEST_RESU(
    CHAM_ELEM=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-08,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_SIG1,
            VALE_CALC=1.5063505998114124e-12,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-08,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_SIG1,
            VALE_CALC=1.7053025658242404e-12,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            ORDRE_GRANDEUR=5.0e-03,
            PRECISION=1.0e-08,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_VAR1,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            ORDRE_GRANDEUR=5.0e-3,
            PRECISION=1.0e-08,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_VAR1,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

# =========================================================
#            SOLVEUR NON LINEAIRE SNES
# =========================================================

if config["ASTER_HAVE_PETSC4PY"]:
    myOptions = "-pc_type lu -pc_factor_mat_solver_type mumps -ksp_type fgmres -snes_linesearch_type basic  -snes_max_it 10 -snes_view "
    SOLU2 = MECA_NON_LINE(
        MODELE=model,
        CHAM_MATER=mater,
        EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
        COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="GDEF_LOG"),
        NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
        METHODE="SNES",
        CONVERGENCE=_F(RESI_GLOB_RELA=1e-8),
        INCREMENT=_F(LIST_INST=LIST),
        SOLVEUR=_F(METHODE="PETSC", OPTION_PETSC=myOptions),
        INFO=1,
    )

    SIGMA2 = CREA_CHAMP(
        OPERATION="EXTR", TYPE_CHAM="ELGA_SIEF_R", NOM_CHAM="SIEF_ELGA", RESULTAT=SOLU2, INST=1.0
    )

    VARI2 = CREA_CHAMP(
        OPERATION="EXTR", TYPE_CHAM="ELGA_VARI_R", NOM_CHAM="VARI_ELGA", RESULTAT=SOLU2, INST=1.0
    )

    DIF_SIG2 = SIGMA_REF - SIGMA2
    DIF_VAR2 = VARI_REF - VARI2

    TEST_RESU(
        CHAM_ELEM=(
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-08,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_SIG2,
                VALE_CALC=1.5063505998114124e-12,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                PRECISION=1.0e-08,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_SIG2,
                VALE_CALC=1.7053025658242404e-12,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                ORDRE_GRANDEUR=5.0e-03,
                PRECISION=1.0e-08,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_VAR2,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="AUTRE_ASTER",
                ORDRE_GRANDEUR=5.0e-3,
                PRECISION=1.0e-08,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_VAR2,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
        )
    )


FIN()
