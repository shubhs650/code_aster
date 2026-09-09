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
# -------------------------------------------------------------------

################# Modification du test zzzz413b ####################


from code_aster.Commands import *

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="NON"))

MAIL = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH")

MODELE = AFFE_MODELE(
    MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="TUYAU_3M")
)

CAREL = AFFE_CARA_ELEM(
    MODELE=MODELE,
    POUTRE=(
        _F(
            SECTION="CERCLE",
            GROUP_MA="TOUT",
            CARA=("R", "EP"),
            VALE=(0.10, 0.05),
            TUYAU_NSEC=8,
            TUYAU_NCOU=3,
        ),
    ),
)

MAT = DEFI_MATERIAU(ELAS=_F(E=2.1e11, NU=0.3, RHO=7800.0))

CHMAT = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MAT))

CHAR = AFFE_CHAR_MECA(
    MODELE=MODELE, FORCE_NODALE=_F(GROUP_NO="CHG", FX=-1000.0, FY=5000.0, FZ=10000.0)
)
CHAR.debugPrint()

BLOCAGE = AFFE_CHAR_CINE(
    MECA_IMPO=_F(DX=0.0, DY=0.0, DZ=0.0, DRX=0.0, DRY=0.0, DRZ=0.0, GROUP_NO=("ENC",)),
    MODELE=MODELE,
)

L_INIT = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1, NOMBRE=10))

FOMULT = DEFI_FONCTION(
    NOM_PARA="INST", VALE=(0.0, 0.0, 1.0, 1.0), PROL_DROITE="EXCLU", PROL_GAUCHE="EXCLU"
)

RES = STAT_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    CARA_ELEM=CAREL,
    EXCIT=(_F(CHARGE=BLOCAGE, FONC_MULT=FOMULT), _F(CHARGE=CHAR, FONC_MULT=FOMULT)),
    COMPORTEMENT=(_F(TOUT="OUI", RELATION="ELAS", DEFORMATION="PETIT"),),
    INCREMENT=_F(LIST_INST=L_INIT),
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-8, ITER_GLOB_ELAS=25, ITER_GLOB_MAXI=10),
    NEWTON=_F(PREDICTION="ELASTIQUE", REAC_INCR=1, MATRICE="TANGENTE", REAC_ITER=3),
    SOLVEUR=_F(METHODE="MUMPS"),
)

RES_N = MECA_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    CARA_ELEM=CAREL,
    EXCIT=(_F(CHARGE=BLOCAGE, FONC_MULT=FOMULT), _F(CHARGE=CHAR, FONC_MULT=FOMULT)),
    COMPORTEMENT=(_F(TOUT="OUI", RELATION="ELAS", DEFORMATION="PETIT"),),
    INCREMENT=_F(LIST_INST=L_INIT),
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-8, ITER_GLOB_ELAS=25, ITER_GLOB_MAXI=10),
    NEWTON=_F(PREDICTION="ELASTIQUE", REAC_INCR=1, MATRICE="TANGENTE", REAC_ITER=3),
    SOLVEUR=_F(METHODE="MUMPS"),
)

nbIndexes = RES.getNumberOfIndexes()

DEPL_REF = RES.getField("DEPL", nbIndexes - 1)
SIGMA_REF = RES.getField("SIEF_ELGA", nbIndexes - 1)
VARI_REF = RES.getField("VARI_ELGA", nbIndexes - 1)

DEPL = RES_N.getField("DEPL", nbIndexes - 1)
SIGMA = RES_N.getField("SIEF_ELGA", nbIndexes - 1)
VARI = RES_N.getField("VARI_ELGA", nbIndexes - 1)

DIF_DEPL = DEPL_REF - DEPL
DIF_SIG = SIGMA_REF - SIGMA
DIF_VAR = VARI_REF - VARI

TEST_RESU(
    CHAM_ELEM=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=100.0,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_SIG,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=100.0,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_SIG,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=1.0,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=1.0,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

TEST_RESU(
    CHAM_NO=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=10e-8,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_DEPL,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=10e-8,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_DEPL,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

FIN()
