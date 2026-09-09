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

from code_aster.Commands import *
from code_aster import CA

import numpy as np


def testLoi(nombrePas=1, nomComportement="INTERF_POU_ELAS", optionSurfCharge=False):

    # ---------------------------------------------------------------------------- #
    CA.init("--test", "--abort")

    test = CA.TestCase()

    ################################################################################

    #  MAILLAGE
    #

    Grma27 = "Hexa27"

    # vecteur d'orientation de l'axe de la poutre
    vOrie = [1.0, 0.0, 0.0]

    # lire maillage
    MA = CA.Mesh()
    MA.readAsterFile("zzzz366a.mail")

    MA = MODI_MAILLAGE(reuse=MA, MAILLAGE=MA, ORIE_INTERF_POU=_F(GROUP_MA=Grma27, VECT_ORIE=vOrie))

    ################################################################################

    #  CARACTERISTIQUES GEOMETRIQUES
    #

    # pieu
    hy = 1.0
    hz = 1.0
    anglvril = -np.pi / 3

    ################################################################################

    #  CARACTERISTIQUES MATERIAU
    #

    # pieu
    E_p = 1
    nu_p = 0.25

    # interface
    k_n = 5
    k_t = 3
    f_ny = k_n
    f_ty = k_t
    g_n = k_n / 3
    g_t = k_t / 3

    # ATTENTION :
    #  Les valeurs des paramètres matériau n'ont pas de sens physique.
    #  Elles servent uniquement à des fins de vérification.

    ################################################################################

    #  DEFINITION DES MATERIAUX
    #

    # pieu

    POU = DEFI_MATERIAU(ELAS=_F(E=E_p, NU=nu_p))

    # interface

    d_parel = dict(K_N=k_n, K_T=k_t)
    d_parpl = dict(F_NY=f_ny, F_TY=f_ty, G_N=g_n, G_T=g_t)

    d_parel_fo = dict()
    d_parpl_fo = dict()

    for par in [d_parel, d_parpl]:
        for key in par.keys():
            l_vale = (-1.0, par[key], 1.0, par[key])
            par[key] = DEFI_FONCTION(
                NOM_PARA="X", PROL_DROITE="CONSTANT", PROL_GAUCHE="CONSTANT", VALE=l_vale
            )

    INT = DEFI_MATERIAU(INTERF_POU_ELAS_FO=d_parel, INTERF_POU_CINE_FO=d_parpl)

    ################################################################################

    #  AFFECTATIONS AU MODELE
    #

    MODELE = AFFE_MODELE(
        MAILLAGE=MA,
        AFFE=(_F(GROUP_MA="Hexa27", PHENOMENE="MECANIQUE", MODELISATION="3D_INTERF_POU"),),
    )

    CAREL = AFFE_CARA_ELEM(
        MODELE=MODELE,
        POUTRE=(_F(GROUP_MA="Hexa27", SECTION="RECTANGLE", CARA=("HY", "HZ"), VALE=(hy, hz)),),
        ORIENTATION=(_F(GROUP_MA="Hexa27", CARA="ANGL_VRIL", VALE=np.rad2deg(anglvril)),),
    )

    CHMAT = AFFE_MATERIAU(MODELE=MODELE, AFFE=(_F(GROUP_MA="Hexa27", MATER=(INT, POU)),))

    ############################################

    #  CHARGEMENT
    #

    l_inst = [0.0, 1.0, 2.0, 3.0, 4.0]
    l_valy = [0.0, 1.0, 1.25, -0.75, -1.0]
    l_inst2 = l_inst[: nombrePas + 1]
    l_valy2 = l_valy[: nombrePas + 1]

    FONCT = DEFI_FONCTION(
        NOM_PARA="INST", VALE=[item for pair in zip(l_inst2, l_valy2) for item in pair]
    )

    LINST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=nombrePas, NOMBRE=nombrePas),))

    ################################################################################

    #  CONDITIONS LIMITES
    #

    CLF = AFFE_CHAR_MECA(MODELE=MODELE, FORCE_NODALE=(_F(GROUP_NO=["N23", "N25"], FY=-f_ny),))

    CLD = AFFE_CHAR_CINE(
        MODELE=MODELE,
        MECA_IMPO=(
            _F(GROUP_NO=["N23", "N25"], DX=0.0, DZ=0.0, DRX=0.0, DRY=0.0, DRZ=0.0),
            _F(GROUP_NO=["N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"], DX=0.0, DY=0.0, DZ=0.0),
        ),
    )

    ################################################################################

    #  RESOLUTION NON LIN
    #

    RES = STAT_NON_LINE(
        MODELE=MODELE,
        CHAM_MATER=CHMAT,
        CARA_ELEM=CAREL,
        EXCIT=(_F(CHARGE=CLF, FONC_MULT=FONCT), _F(CHARGE=CLD)),
        COMPORTEMENT=_F(
            RELATION=nomComportement, GROUP_MA=Grma27, RESI_INTE=1e-8, ITER_INTE_MAXI=15
        ),
        INCREMENT=_F(LIST_INST=LINST),
    )

    RES = CALC_CHAMP(reuse=RES, RESULTAT=RES, DEPLACEMENT=("SAUT_ELNO"), VARI_INTERNE=("VARI_ELNO"))

    ################################################################################

    #  TESTS_RESU
    #

    # Matrice tangente

    preci = 1e-5
    # résultats analytiques
    l_res = [0.0, 1.0, 2.0, 0.0, -1.0]
    # valeurs de non-régression
    l_nreg = [0.0, 0.9999999999999996, 1.9999999999999982, 5.4225157697374e-10, -1.0000000016336237]

    for i in range(len(l_inst2)):

        crit = "RELATIF"
        if abs(l_res[i]) < preci:
            crit = "ABSOLU"

        TEST_RESU(
            RESU=_F(
                INST=l_inst2[i],
                GROUP_MA=Grma27,
                POINT=23,
                REFERENCE="ANALYTIQUE",
                RESULTAT=RES,
                NOM_CHAM="SAUT_ELNO",
                NOM_CMP="DY",
                VALE_CALC=l_nreg[i],
                VALE_REFE=l_res[i],
                CRITERE=crit,
                PRECISION=preci,
            )
        )

    # Ecrouissage de la surface de charge

    if optionSurfCharge:

        TAB = POST_RELEVE_T(
            ACTION=_F(
                INTITULE="TABLE",
                OPERATION="EXTRACTION",
                GROUP_NO="N23",
                RESULTAT=RES,
                NOM_CHAM="VARI_ELNO",
                TOUT_CMP="OUI",
            )
        )

        tabext = TAB.EXTR_TABLE()
        tabval = tabext.values()
        tabval2 = {key: tabval[key] for key in ["INST", "V1", "V2", "V3"]}

        assert len(l_inst2) == len(tabval2["INST"])

        for i in range(len(tabval2["INST"])):
            f_n = l_valy2[i] * f_ny
            upl_n = np.array([tabval2["V1"][i], tabval2["V2"][i], tabval2["V3"][i]])
            G = np.array([g_t, g_n, g_n])
            r = G * upl_n
            r_n = r[1]
            f_n_2 = np.sign(f_n - r_n) * f_ny + r_n
            #
            test.assertAlmostEqual(f_n, f_n_2)

    # ---------------------------------------------------------------------------- #
    CA.close()


# ---------------------------------------------------------------------------- #
testLoi(nombrePas=4, nomComportement="INTERF_POU_CINE", optionSurfCharge=True)
