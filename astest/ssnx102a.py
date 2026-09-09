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

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"))


def getConnectivityOfGroup(grma, mesh):
    """
    Get connectivity of GROUP_MA
    """
    conn = mesh.getConnectivity()
    gr_cells = mesh.getCells(grma)
    con = [item for cell in gr_cells for item in conn[cell]]
    con.insert(0, None)
    return con


def nautelemorie(gx):
    """
    Get nautical angles (in degrees) from a 3D-vector.

    The choice of angles is fully constrained and non-ambiguous.

    @angvx.F90
    """
    import math

    tst = 1e-9

    if abs(gx[1]) <= tst and abs(gx[0]) <= tst:
        alpha = 0.0
    else:
        alpha = math.atan2(gx[1], gx[0])

    p = math.sqrt(gx[0] ** 2 + gx[1] ** 2)
    if abs(gx[2]) <= tst and abs(p) <= tst:
        beta = 0.0
    else:
        beta = -math.atan2(gx[2], p)

    alpha = np.rad2deg(alpha)
    beta = np.rad2deg(beta)

    return alpha, beta


def rotmat(angl_naut):
    """
    Construct a rotation matrix from nautical angles.

    Parameters:
    angl_naut (list or np.array): Nautical angles (Alpha, Beta, Gamma)

    Returns:
    np.array: Rotation matrix (3x3)
    """

    # Extract angles
    alpha = np.deg2rad(angl_naut[0])
    beta = np.deg2rad(angl_naut[1])
    gamma = np.deg2rad(angl_naut[2])

    # Calculate sine and cosine of angles
    cosa = np.cos(alpha)
    sina = np.sin(alpha)
    cosb = np.cos(beta)
    sinb = np.sin(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)

    # Initialize the rotation matrix
    pgl = np.zeros((3, 3))

    # Fill the rotation matrix
    pgl[0, 0] = cosb * cosa
    pgl[1, 0] = sing * sinb * cosa - cosg * sina
    pgl[2, 0] = sing * sina + cosg * sinb * cosa
    pgl[0, 1] = cosb * sina
    pgl[1, 1] = cosg * cosa + sing * sinb * sina
    pgl[2, 1] = cosg * sinb * sina - cosa * sing
    pgl[0, 2] = -sinb
    pgl[1, 2] = sing * cosb
    pgl[2, 2] = cosg * cosb

    return pgl


def carapou(sectyp="RECTANGLE", grma="Hexa27"):
    """
    Definition du type et caractéristiques de section du pieu
    """

    assert sectyp in ["RECTANGLE", "CERCLE"], "The section type should be 'RECTANGLE' or 'CERCLE'."

    if sectyp == "RECTANGLE":
        hy = 0.5
        hz = 1.0
        facts_loc = (2 * (hy + hz), hy, hz)

        carsec = (_F(GROUP_MA=grma, SECTION="RECTANGLE", CARA=("HY", "HZ"), VALE=(hy, hz)),)

    elif sectyp == "CERCLE":
        r = 0.2
        d = 2 * r
        facts_loc = (np.pi * d, d, d)

        carsec = (_F(GROUP_MA=grma, SECTION="CERCLE", CARA=("R"), VALE=(r,)),)

    return carsec, facts_loc


def calcSRef(uimp, vOrie, LX, anglvril, kddll):
    """
    Calculate reference solutions.
    """

    # orientation élément
    alpha, beta = nautelemorie(vOrie)
    matrotEl = rotmat([alpha, beta, 0.0])
    tmatrotEl = np.transpose(matrotEl)
    # orientation section
    matrotSec = rotmat([0.0, 0.0, anglvril])
    tmatrotSec = np.transpose(matrotSec)

    # Ajustement solution pour les rotations imposées
    switchrot2tran = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])

    # déplacements imposés
    # uimpPR = np.matmul(switchrot2tran, uimp)
    uimpPT = uimp
    uimpST = uimp

    # facteurs de réponse analytiques
    afactanaPR = np.ones(3) * LX**2 / 12
    afactanaPR[0] = 0.0
    afactanaPT = np.ones(3) * 1.0 / 2 * LX
    afactanaST = np.ones(3) * 1.0 / 2 * LX / 4

    # calcul des efforts solution
    mId3 = np.identity(3)

    func_Uloc = lambda uimp, switchrot2tran: np.matmul(
        matrotSec, np.matmul(switchrot2tran, np.matmul(matrotEl, uimp))
    )
    func_Sloc = lambda Uloc: np.multiply(kddll, Uloc)
    func_Ssec = lambda Sloc: np.matmul(tmatrotSec, Sloc)
    func_Sglob = lambda Ssec: np.matmul(tmatrotEl, Ssec)
    func_Fglob = lambda Ssec, afactana: np.matmul(tmatrotEl, np.multiply(afactana, Ssec))
    UlocST = func_Uloc(uimpST, mId3)
    SlocST = func_Sloc(UlocST)
    SsecST = func_Ssec(SlocST)
    SglobST = func_Sglob(SsecST)
    FglobST = func_Fglob(SsecST, afactanaST)

    FglobPR = func_Fglob(func_Ssec(func_Sloc(func_Uloc(uimpST, switchrot2tran))), afactanaPR)
    FglobPT = func_Fglob(func_Ssec(func_Sloc(func_Uloc(uimpPT, mId3))), afactanaPT)

    return FglobPR, FglobPT, FglobST, SglobST


def getValReg(testnum, restype="STAT_NON_LINE"):
    """
    Get numerical reference values for non-regression test.
    """

    if testnum == 0:
        lvalReg = [
            0.371230501935536,
            18.43161917689161,
            -2.2731312297060983e-17,
            23.576469276950043,
            10.818197573087605,
            22.05,
            5.894117319237511,
            2.704549393271901,
            5.512500000000001,
            # 9.0,
            # 2.5,
            # 5.0,
            1.423512627623142,
            0.6531868988741312,
            1.331346652052679,
            5.312621451505441,
            2.4377266934366126,
            4.968653347947321,
            0.0,
            -4.440892098500626e-16,
            0.0,
            6.736134079128583,
            3.090913592310744,
            6.3,
            -1.1102230246251565e-16,
            -5.551115123125783e-17,
            0.0,
            1.4999999999999998,
            0.39999999999999997,
            0.6999999999999997,
            -4.04145188432738,
            4.04145188432738,
        ]

    elif testnum == 1:
        lvalReg = [
            0.20999999999999974,
            -1.2858791391047208e-17,
            -4.8000000000000025,
            9.599999999999998,
            9.613296,
            0.42000000000000026,
            2.4000000000000004,
            2.4033239999999996,
            0.10499999999999995,
            # 3.76992,
            # 2.0000000000000004,
            # 2.0000000000000004,
            1.3524791385931978,
            1.3543523222001492,
            0.05917096231345237,
            5.047520861406803,
            5.054511677799851,
            0.22082903768654738,
            -4.440892098500626e-16,
            0.0,
            4.163336342344337e-17,
            6.4,
            6.4088639999999995,
            0.2799999999999997,
            -1.1102230246251565e-16,
            -5.551115123125783e-17,
            2.0816681711721685e-17,
            3.1999999999999997,
            1.6999999999999997,
            0.13999999999999985,
            -4.04145188432738,
            4.04145188432738,
        ]

    return lvalReg


def checkOptions(idOrie, sectype, restype, postOptions, nomComportement):
    """
    Check Options
    """
    assert idOrie in [0, 1, 2], "The reinforcement orientation axis should be 0(x), 1(y), or 2(z)."
    assert sectype in ["RECTANGLE", "CERCLE"], "The section type should be 'RECTANGLE' or 'CERCLE'."
    assert restype in [
        "STAT_NON_LINE",
        "MECA_STATIQUE",
    ], "The resolution operattion should be 'STAT_NON_LINE' or 'MECA_STATIQUE'."
    for option in postOptions:
        assert option in [
            "FORC_NODA",
            "SIEF_ELGA",
            "SIEF_ELNO",
            "SAUT_ELNO",
            "COOR_ELGA",
        ], "The postprocessing fields should be in ['FORC_NODA', 'SIEF_ELGA', 'SIEF_ELNO', 'SAUT_ELNO', 'COOR_ELGA']."
    assert nomComportement in [
        "INTERF_POU_ELAS",
        "INTERF_POU_CINE",
    ], "The behaviour for STAT_NON_LINE should be 'INTERF_POU_ELAS' or 'INTERF_POU_CINE'."


def faireTest(
    idOrie, anglvril, sectype, nomComportement, restype, dimp, testOptions, valRegression
):
    """
    Testcase for 3D_INTERF_POU element.

    testing :
    - element orientation on X (idOrie=0), Y (idOrie=1) or Z (idOrie=2) : MODI_MAILLAGE(ORIE_INTERF_POU)
    - section orientation : anglvril in degrees. ORIENTATION(CARA="ANGL_VRIL")
    - section type : "RECTANGLE" or "CERCLE" (with fixed geometrical properties)
    - calculation options : STAT_NON_LINE, RIGI_MECA, FORC_NODA, SIEF_ELGA, SIEF_ELNO, SAUT_ELGA, SAUT_ELNO
    """
    ################################################################################

    #  VÉRIFICATION DES ENTRÉES
    #

    checkOptions(idOrie, sectype, restype, testOptions, nomComportement)

    ################################################################################

    #  MAILLAGE
    #

    Grma27 = "Hexa27"

    # vecteur d'orientation de l'axe de la poutre
    vOrie = [0.0, 0.0, 0.0]
    vOrie[idOrie] = 1.0

    MA = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=20)

    MA = MODI_MAILLAGE(reuse=MA, MAILLAGE=MA, ORIE_INTERF_POU=_F(GROUP_MA=Grma27, VECT_ORIE=vOrie))

    a_conn = np.asarray(getConnectivityOfGroup(Grma27, MA))

    # longueur du pieu
    l_coor = MA.getCoordinates()
    l_coorPy = l_coor.toNumpy()
    LX = abs(l_coorPy[a_conn[25]][idOrie] - l_coorPy[a_conn[23]][idOrie])

    # noeuds pour les conditions limites
    MA.setGroupOfNodes("GSXM", a_conn[[1, 4, 5, 8]])
    MA.setGroupOfNodes("GSXP", a_conn[[2, 3, 6, 7]])
    MA.setGroupOfNodes("GS", a_conn[range(1, 9)])
    MA.setGroupOfNodes("GPXM", a_conn[[25]])
    MA.setGroupOfNodes("GPXP", a_conn[[23]])
    MA.setGroupOfNodes("GP", a_conn[[25, 23]])

    # noeuds pour les posttraitements
    MA.setGroupOfNodes("GSXM1", a_conn[[1]])
    MA.setGroupOfNodes("GSXP1", a_conn[[2]])
    MA = DEFI_GROUP(
        MAILLAGE=MA,
        reuse=MA,
        CREA_GROUP_NO=(
            _F(NOM="GRNP", UNION=("GPXM", "GPXP")),
            _F(NOM="GRNS", UNION=["GSXM1", "GSXP1"]),
        ),
    )

    ################################################################################

    #  CARACTERISTIQUES GEOMETRIQUES
    #

    # pieu
    carsec, facts_loc = carapou(sectyp=sectype, grma=Grma27)

    ################################################################################

    #  CARACTERISTIQUES MATERIAU
    #

    # pieu
    E_p = 1
    nu_p = 0.25

    # interface
    k_n = 5
    k_t = 3
    f_ny = k_n * 100
    f_ty = k_t * 100
    g_n = 0.0
    g_t = 0.0

    # ATTENTION :
    #  Les valeurs des paramètres matériau n'ont pas de sens physique.
    #  Elles servent uniquement à des fins de vérification.

    ################################################################################

    #  CHARGEMENT
    #

    ldimp = dimp

    ################################################################################

    #  DEFINITION DES MATERIAUX
    #

    POU = DEFI_MATERIAU(ELAS=_F(E=E_p, NU=nu_p))

    INT = DEFI_MATERIAU(
        INTERF_POU_ELAS=_F(K_N=k_n, K_T=k_t),
        INTERF_POU_CINE=_F(F_NY=f_ny, F_TY=f_ty, G_N=g_n, G_T=g_t),
    )

    ################################################################################

    #  AFFECTATIONS AU MODELE
    #

    MODELE = AFFE_MODELE(
        MAILLAGE=MA,
        AFFE=(_F(GROUP_MA=Grma27, PHENOMENE="MECANIQUE", MODELISATION="3D_INTERF_POU"),),
    )

    CAREL = AFFE_CARA_ELEM(
        MODELE=MODELE,
        POUTRE=carsec,
        ORIENTATION=(_F(GROUP_MA=Grma27, CARA="ANGL_VRIL", VALE=anglvril),),
    )

    CHMAT = AFFE_MATERIAU(MODELE=MODELE, AFFE=(_F(GROUP_MA=Grma27, MATER=(INT, POU)),))

    ############################################

    #  CHARGEMENT
    #

    FONCT = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1.0, 1.0))

    LINST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=1.0, NOMBRE=1),))

    ################################################################################

    #  CONDITIONS LIMITES
    #

    # déplacement imposé sur la poutre
    CLPT = AFFE_CHAR_CINE(
        MODELE=MODELE,
        MECA_IMPO=(
            _F(GROUP_NO="GS", DX=0.0, DY=0.0, DZ=0.0),
            _F(GROUP_NO="GPXM", DX=0.0, DY=0.0, DZ=0.0, DRX=0.0, DRY=0.0, DRZ=0.0),
            _F(GROUP_NO="GPXP", DX=ldimp[0], DY=ldimp[1], DZ=ldimp[2], DRX=0.0, DRY=0.0, DRZ=0.0),
        ),
    )

    # rotation imposée sur la poutre
    CLPR = AFFE_CHAR_CINE(
        MODELE=MODELE,
        MECA_IMPO=(
            _F(GROUP_NO="GS", DX=0.0, DY=0.0, DZ=0.0),
            _F(GROUP_NO="GPXM", DX=0.0, DY=0.0, DZ=0.0, DRX=0.0, DRY=0.0, DRZ=0.0),
            _F(GROUP_NO="GPXP", DX=0.0, DY=0.0, DZ=0.0, DRX=ldimp[0], DRY=ldimp[1], DRZ=ldimp[2]),
        ),
    )

    # déplacement imposé sur le sol
    CLST = AFFE_CHAR_CINE(
        MODELE=MODELE,
        MECA_IMPO=(
            _F(GROUP_NO="GP", DX=0.0, DY=0.0, DZ=0.0, DRX=0.0, DRY=0.0, DRZ=0.0),
            _F(GROUP_NO="GSXM", DX=0.0, DY=0.0, DZ=0.0),
            _F(GROUP_NO="GSXP", DX=ldimp[0], DY=ldimp[1], DZ=ldimp[2]),
        ),
    )

    l_CL = [CLPR, CLPT, CLST]

    ################################################################################

    #  RESOLUTION LIN - NON LIN
    #

    l_RES = []

    for iCL in l_CL:

        if restype == "MECA_STATIQUE":

            RES = MECA_STATIQUE(
                MODELE=MODELE,
                CHAM_MATER=CHMAT,
                CARA_ELEM=CAREL,
                EXCIT=(_F(CHARGE=iCL, FONC_MULT=FONCT),),
                LIST_INST=LINST,
            )

            RES = CALC_CHAMP(reuse=RES, RESULTAT=RES, CONTRAINTE=("SIEF_ELGA",))

        elif restype == "STAT_NON_LINE":

            RES = STAT_NON_LINE(
                MODELE=MODELE,
                CHAM_MATER=CHMAT,
                CARA_ELEM=CAREL,
                EXCIT=(_F(CHARGE=iCL, FONC_MULT=FONCT),),
                COMPORTEMENT=_F(
                    RELATION=nomComportement, GROUP_MA=Grma27, RESI_INTE=1e-8, ITER_INTE_PAS=-4
                ),
                INCREMENT=_F(LIST_INST=LINST),
            )

        RES = CALC_CHAMP(
            reuse=RES,
            RESULTAT=RES,
            FORCE=("FORC_NODA",),
            DEPLACEMENT=("SAUT_ELGA", "SAUT_ELNO"),
            CONTRAINTE=("SIEF_ELNO", "SIEF_NOEU"),
        )

        l_RES.append(RES)

    ################################################################################

    #  TESTS_RESU
    #

    # Solutions de référence de non-régression
    l_valCalc = valRegression
    iValCalc = 0

    # Solutions de référence analytiques
    kddll = [k_t * facts_loc[0], k_n * facts_loc[1], k_n * facts_loc[2]]
    FresPR, FresPT, FresST, SresST = calcSRef(ldimp, vOrie, LX, anglvril, kddll)

    if "FORC_NODA" in testOptions:

        ####################

        # ROTATION IMPOSEE sur 1 noeud du PIEU

        iRES = l_RES[0]

        TAB = POST_RELEVE_T(
            ACTION=(
                _F(
                    INTITULE="TABLE",
                    OPERATION="EXTRACTION",
                    GROUP_NO="GRNP",
                    RESULTAT=iRES,
                    NOM_CHAM="FORC_NODA",
                    TOUT_CMP="OUI",
                ),
            )
        )

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresPR[0],
            NOM_PARA="DX",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(
                _F(NOM_PARA="INST", VALE=1.0),
                # _F(NOM_PARA="NOM_CHAM", VALE_K="FORC_NODA"),
            ),
        )

        iValCalc += 1

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresPR[1],
            NOM_PARA="DY",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(_F(NOM_PARA="INST", VALE=1.0),),
        )

        iValCalc += 1

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresPR[2],
            NOM_PARA="DZ",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(_F(NOM_PARA="INST", VALE=1.0),),
        )

        iValCalc += 1

        ####################

        # DEPLACEMENT IMPOSE sur 1 noeud du PIEU

        iRES = l_RES[1]

        TAB = POST_RELEVE_T(
            ACTION=(
                _F(
                    INTITULE="TABLE",
                    OPERATION="EXTRACTION",
                    GROUP_NO="GRNP",
                    RESULTAT=iRES,
                    NOM_CHAM="FORC_NODA",
                    TOUT_CMP="OUI",
                ),
            )
        )

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresPT[0],
            NOM_PARA="DX",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(_F(NOM_PARA="INST", VALE=1.0),),
        )

        iValCalc += 1

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresPT[1],
            NOM_PARA="DY",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(_F(NOM_PARA="INST", VALE=1.0),),
        )

        iValCalc += 1

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresPT[2],
            NOM_PARA="DZ",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(_F(NOM_PARA="INST", VALE=1.0),),
        )

        iValCalc += 1

        ####################

        # DEPLACEMENT IMPOSE sur 4 noeuds du SOL

        iRES = l_RES[2]

        TAB = POST_RELEVE_T(
            ACTION=(
                _F(
                    INTITULE="TABLE",
                    OPERATION="EXTRACTION",
                    GROUP_NO="GRNS",
                    RESULTAT=iRES,
                    NOM_CHAM="FORC_NODA",
                    TOUT_CMP="OUI",
                ),
            )
        )

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresST[0],
            NOM_PARA="DX",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(_F(NOM_PARA="INST", VALE=1.0),),
        )

        iValCalc += 1

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresST[1],
            NOM_PARA="DY",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(_F(NOM_PARA="INST", VALE=1.0),),
        )

        iValCalc += 1

        TEST_TABLE(
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            CRITERE="RELATIF",
            VALE_CALC=l_valCalc[iValCalc],
            VALE_REFE=FresST[2],
            NOM_PARA="DZ",
            TYPE_TEST="SOMM",
            TABLE=TAB,
            FILTRE=(_F(NOM_PARA="INST", VALE=1.0),),
        )

        iValCalc += 1

    FF1 = lambda xi: (1 - xi) / 2
    FF2 = lambda xi: (1 + xi) / 2

    if "SIEF_ELGA" in testOptions:

        # SIEF_ELGA : repère local

        iRES = l_RES[2]
        iValCalc = 9

        for i, xi in enumerate([-np.sqrt(1 / 3), np.sqrt(1 / 3)]):

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=i + 1,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SIEF_ELGA",
                    NOM_CMP="FX",
                    VALE_CALC=l_valCalc[iValCalc],
                    # VALE_REFE=(FF1(xi) * 0. + FF2(xi) * uimprST[0]) * kddll[0],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * SresST[0]),
                    CRITERE="RELATIF",
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=i + 1,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SIEF_ELGA",
                    NOM_CMP="FY",
                    VALE_CALC=l_valCalc[iValCalc],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * SresST[1]),
                    CRITERE="RELATIF",
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=i + 1,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SIEF_ELGA",
                    NOM_CMP="FZ",
                    VALE_CALC=l_valCalc[iValCalc],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * SresST[2]),
                    CRITERE="RELATIF",
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

    if "SIEF_ELNO" in testOptions:

        # SIEF_ELNO : repère local

        iRES = l_RES[2]
        iValCalc = 15

        for i, xi in enumerate([-1.0, 1.0]):

            if i == 0:
                critere = "ABSOLU"
            elif i == 1:
                critere = "RELATIF"

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=-2 * i + 25,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SIEF_ELNO",
                    NOM_CMP="FX",
                    VALE_CALC=l_valCalc[iValCalc],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * SresST[0]),
                    CRITERE=critere,
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=-2 * i + 25,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SIEF_ELNO",
                    NOM_CMP="FY",
                    VALE_CALC=l_valCalc[iValCalc],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * SresST[1]),
                    CRITERE=critere,
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=-2 * i + 25,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SIEF_ELNO",
                    NOM_CMP="FZ",
                    VALE_CALC=l_valCalc[iValCalc],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * SresST[2]),
                    CRITERE=critere,
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

    if "SAUT_ELNO" in testOptions:

        # SAUT_ELNO

        iRES = l_RES[2]
        iValCalc = 21

        for i, xi in enumerate([-1.0, 1.0]):

            if i == 0:
                critere = "ABSOLU"
            elif i == 1:
                critere = "RELATIF"

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=-2 * i + 25,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SAUT_ELNO",
                    NOM_CMP="DX",
                    VALE_CALC=l_valCalc[iValCalc],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * ldimp[0]),
                    CRITERE=critere,
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=-2 * i + 25,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SAUT_ELNO",
                    NOM_CMP="DY",
                    VALE_CALC=l_valCalc[iValCalc],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * ldimp[1]),
                    CRITERE=critere,
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

            TEST_RESU(
                RESU=_F(
                    INST=1.0,
                    GROUP_MA=Grma27,
                    POINT=-2 * i + 25,
                    REFERENCE="ANALYTIQUE",
                    RESULTAT=iRES,
                    NOM_CHAM="SAUT_ELNO",
                    NOM_CMP="DZ",
                    VALE_CALC=l_valCalc[iValCalc],
                    VALE_REFE=(FF1(xi) * 0.0 + FF2(xi) * ldimp[2]),
                    CRITERE=critere,
                    PRECISION=1e-5,
                )
            )

            iValCalc += 1

    if "COOR_ELGA" in testOptions:

        # COOR_ELGA

        iRES = l_RES[2]
        iValCalc = 27

        # coordonnées points de Gauss
        resu_COGA = CALC_CHAM_ELEM(MODELE=MODELE, GROUP_MA=Grma27, OPTION="COOR_ELGA")

        TAB_COGA = CREA_TABLE(
            RESU=_F(INTITULE="TABLE", GROUP_MA=Grma27, CHAM_GD=resu_COGA, TOUT_CMP="OUI")
        )

        stringCoor = "XYZ"

        for i, xi in enumerate([-np.sqrt(1 / 3), np.sqrt(1 / 3)]):

            TEST_TABLE(
                REFERENCE="ANALYTIQUE",
                PRECISION=1e-5,
                CRITERE="RELATIF",
                VALE_CALC=l_valCalc[iValCalc],
                VALE_REFE=xi * LX,
                NOM_PARA="COOR_" + stringCoor[idOrie],
                TABLE=TAB_COGA,
                FILTRE=(_F(NOM_PARA="POINT", VALE_I=i + 1),),
            )

            iValCalc += 1


# ------------------------------------------------------------------------------ #

# Tests

faireTest(
    idOrie=2,
    anglvril=35.0,
    sectype="RECTANGLE",
    nomComportement="INTERF_POU_CINE",
    restype="STAT_NON_LINE",
    dimp=np.array([1.5, 0.4, 0.7]),
    testOptions=["FORC_NODA", "SIEF_ELGA", "SIEF_ELNO", "SAUT_ELNO", "COOR_ELGA"],
    valRegression=getValReg(0, restype="STAT_NON_LINE"),
)

FIN()
