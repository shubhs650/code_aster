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
from code_aster.Messages import MessageLog

CA.init("--test", IGNORE_ALARM="MATERIAL2_20")

test = CA.TestCase()

mesh = CA.Mesh.buildCube(refine=0)

trac = DEFI_FONCTION(
    NOM_PARA="EPSI",
    ABSCISSE=(0.002, 0.1, 0.15),
    ORDONNEE=(300.0, 300.0, 300.0),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

cst1 = DEFI_CONSTANTE(VALE=1.0)
cst2 = DEFI_CONSTANTE(VALE=2.0)

mater = CA.Material()
mater.addProperties("ELAS", E=3.7272000000e10, NU=0.0, RHO=2400.0)
mater.addProperties("TRACTION", SIGM=trac)
mater.addProperties("MFRONT", LISTE_COEF=(3.69e10, 0.3, 151.0, 87.0, 2.3))
mater.addProperties("UMAT_FO", LISTE_COEF=(cst1, cst2))
mater.addProperties("DRUCK_PRAGER", ALPHA=1.0, SY=2.0, P_ULTM=3.0, ECROUISSAGE="LINEAIRE", H=4.0)
mater.addProperties("ROUSSELIER_FO", D=cst1, SIGM_1=cst2, PORO_INIT=cst1)

chmat = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=mater))

# check for Material API
print("\nchecking material assignment...")
list_affe = chmat.getVectorOfPartOfMaterialField()
test.assertEqual(len(list_affe), 1, msg="number of zones")

affe = list_affe[0]
list_mate = list_affe[0].getVectorOfMaterial()
test.assertEqual(len(list_mate), 1, msg="number of materials")

mat = list_mate[0]
print("\nchecking material...")
test.assertEqual(mater, mat, msg="same objects")

# MaterialProperty == factor keyword in DEFI_MATERIAU
test.assertEqual(mat.size(), 6, msg="number of material properties")
test.assertCountEqual(
    mat.getMaterialNames(), ["ELAS", "TRACTION", "MFRONT", "DRUCK_PRAGER", "ROUSSELIER", "UMAT"]
)

print("\nchecking ELAS...")
test.assertIsNone(mat.getFunction("ELAS", "E"), msg="no function")
test.assertIsNone(mat.getFunction("ELAS", "SIGM"), msg="no function")

print("\nchecking TRACTION...")
fSIGM = mat.getFunction("TRACTION", "SIGM")
test.assertEqual(fSIGM, trac, msg="retreive 'trac'")

Kinv = 3.2841e-4
Kv = 1.0 / Kinv
SY = 437.0
Rinf = 758.0
Qzer = 758.0 - 437.0
Qinf = Qzer + 100.0
b = 2.3
C1inf = 63767.0 / 2.0
C2inf = 63767.0 / 2.0
Gam1 = 341.0
Gam2 = 341.0
C_Pa = 1.0e6

matvisco = CA.Material()
matvisco.addProperties("ELAS", E=3.7272000000e10, NU=0.0, RHO=2400.0)
matvisco.addProperties(
    "VISCOCHAB",
    K=SY * C_Pa,
    B=b,
    MU=10,
    Q_M=Qinf * C_Pa,
    Q_0=Qzer * C_Pa,
    C1=C1inf * C_Pa,
    C2=C2inf * C_Pa,
    G1_0=Gam1,
    G2_0=Gam2,
    K_0=Kv * C_Pa,
    N=11,
    A_K=1.0,
)

test.assertEqual(matvisco.size(), 2, msg="number of material properties")
test.assertCountEqual(matvisco.getMaterialNames(), ["ELAS", "VISCOCHAB"])

E = 2.0e05
NU = 0.3
H = 1.0
C_MEM = E * H / (1.0 - (NU * NU))
C_FLE = (E * H * H * H) / (12.0 * (1.0 - (NU * NU)))
C_CIS = E * H / (1.0 + NU)
C1111 = C_MEM * 1.0
C1112 = C_MEM * NU
C2222 = C_MEM * 1.0
C1212 = C_MEM * (1.0 - NU) / 2.0
D1111 = C_FLE * 1.0
D1112 = C_FLE * NU
D2222 = C_FLE * 1.0
D1212 = C_FLE * (1.0 - NU) / 2.0
G11 = C_CIS * 5.0 / 12.0
G22 = C_CIS * 5.0 / 12.0

matord = CA.Material()
matord.addProperties(
    "ELAS_COQUE",
    MEMB_L=C1111,
    MEMB_LT=C1112,
    MEMB_T=C2222,
    MEMB_G_LT=C1212,
    FLEX_L=D1111,
    FLEX_LT=D1112,
    FLEX_T=D2222,
    FLEX_G_LT=D1212,
    CISA_L=G11,
    CISA_T=G22,
    RHO=8.0e-6,
    ALPHA=1.0e-05,
)

test.assertEqual(matord.size(), 1, msg="number of material properties")
test.assertCountEqual(matord.getMaterialNames(), ["ELAS_COQUE"])

matord.addProperties("MONO_VISC1", N=0.0, K=1.0e11, C=0.0)

test.assertEqual(matord.size(), 2, msg="number of material properties")
test.assertCountEqual(matord.getMaterialNames(), ["ELAS_COQUE", "MONO_VISC1"])

BIDON = DEFI_CONSTANTE(VALE=1.0)

UN = DEFI_CONSTANTE(VALE=1.0)

ZERO = DEFI_CONSTANTE(VALE=0.0)

VISCOLIQ = DEFI_CONSTANTE(VALE=1.0e-3)

VISCOGAZ = DEFI_CONSTANTE(VALE=1.8e-5)

DVISCOL = DEFI_CONSTANTE(VALE=0.0)

DVISCOG = DEFI_CONSTANTE(VALE=0.0)

KINT = DEFI_CONSTANTE(VALE=1.0e-18)

THMALP1 = DEFI_CONSTANTE(VALE=0.000100)

matthm = CA.Material()
matthm.addProperties("ELAS", E=22.4e6, NU=0.3, RHO=2500.0, ALPHA=1.0e-5)
matthm.addProperties("CJS", BETA_CJS=-0.03, GAMMA_CJS=0.82, RM=0.289, PA=-100.0e3)

# context=_F(COMP_THM="LIQU_SATU"),
matthm.addProperties(
    "THM_LIQU",
    context=dict(COMP_THM="LIQU_SATU"),
    RHO=1000.0,
    UN_SUR_K=0.0,
    ALPHA=THMALP1,
    CP=4180.0,
    VISC=VISCOLIQ,
    D_VISC_TEMP=DVISCOL,
)
matthm.addProperties(
    "THM_GAZ",
    context=_F(COMP_THM="LIQU_SATU"),
    MASS_MOL=28.96e-3,
    CP=1000.0,
    VISC=VISCOGAZ,
    D_VISC_TEMP=DVISCOG,
)
matthm.addProperties(
    "THM_VAPE_GAZ",
    context=_F(COMP_THM="LIQU_SATU"),
    MASS_MOL=18.0e-3,
    CP=1870.0,
    VISC=BIDON,
    D_VISC_TEMP=BIDON,
)
matthm.addProperties(
    "THM_DIFFU",
    context=_F(COMP_THM="LIQU_SATU"),
    R_GAZ=8.315,
    RHO=2400.0,
    CP=800.0,
    BIOT_COEF=1.0,
    SATU_PRES=UN,
    D_SATU_PRES=ZERO,
    PESA_X=0.0,
    PESA_Y=0.0,
    PESA_Z=0.0,
    PERM_IN=KINT,
    PERM_LIQU=UN,
    D_PERM_LIQU_SATU=ZERO,
    PERM_GAZ=UN,
    D_PERM_SATU_GAZ=ZERO,
    D_PERM_PRES_GAZ=ZERO,
)
matthm.addProperties(
    "THM_INIT",
    context=_F(COMP_THM="LIQU_SATU"),
    TEMP=293.0,
    PRE1=1.0,
    PRE2=0.1e6,
    PORO=0.14,
    PRES_VAPE=2269.8,
)

test.assertEqual(matthm.size(), 7, msg="number of material properties")
test.assertCountEqual(
    matthm.getMaterialNames(),
    ["ELAS", "CJS", "THM_LIQU", "THM_GAZ", "THM_VAPE_GAZ", "THM_DIFFU", "THM_INIT"],
)

matrag1 = CA.Material()
matrag1.addProperties("ELAS", E=32000.0e06, NU=0.25)
matrag1.addProperties(
    "BETON_RAG",
    COMP_BETON="ENDO",
    ENDO_MC=1.95,
    ENDO_MT=2.00,
    ENDO_SIGUC=-35.00e06,
    ENDO_SIGUT=-3.18e06,
    ENDO_DRUPRA=0.15,
)

test.assertEqual(matrag1.size(), 2, msg="number of material properties")
test.assertCountEqual(matrag1.getMaterialNames(), ["ELAS", "BETON_RAG"])

matrag2 = CA.Material()
matrag2.addProperties("ELAS", E=32000.0e06, NU=0.25)
matrag2.addProperties(
    "BETON_RAG",
    COMP_BETON="ENDO_FLUA",
    # unités : Pa
    ENDO_MC=1.95,
    ENDO_MT=2.00,
    ENDO_SIGUC=-35.00e06,
    ENDO_SIGUT=-3.18e06,
    ENDO_DRUPRA=0.15,
    # Unités : Pa Jour
    FLUA_SPH_KR=89000.0e06,
    FLUA_SPH_NR=156000.0e06,
    FLUA_SPH_NI=410000.0e06,
    FLUA_DEV_KR=42000.0e06,
    FLUA_DEV_NR=117000.0e06,
    FLUA_DEV_NI=840000.0e06,
)

test.assertEqual(matrag2.size(), 2, msg="number of material properties")
test.assertCountEqual(matrag2.getMaterialNames(), ["ELAS", "BETON_RAG"])

matdis = CA.Material()
matdis.addProperties(
    "DIS_CONTACT",
    COULOMB=0.5,
    DIST_1=0.0051,
    INST_COMP_INIT=(-1.0, 0.0),
    RIGI_NOR=1000000.0,
    RIGI_TAN=1000000.0,
)

test.assertEqual(matdis.size(), 1, msg="number of material properties")
test.assertCountEqual(matdis.getMaterialNames(), ["DIS_CONTACT"])

G = 3176517.07 + 2344699.5j
K = complex(2.22e9, 0)
nu = (3 * K - 2 * G) / (3 * K + 2 * G) * 0.5
rho_d = 1460

matcmplx = CA.Material()
matcmplx.addProperties("ELAS_VISCO", G=G, NU=nu, RHO=rho_d)

test.assertEqual(matcmplx.size(), 1, msg="number of material properties")
test.assertCountEqual(matcmplx.getMaterialNames(), ["ELAS_VISCO"])

TRCMNDA = DEFI_TRC(
    HIST_EXP=(
        _F(
            VALE=(
                -1.000e-01,
                1.000e01,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                8.300e02,
                0.000e00,
                0.000e00,
                0.000e00,
                7.591e02,
                1.000e-02,
                0.000e00,
                0.000e00,
                7.550e02,
                6.700e-01,
                0.000e00,
                0.000e00,
                6.200e02,
                6.800e-01,
                0.000e00,
                0.000e00,
                6.159e02,
                6.800e-01,
                0.000e00,
                0.000e00,
                5.247e02,
                6.800e-01,
                0.000e00,
                1.000e-02,
                5.150e02,
                6.800e-01,
                0.000e00,
                3.100e-01,
                3.700e02,
                6.800e-01,
                0.000e00,
                3.200e-01,
                3.603e02,
            )
        ),
        _F(
            VALE=(
                -1.000e00,
                1.000e01,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                8.300e02,
                0.000e00,
                0.000e00,
                0.000e00,
                7.586e02,
                1.000e-02,
                0.000e00,
                0.000e00,
                7.500e02,
                2.900e-01,
                0.000e00,
                0.000e00,
                6.300e02,
                3.000e-01,
                0.000e00,
                0.000e00,
                6.214e02,
                3.000e-01,
                0.000e00,
                0.000e00,
                5.853e02,
                3.000e-01,
                0.000e00,
                1.000e-02,
                5.800e02,
                3.000e-01,
                0.000e00,
                6.900e-01,
                4.000e02,
                3.000e-01,
                0.000e00,
                7.000e-01,
                3.947e02,
            )
        ),
        _F(
            VALE=(
                -1.000e01,
                1.000e01,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                8.300e02,
                0.000e00,
                0.000e00,
                0.000e00,
                5.994e02,
                0.000e00,
                0.000e00,
                1.000e-02,
                5.950e02,
                0.000e00,
                0.000e00,
                9.000e-01,
                4.000e02,
                0.000e00,
                0.000e00,
                9.100e-01,
                3.956e02,
            )
        ),
        _F(
            VALE=(
                -6.000e01,
                1.000e01,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                8.300e02,
                0.000e00,
                0.000e00,
                0.000e00,
                5.094e02,
                0.000e00,
                0.000e00,
                1.000e-02,
                5.000e02,
                0.000e00,
                0.000e00,
                1.900e-01,
                4.150e02,
                0.000e00,
                0.000e00,
                2.000e-01,
                4.056e02,
            )
        ),
    ),
    TEMP_MS=_F(SEUIL=4.500e-01, AKM=-3.125e01, BKM=1.406e01, TPLM=-3.497e00),
    GRAIN_AUST=_F(DREF=11.00e-6, A=11200.0),
)

matmeta = CA.Material()
matmeta.addProperties("THER", RHO_CP=5260000.0, LAMBDA=33.5)
matmeta.addProperties(
    "META_ACIER",
    TRC=TRCMNDA,
    AR3=830.0,
    ALPHA=-0.0306,
    MS0=400.0,
    AC1=724.0,
    AC3=846.0,
    TAUX_1=0.034,
    TAUX_3=0.034,
    LAMBDA0=0.117,
    QSR_K=37500.0,
    D10=3.31,
    WSR_K=12860.0,
)

test.assertEqual(matmeta.size(), 2, msg="number of material properties")
test.assertCountEqual(matmeta.getMaterialNames(), ["THER", "META_ACIER"])

Young = 190000.0e6
fprg = 1800.0e6
kecoul = 0.800646195576
necoul = 8.50471392583
necrou = 1.45855523878
becrou = 49503.9155816
cecrou = 33211.7441074
TempeRef = 20.0

Fkecoul = DEFI_FONCTION(
    NOM_PARA="TEMP",
    VALE=(TempeRef, kecoul, 300.0, kecoul),
    PROL_DROITE="CONSTANT",
    PROL_GAUCHE="CONSTANT",
)

matbpel = CA.Material()
matbpel.addProperties("ELAS", E=Young, NU=0.10, ALPHA=20.0e-05)
matbpel.addProperties(
    "RELAX_ACIER",
    F_PRG=DEFI_CONSTANTE(VALE=fprg),
    ECOU_K=Fkecoul,
    ECOU_N=DEFI_CONSTANTE(VALE=necoul),
    ECRO_N=DEFI_CONSTANTE(VALE=necrou),
    ECRO_B=DEFI_CONSTANTE(VALE=becrou),
    ECRO_C=DEFI_CONSTANTE(VALE=cecrou),
)
matbpel.addProperties("BPEL_ACIER", F_PRG=fprg * 1.10)

test.assertEqual(matbpel.size(), 3, msg="number of material properties")
test.assertCountEqual(matbpel.getMaterialNames(), ["ELAS", "RELAX_ACIER", "BPEL_ACIER"])

TREF = 20.0
T0 = TREF
Tmax = 200.0

ZERO = DEFI_CONSTANTE(VALE=0.0)

YOUN = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 150000.0, Tmax, 100000.0))

ALPH = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 1.0e-5, Tmax, 2.0e-5))

# Parametre du modele elasto-viscoplastique de Chaboche
K_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 25, Tmax, 40))

aphak_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 1.0, Tmax, 1.0))

alphar_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 0.50, Tmax, 0.80))

K0_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 60, Tmax, 80))

N_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 30, Tmax, 15))

alpha_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 0.0, Tmax, 0.0))

b_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 15, Tmax, 15))

mr_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 2.0, Tmax, 2.0))

gamar_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 2.5e-7, Tmax, 1.5e-7))

mu_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 22, Tmax, 16))

Q0_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 40, Tmax, 45))

Qm_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 500, Tmax, 400))

Qr0_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 150, Tmax, 250))

eta_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 0.06, Tmax, 0.03))

C1_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 1600, Tmax, 1600))

m1_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 3, Tmax, 5))

d1_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 0.360e-3, Tmax, 0.420e-3))

gx1_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 2.5e-13, Tmax, 1.5e-13))

g10_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 40, Tmax, 60))
C2_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 55000, Tmax, 55000))

m2_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 5, Tmax, 3.5))

d2_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 0.50e-1, Tmax, 0.60e-1))

gx2_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 0.8e-12, Tmax, 1.5e-12))

g20_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 1500, Tmax, 1000))

ainfi_F = DEFI_FONCTION(NOM_PARA="TEMP", VALE=(T0, 0.41, Tmax, 0.56))

matcomp0 = CA.Material()
matcomp0.addProperties(
    "ELAS_FO",
    ALPHA=ALPH,
    B_ENDOGE=0.0,
    COEF_AMOR=1.0,
    E=YOUN,
    K_DESSIC=0.0,
    NU=ZERO,
    PRECISION=1.0,
    TEMP_DEF_ALPHA=20.0,
)
matcomp0.addProperties(
    "VISCOCHAB_FO",
    ALP=alpha_F,
    A_I=ainfi_F,
    A_K=aphak_F,
    A_R=alphar_F,
    B=b_F,
    C1=C1_F,
    C2=C2_F,
    D1=d1_F,
    D2=d2_F,
    ETA=eta_F,
    G1_0=g10_F,
    G2_0=g20_F,
    G_R=gamar_F,
    G_X1=gx1_F,
    G_X2=gx2_F,
    K=K_F,
    K_0=K0_F,
    MU=mu_F,
    M_1=m1_F,
    M_2=m2_F,
    M_R=mr_F,
    N=N_F,
    QR_0=Qr0_F,
    Q_0=Q0_F,
    Q_M=Qm_F,
)

test.assertEqual(matcomp0.size(), 2, msg="number of material properties")
test.assertCountEqual(matcomp0.getMaterialNames(), ["ELAS", "VISCOCHAB"])

matcomp = CA.Material()
matcomp.addProperties(
    "ELAS", ALPHA=0.0, B_ENDOGE=0.0, COEF_AMOR=1.0, E=149500.0, K_DESSIC=0.0, NU=0.0
)
matcomp.addProperties(
    "VISCOCHAB",
    ALP=0.0,
    A_I=0.4115,
    A_K=1.0,
    A_R=0.503,
    B=15.0,
    C1=1600.0,
    C2=55000.0,
    D1=0.00036060000000000004,
    D2=0.050100000000000006,
    ETA=0.059699999999999996,
    G1_0=40.2,
    G2_0=1495.0,
    G_R=2.4899999999999997e-07,
    G_X1=2.49e-13,
    G_X2=8.07e-13,
    K=25.15,
    K_0=60.2,
    MU=21.94,
    M_1=3.02,
    M_2=4.985,
    M_R=2.0,
    N=29.85,
    QR_0=151.0,
    Q_0=40.05,
    Q_M=499.0,
)

test.assertEqual(matcomp.size(), 2, msg="number of material properties")
test.assertCountEqual(matcomp.getMaterialNames(), ["ELAS", "VISCOCHAB"])

matmeta2 = CA.Material()
matmeta2.addProperties(
    "ELAS_META_FO",
    E=YOUN,
    NU=ZERO,
    F_ALPHA=ALPH,
    C_ALPHA=ALPH,
    TEMP_DEF_ALPHA=20.0,
    PHASE_REFE="FROID",
    EPSF_EPSC_TREF=1.0e-2,
    F1_SY=cst2,
    F2_SY=cst2,
    F3_SY=cst2,
    F4_SY=cst2,
    C_SY=cst2,
    SY_MELANGE=cst1,
)
matmeta2.addProperties(
    "META_TRACTION", SIGM_F1=trac, SIGM_F2=trac, SIGM_F3=trac, SIGM_F4=trac, SIGM_C=trac
)

test.assertEqual(matmeta2.size(), 2, msg="number of material properties")
test.assertCountEqual(matmeta2.getMaterialNames(), ["ELAS_META", "META_TRACTION"])
fE = matmeta2.getFunction("ELAS_META", "E")
test.assertEqual(fE.userName, YOUN.userName)
fSY = matmeta2.getFunction("ELAS_META", "SY_MELANGE")
test.assertEqual(fSY.userName, cst1.userName)

mater = CA.Material()
mater.addProperties("ELAS", E=3.7272000000e10, NU=0.0, RHO=2400.0)

with test.assertRaisesRegex(RuntimeError, "already defined"):
    mater.addProperties("ELAS", E=1.0, NU=0.0, RHO=2.0)

test.assertEqual(mater.getValueReal("ELAS", "E"), 3.7272e10, msg="get value of E")
with test.assertRaisesRegex(RuntimeError, "property not found"):
    mater.getValueComplex("ELAS", "E")
with test.assertRaisesRegex(RuntimeError, "property not found"):
    mater.getValueReal("ELAS", "unknown")
with test.assertRaisesRegex(RuntimeError, "material not found"):
    mater.getValueReal("material_name", "invalid")

G = 3176517.07 + 2344699.5j
K = complex(2.22e9, 0)
nu = (3 * K - 2 * G) / (3 * K + 2 * G) * 0.5
rho_d = 1460

test.assertEqual(matcmplx.getValueReal("ELAS_VISCO", "RHO"), 1460.0, msg="get value of RHO")
test.assertEqual(matcmplx.getValueComplex("ELAS_VISCO", "G"), G, msg="get value of G")
test.assertEqual(matcmplx.getValueComplex("ELAS_VISCO", "NU"), nu, msg="get value of NU")

with test.assertRaisesRegex(RuntimeError, "property not found"):
    matcmplx.getValueReal("ELAS_VISCO", "G")

EL = 2.1e5
ET = 4.0e5
EN = 2.1e5
GLT = 0.45e5
GTN = 0.45e5
GLN = 0.35e5
NULT = 0.075
NULN = 0.075
NUTN = 0.0142857143

matorth = CA.Material()
matorth.addProperties(
    "ELAS_ORTH",
    E_L=EL,
    E_T=ET,
    E_N=EN,
    G_LT=GLT,
    G_TN=GTN,
    G_LN=GLN,
    NU_LT=NULT,
    NU_LN=NULN,
    NU_TN=NUTN,
)

test.assertEqual(matorth.size(), 1, msg="number of material properties")
test.assertCountEqual(matorth.getMaterialNames(), ["ELAS_ORTH"])

infos = MessageLog.get_info_alarm()
emitted = [alr[1] for alr in infos if alr[0] == "MATERIAL2_20"]
test.assertEqual(len(emitted), 0, msg="check for warning MATERIAL2_20: off")

# this will raised a warning on checking eigen value of the Hooke matrix
NULN = 0.075 + 1

matorth = CA.Material()
matorth.addProperties(
    "ELAS_ORTH",
    E_L=EL,
    E_T=ET,
    E_N=EN,
    G_LT=GLT,
    G_TN=GTN,
    G_LN=GLN,
    NU_LT=NULT,
    NU_LN=NULN,
    NU_TN=NUTN,
)

test.assertEqual(matorth.size(), 1, msg="number of material properties")
test.assertCountEqual(matorth.getMaterialNames(), ["ELAS_ORTH"])

infos = MessageLog.get_info_alarm()
emitted = [alr[1] for alr in infos if alr[0] == "MATERIAL2_20"]
test.assertEqual(emitted[0] if emitted else 0, 1, msg="check for warning MATERIAL2_20: on")

test.printSummary()

FIN()
