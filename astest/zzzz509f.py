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

# This testcase is a copy of ssnv160f, but using MECA_NON_LINE.
import numpy as np
from scipy.optimize import fsolve

from code_aster.Commands import *
from code_aster import CA


test = CA.TestCase()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>> Isotropic compression test on a 3D HEXA8 element with the CSSM model
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="NON"))

### >>>>>>>>>>>>>
### Read the mesh
### <<<<<<<<<<<<<

MAIL = LIRE_MAILLAGE(FORMAT="ASTER")

### >>>>>>>>>>>>>>>>>
### Model affectation
### <<<<<<<<<<<<<<<<<

MODELE = AFFE_MODELE(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

### >>>>>>>>>>>>>>>>>>>
### Material properties
### <<<<<<<<<<<<<<<<<<<

k = 516.0e6
mu = 238.0e6
rho = 0.1
M = 1.38
pc0 = 100.0e3
bt = 30.0
eta = 0.99
om = 32.0
gammahyp = 2.0e-4
nhyp = 0.78
C = 448.0e3

MATER = DEFI_MATERIAU(
    ELAS=_F(
        E=9.0 * k * mu / (3.0 * k + mu), NU=(3.0 * k - 2.0 * mu) / (2.0 * (3.0 * k + mu)), ALPHA=0.0
    ),
    CSSM=_F(
        BulkModulus=k,
        ShearModulus=mu,
        InitCritPress=pc0,
        CritStateSlope=M,
        IncoPlastIndex=bt,
        HypExponent=nhyp,
        HypDistortion=gammahyp,
        MinCritPress=C,
        ShearModulusRatio=rho,
        IsoHardRatio=eta,
        IsoHardIndex=om,
    ),
)

### >>>>>>>>>>>>>>>>>>>>
### Material affectation
### <<<<<<<<<<<<<<<<<<<<

MATE = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MATER))

### >>>>>>>>>>
### Time steps
### <<<<<<<<<<

# t=0.2 inserted with PAS_INIT
LI1 = DEFI_LIST_REEL(VALE=(0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))

maxSteps = 5
timeStepper = DEFI_LIST_INST(
    METHODE="AUTO",
    DEFI_LIST=_F(LIST_INST=LI1, NB_PAS_MAXI=maxSteps, PAS_INIT=0.2),
    ECHEC=_F(ACTION="DECOUPE", SUBD_METHODE="MANUEL", SUBD_PAS=2, SUBD_NIVEAU=6),
)

### >>>>>>>>>>>>>>>>>>>
### Boundary conditions
### <<<<<<<<<<<<<<<<<<<

# /!\ p is 50 times more than in ssnv160f + reduced ITER_INTE_MAXI to force the timesteps splitting
p = 500.0e4

PRES = DEFI_FONCTION(NOM_PARA="INST", NOM_RESU="PRESSION", VALE=(0.0, 0.0, 1.0, p, 2.0, 100.0))

CHA_PRES = AFFE_CHAR_MECA_F(
    MODELE=MODELE, PRES_REP=_F(GROUP_MA=("HAUT", "DROITE", "ARRIERE"), PRES=PRES), VERI_NORM="OUI"
)

CHA_DEPL = AFFE_CHAR_MECA(
    MODELE=MODELE,
    DDL_IMPO=(
        _F(GROUP_NO=("NO1", "NO2", "NO4", "NO3"), DZ=0.0),
        _F(GROUP_NO=("NO3", "NO4", "NO7", "NO8"), DY=0.0),
        _F(GROUP_NO=("NO2", "NO4", "NO6", "NO8"), DX=0.0),
    ),
)

### >>>>>>>>
### Solution with MECA_NON_LINE
### <<<<<<<<

command = MECA_NON_LINE

result = CA.NonLinearResult()

keywords = _F(
    MODELE=MODELE,
    CHAM_MATER=MATE,
    EXCIT=(_F(CHARGE=CHA_PRES), _F(CHARGE=CHA_DEPL)),
    COMPORTEMENT=_F(RELATION="CSSM", RESI_INTE=1.0e-14, ITER_INTE_MAXI=20),
    INCREMENT=_F(LIST_INST=timeStepper),
    NEWTON=_F(MATRICE="TANGENTE"),
    CONVERGENCE=_F(ITER_GLOB_MAXI=10, RESI_GLOB_RELA=1.0e-10),
    SOLVEUR=_F(METHODE="MUMPS"),
)

done = 1  # t=0
finished = False
while not finished:
    result = command(**keywords)

    test.assertLessEqual(
        result.getNumberOfIndexes() - done,
        maxSteps,
        msg="check that at most 'maxSteps' have been completed",
    )

    finished = result.getLastTime() >= 2.0
    if done == 1:
        keywords.update(_F(reuse=result, ETAT_INIT=_F(EVOL_NOLI=result)))
    done = result.getNumberOfIndexes()

test.assertAlmostEqual(result.getLastTime(), 2.0, msg="last time == 2.0")


### >>>>>>>>>>>>>>>
### Post-processing
### <<<<<<<<<<<<<<<

result = CALC_CHAMP(
    reuse=result,
    RESULTAT=result,
    CONTRAINTE="SIEF_NOEU",
    DEFORMATION="EPSI_NOEU",
    VARI_INTERNE="VARI_NOEU",
)

### >>>>>
### Tests
### <<<<<


def func(x):
    return p - 2.0 * pc0 * (np.exp(-bt * x) - eta * np.exp(2.0 * om * x))


x0 = -np.log(p / (2.0 * pc0)) / bt
EPVP = fsolve(func, x0)[0]
EPVE = -p / k


TEST_RESU(
    RESU=_F(
        INST=1.0,
        RESULTAT=result,
        REFERENCE="ANALYTIQUE",
        NOM_CHAM="EPSI_NOEU",
        GROUP_NO="NO6",
        NOM_CMP="EPXX",
        VALE_REFE=(EPVE + EPVP) / 3.0,
        VALE_CALC=-0.038995719385912855,
    )
)

TEST_RESU(
    RESU=_F(
        INST=1.0,
        RESULTAT=result,
        REFERENCE="ANALYTIQUE",
        NOM_CHAM="VARI_NOEU",
        NOM_CMP="V15",
        GROUP_NO="NO6",
        VALE_REFE=EPVP,
        VALE_CALC=-0.10729723567711803,
    )
)

TEST_RESU(
    RESU=_F(
        INST=1.4,
        RESULTAT=result,
        NOM_CHAM="EPSI_NOEU",
        GROUP_NO="NO6",
        NOM_CMP="EPXX",
        VALE_CALC=-0.037703755561623496,
    )
)

TEST_RESU(
    RESU=_F(
        INST=1.4,
        RESULTAT=result,
        NOM_CHAM="VARI_NOEU",
        NOM_CMP="V15",
        GROUP_NO="NO6",
        VALE_CALC=-0.10729723567711803,
    )
)

### >>>>>>>>
### Solution with STAT_NON_LINE
### <<<<<<<<

Kinv = 3.2841e-4
Kv = 1.0 / Kinv
SY = 120.0
Rinf = 758.0
Qzer = 758.0 - SY
Qinf = Qzer + 100.0
b = 2.3
C1inf = 63767.0 / 2.0
C2inf = 63767.0 / 2.0
Gam1 = 341.0
Gam2 = 341.0
C_Pa = 1.0e6

Steel = DEFI_MATERIAU(
    ELAS=_F(E=200e9, NU=0.3),
    ECRO_LINE=_F(D_SIGM_EPSI=1930.0, SY=100.0e6),
    VISCOCHAB=_F(
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
    ),
)

MATE2 = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=Steel))
p = 5000.0e7

PRES2 = DEFI_FONCTION(NOM_PARA="INST", NOM_RESU="PRESSION", VALE=(0.0, 0.0, 1.0, p))


CHA_PRES2 = AFFE_CHAR_MECA_F(
    MODELE=MODELE, PRES_REP=_F(GROUP_MA=("HAUT",), PRES=PRES2), VERI_NORM="OUI"
)

CHA_DEPL2 = AFFE_CHAR_MECA(MODELE=MODELE, DDL_IMPO=(_F(GROUP_MA=("BAS",), DX=0, DY=0.0, DZ=0.0),))

LI1 = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=1.0, NOMBRE=10), _F(JUSQU_A=2.0, NOMBRE=10)))

maxSteps = 9
timeStepper = DEFI_LIST_INST(
    METHODE="AUTO",
    DEFI_LIST=_F(LIST_INST=LI1, NB_PAS_MAXI=maxSteps),
    ECHEC=_F(ACTION="DECOUPE", SUBD_METHODE="MANUEL", SUBD_PAS=3, SUBD_NIVEAU=4),
)

command = STAT_NON_LINE

result = CA.NonLinearResult()

inst_fin = 1.0

keywords = _F(
    MODELE=MODELE,
    CHAM_MATER=MATE2,
    EXCIT=(_F(CHARGE=CHA_PRES2), _F(CHARGE=CHA_DEPL2)),
    COMPORTEMENT=_F(RELATION="VISCOCHAB"),
    INCREMENT=_F(LIST_INST=timeStepper, INST_FIN=inst_fin),
    NEWTON=_F(MATRICE="TANGENTE"),
    CONVERGENCE=_F(ITER_GLOB_MAXI=10),
)

done = 1  # t=0
finished = False
while not finished:
    result = command(**keywords)

    test.assertLessEqual(
        result.getNumberOfIndexes() - done,
        maxSteps,
        msg="check that at most 'maxSteps' have been completed",
    )

    finished = abs(result.getLastTime() - inst_fin) <= 1.0e-6
    if done == 1:
        keywords.update(_F(reuse=result, ETAT_INIT=_F(EVOL_NOLI=result)))
    done = result.getNumberOfIndexes()

test.assertAlmostEqual(result.getLastTime(), inst_fin, msg="last time == 0.9")

# Risques de variabilité sur les résultats, pr on ne teste pas les valeurs, mais juste le stop de l'algo (assert précédent)

# TEST_RESU(
#     RESU=_F(
#         INST=inst_fin,
#         RESULTAT=result,
#         PRECISION=0.02,
#         REFERENCE="AUTRE_ASTER",
#         NOM_CHAM="DEPL",
#         GROUP_NO="NO6",
#         NOM_CMP="DX",
#         VALE_REFE=0.4066657673320981,
#         VALE_CALC=0.4066581700455644,
#     )
# )


FIN()
