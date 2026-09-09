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


from cataelem.Tools.base_objects import Attribute, objects_from_context

# --------------------------------------------------------------------------------------
# Remarques importantes :
# ------------------------
# Les 10 attributs suivants (qui ont "auto=True") sont "automatiques".
# Ils se déduisent d'informations données dans le catalogue phenomenons_modelisations.py
# et dans le catalogue mesh_types.py
# Il NE FAUT PAS les redéfinir dans les catalogues.

# Pour certains de ces attributs : ALIAS8, MODELI,
# la liste des valeurs autorisées (donnée ci-dessous) n'est pas complète car cette liste serait
# trop fastidieuse à maintenir.
#
# Les commentaires (comment=...) sont très importants. Ce sont eux qui "définissent" les attributs.
# Pour l'instant, ces commentaires ne sont pas stockés dans les objets jeveux '&CATA....'
#
# --------------------------------------------------------------------------------------

DIM_TOPO_MAILLE = Attribute(
    auto=True,
    value=("0", "1", "2", "3"),
    comment="""
  DIM_TOPO_MAILLE : dimension topologique de la maille.
                     POI1 : 0 ; SEG: 1 ; TRIA/QUAD : 2 ; HEXA/PENTA/TETRA/PYRA:3
  Cet attribut est obtenu a partir du parametre "dim" du catalogue mesh_types

""",
)


PHENO = Attribute(
    auto=True,
    value=("AC", "ME", "PR", "TH"),
    comment="""
PHENO  =  'code' du phenomene
Si un element appartient a plusieurs phenomenes  : PHENO='##'
""",
)

MODELI = Attribute(
    auto=True,
    value=(
        "2FA",
        "2FL",
        "2FP",
        "3FA",
        "3FI",
        "3FL",
        "AFI",
        "AXF",
        "CL1",
        "CL2",
        "D2D",
        "D3D",
        "DIT",
        "FLS",
        "FS2",
        "FSA",
        "PFI",
        # ...
    ),
    comment="""
  MODELI :  'code' de la modelisaton
  Si un element appartient a plusieurs modelisations : MODELI='###'
""",
)

TYPMA = Attribute(
    auto=True,
    value=(
        "POI",
        "SE2",
        "TE4",
        # ...
    ),
    comment="""
TYPMA  =  'code' du type_maille
""",
)

ALIAS8 = Attribute(
    auto=True,
    value=("MEDTRSE2",),
    comment="""
  ALIAS8 : chaine formee par concatenation de 3 codes  :
           ALIAS(1:2) : code du phenomene
           ALIAS(3:5) : code de la modelisation
           ALIAS(6:8) : code du type_maille

           exemple :   MEFL_HEXA8 => ALIAS8=ME3FLHE8
""",
)

DIM_COOR_MODELI = Attribute(
    auto=True,
    value=("2", "3"),
    comment="""
  DIM_COOR_MODELI : dimension de l'espace du maillage
                    (nombre de coordonnees des points de l'espace) : 2 ou 3
  Cet attribut est obtenu a partir du parametre dim(2) de la modelisation
  remarque : DIM_COOR_MODELI >=  DIM_TOPO_MODELI
""",
)


DIM_TOPO_MODELI = Attribute(
    auto=True,
    value=("0", "1", "2", "3"),
    comment="""
  DIM_TOPO_MODELI : dimension "topologique" de la modelisation a laquelle appartient
  l'element. Exemples :
    - 3 pour les elements isoparametriques 3D
    - 2 pour les elements isoparametriques 2D
    - 2 pour plaques et coques en 3D
    - 1 pour les poutres ou pour les coques en 2D
    - ...
  Cet attribut est obtenu a partir du parametre dim(1) de la modelisation
  Convention importante : Pour les modelisations "discretes" ([2D_]DIS_xxx),
   on choisit DIM_TOPO_MODELI=-1
""",
)


PRINCIPAL = Attribute(
    auto=True,
    value=("OUI",),
    comment="""
  PRINCIPAL = 'OUI' :
  L'element est "principal" pour la modelisation (i.e. ce n'est pas un element de "bord").
  La dimension toplogique de sa maille (DIM_TOPO_MAILLE) est identique a celle de sa
  modelisation (DIM_TOPO_MODELI)
  Pour les modelisations "discretes" ([2D_]DIS_xxx), TOUS les elements sont principaux.
""",
)


BORD = Attribute(
    auto=True,
    value=("0", "-1", "-2"),
    comment="""
  / '0'  : l'element est principal : DIM_TOPO_MAILLE = DIM_TOPO_MODELI
  / '-1' : l'element est tel que   : DIM_TOPO_MAILLE = DIM_TOPO_MODELI - 1
  / '-2' : l'element est tel que   : DIM_TOPO_MAILLE = DIM_TOPO_MODELI - 2
  / '-3' : l'element est tel que   : DIM_TOPO_MAILLE = DIM_TOPO_MODELI - 3

  Un element de "bord" ('-1' et '-2') sert en general a appliquer des efforts sur l'element principal qu'il
  borde. Il fait partie de la meme modelisation que l'element principal.
  Il n'a pas de sens sans un element principal a cote de lui.
  Ses noeuds sont en principe ceux d'une face (ou d'une arrete) de son element principal.
  Il ne calcule pas de "rigidite".

  Un element de coque "plaque" sur un element de volume n'est pas un element de bord.
""",
)


DISCRET = Attribute(
    auto=True,
    value=("OUI",),
    comment="""
  DISCRET = 'OUI'
  La modelisation de l'element est "discrete" ([2D_]DIS_xxx)
  On reconnait les modelisations "discretes" au fait que DIM__(1)=-1
  Les elements d'une modelisations discrete sont TOUS "principaux" (voir PRINCIPAL).
""",
)
# --------------------------------------------------------------------------------------


ABSO = Attribute(
    value=("OUI",),
    comment="""
  ABSO = 'OUI' : l'element est un element de frontiere absorbante (R4.02.05)
""",
)

AXIS = Attribute(
    value=("OUI",),
    comment="""
  AXIS =  'OUI' : l'element est axisymetrique.
          C'est a dire qu'il est maille en 2D mais qu'il represente en realite un "anneau" de
          revolution autour de l'axe OY.
          Le champ de deplacement dans un plan radial est independant de l'azimut (theta).

""",
)


BORD_ISO = Attribute(
    value=("OUI",),
    comment="""
  BORD_ISO = 'OUI' => l'element est l'element de bord "standard" partage par les elements isoparametriques 2D ou 3D.
""",
)

CABLE = Attribute(
    value=("OUI",),
    comment="""
  CABLE =  'OUI' :  l'element est un element de cable (sans rigidité de compression)
""",
)

CL_DUAL = Attribute(
    value=("OUI",),
    comment="""
CL_DUAL='OUI' : l'element est un element utilise pour dualiser les conditions aux limites.
""",
)


CONTACT = Attribute(
    value=("OUI",),
    comment="""
  CONTACT = 'OUI' : l'element est utilise (en sous-terrain) pour la mise en oeuvre du contact
""",
)


COQUE = Attribute(
    value=("OUI",),
    comment="""
  COQUE  =  'OUI' :  l'element est un element de coque, de plaque, de grille, ... (element de structure surfacique).
""",
)


C_PLAN = Attribute(
    value=("OUI",),
    comment="""
  C_PLAN =  'OUI' : l'element est en "contraintes planes".
          Le maillage est 2D (OXY). L'etat de contraintes est "plan" (SIZZ=SIZX=SIYZ=0)
          L'hypothese des contraintes planes est a priori verifiee pour les structure minces selon OZ (plaque).

""",
)


D_PLAN = Attribute(
    value=("OUI",),
    comment="""
  D_PLAN =  'OUI' : l'element est en "deformations planes".
          Le maillage est 2D (OXY). Le deplacement est independant de Z et le deplacement suivant OZ est nul.
          Un element 2D represente en realite un solide cylindrique infini selon OZ.

""",
)


EFGE = Attribute(
    value=("OUI",),
    comment="""
  EFGE = 'OUI' : l'element mecanique est un element de structure connaissant les efforts generalises.
""",
)


EULER = Attribute(
    value=("OUI",),
    comment="""
  EULER =  'OUI' : l'element est un element de poutre (theorie d'Euler-Bernoulli).
""",
)


FOURIER = Attribute(
    value=("OUI",),
    comment="""
  FOURIER =  'OUI' : l'element est destine a une etude par decomposition en modes de fourier.
          Le maillage est 2D (OXY) (comme pour AXIS='OUI').

""",
)

FLUIDE = Attribute(
    value=("OUI", "NON"),
    comment="""
  FLUIDE =  'OUI' : l'element est destine a une etude d'un fluide
""",
)

FORMULATION = Attribute(
    value=(
        "HHO_CSTE",
        "HHO_LINE",
        "HHO_MLINE",
        "HHO_QUAD",
        "HHO_MQUAD",
        "HHO_CUBI",
        "HHO_MCUBI",
        "HHO_QUAR",
        "HHO_MQUAR",
        "U_P_PHI",
        "U_P",
        "U_PSI",
        "DIL",
        "DIL_INCO",
        "STA",
        "STA_INCO",
    ),
    comment="""
  FORMULATION =  'HHO_CSTE'  : formulation constante for HHO (0/0/0)
  FORMULATION =  'HHO_LINE'  : formulation linear for HHO (1/1/1)
  FORMULATION =  'HHO_MLINE'  : formulation linear for HHO (1/2/1)
  FORMULATION =  'HHO_QUAD'  : formulation quadratic for HHO (2/2/2)
  FORMULATION =  'HHO_MQUAD'  : formulation quadratic for HHO (2/3/2)
  FORMULATION =  'HHO_CUBI'  : formulation cubic for HHO (3/3/3)
  FORMULATION =  'HHO_MCUBI'  : formulation cubic for HHO (3/4/3)
  FORMULATION =  'HHO_QUAR'  : formulation order 4 for HHO (4/4/4)
  FORMULATION =  'HHO_MQUAR'  : formulation order 4 for HHO (4/5/4)
  FORMULATION =  'U_P_PHI' : formulation displacement/pressure/disp potential
  FORMULATION =  'U_P'    : formulation displacement/pressure
  FORMULATION =  'U_PSI'  : formulation displacement/speed potential
  FORMULATION =  'DIL'       : formulation 'old' for DIL
  FORMULATION =  'DIL_INCO'  : formulation incompressible for DIL
  FORMULATION =  'STA'       : formulation standard for MIXTE_STAB
  FORMULATION =  'STA_INCO'  : formulation incompressible for MIXTE_STAB
""",
)

FROTTEMENT = Attribute(
    value=("OUI",),
    comment="""
  FROTTEMENT = 'OUI' => l'element est utilise pour traiter le frottement.
""",
)

FSI = Attribute(
    value=("OUI",),
    comment="""
  FSI =  'OUI' : l'element est destine a une etude IFS
""",
)

GRAND_GLIS = Attribute(
    value=("OUI",),
    comment="""
    ...
""",
)


GRILLE = Attribute(
    value=("OUI",),
    comment="""
  GRILLE =  'OUI' :  l'element est un element de "grille".
""",
)

HHO = Attribute(
    value=("OUI",),
    comment="""
  HHO = 'OUI'.
""",
)


INCO = Attribute(
    value=("C2", "C2O", "C3", "C5GV"),
    comment="""
  INCO   : Type d'elements "incompressibles" (R3.06.08)
           - Le volume de l'element reste constant (trace(epsilon)=0)
           - Il existe des ddls qui "dualisent" la condition precedente (PRES,GONF)
  INCO =  Cxk
          x  = nombre de champ : 3 champs UPG ou 2 champs UP
          k  = stabilisation de l'element : O (OSGS)
               version ne respectant pas la LBB : B
""",
)


INTERFACE = Attribute(
    value=("OUI",),
    comment="""
  INTERFACE = 'OUI' : elements "ecrases" (pour representer un joint par exemple)
           En 3D ce sont des HEXA ou des PENTA
           En 2D ce sont des QUAD
""",
)


INTTHM = Attribute(
    value=("LUM", "RED"),
    comment="""
  INTTHM =  "type" d'integration pour les elements "principaux" des modelisations "THM" :
       /  'RED'      modelisation selective (termes capacitifs integres aux sommets, diffusifs aux points de Gauss)
       /  'LUM'      modelisation lumpee (integration aux sommets)
       Par defaut : INTTHM='CLA' ("classique")
""",
)


LUMPE = Attribute(
    value=("OUI",),
    comment="""
  LUMPE = 'OUI' : elements "lumpes" (R3.06.07):
           - Ce sont des elements de thermique
           - Leur matrice de "masse" (MASS_THER) est diagonale
""",
)


LXFEM = Attribute(
    value=("OUI",),
    comment="""
  LXFEM = 'OUI' : modelisations "XFEM"
          Ce sont des modelisations pour lesquelles on ajoute des ddls d'enrichissement.
          Ce sont des modelisations pour lesquelles, la notion de bord est modifiee puisque
          des fissures peuvent traverser des elements "volumiques".
""",
)


METH_CONTINUE = Attribute(
    value=("OUI",),
    comment="""
    ...
""",
)


NBSIGM = Attribute(
    value=("4", "6"),
    comment="""
  NBSIGM =  /'4'/'6'/...
        Pour les elements de mecanique, c'est le nombre de composantes du tenseur des contraintes.
        En general : '6' pour les elements 3D, '4' pour certains elements 2D
""",
)

PESA = Attribute(
    value=("OUI",),
    comment="""
  PESA =  'OUI' : l'element est destine a une etude fluide à surface libre
""",
)

PLAQUE = Attribute(
    value=("OUI",),
    comment="""
  PLAQUE  =  'OUI' :  l'element est un element de plaque
""",
)


POUTRE = Attribute(
    value=("OUI",),
    comment="""
  POUTRE =  'OUI' :  l'element est un element de poutre, de barre, de cable, de tuyau, ... (element de structure lineique).
""",
)

POUX = Attribute(
    value=("OUI", "NON"),
    comment="""
  POUX  =  'OUI' :  l'élément est de type poutre 'à la POUX' (RDM)
""",
)

RACCORD = Attribute(
    value=("OUI",),
    comment="""
  RACCORD = 'OUI' : l'element est utilise (en sous-terrain) pour la mise en oeuvre de raccord.
""",
)

SECH = Attribute(
    value=("OUI",),
    comment="""
  SECH =  'OUI' : l'element est un élément thermique spécifique du traitement
          du séchage (via SECH_NON_LINE).

""",
)

SIGM = Attribute(
    value=("NON",),
    comment="""
  SIGM = 'NON' : l'element mecanique est un element de structure connaissant les contraintes (SIXX, SIYY, ...)
          Remarque importante : l'attribut SIGM='NON' ne doit etre renseigne que pour les elements ayant
                                l'attribut EFGE='OUI'
""",
)


SOUS_POINT = Attribute(
    value=("OUI",),
    comment="""
  SOUS_POINT = 'OUI' => l'element peut definir des sous-points dans AFFE_CARA_ELEM
""",
)

STRX = Attribute(
    value=("OUI", "NON"),
    comment="""
  STRX : indique si l'élément a besoin du champ STRX
""",
)

TUYAU = Attribute(
    value=("OUI",),
    comment="""
  TUYAU  =  'OUI' :  l'element est un element de "tuyau".
""",
)

MECA = Attribute(
    value=("OUI", "NON"),
    comment="""
  MECA  =  'OUI' :  l'element est en mécanique.
""",
)

THER = Attribute(
    value=("OUI", "NON"),
    comment="""
  THER  =  'OUI' :  l'element est en thermique.
""",
)

HYDR1 = Attribute(
    value=("0", "1", "2"),
    comment="""
 Nombre de phases pour le premier constituant
""",
)

HYDR2 = Attribute(
    value=("0", "1", "2"),
    comment="""
 Nombre de phases pour le second constituant
""",
)

DIL = Attribute(
    value=("OUI", "NON"),
    comment="""
  DIL  =  'OUI' :  l'element est en second gradient de dilatation.
""",
)

TYPE_VOISIN = Attribute(
    value=("A2", "F3"),
    comment="""
  TYPE_VOISIN = typvoi : type du voisinage des elements de type "volume_fini"
                la chaine de caracteres typvoi est formee d'une suite de K2 tels que :
                  'F3' : voisin 3D par une face
                  'A2' : voisin 2D par une arrete
                   ...   (voir routine voiuti.F90)
""",
)


TYPMOD = Attribute(
    value=("AXIS", "0D", "1D", "3D", "C_PLAN", "D_PLAN", "PLAN", "MEMBRANE"),
    comment="""
  TYPMOD : Type de  modelisation utilise pour integrer les lois de comportement (utilise dans NMCOMP)
           Les valeurs possibles de cet attribut sont definies dans le catalogue de chaque comportement
           (bibpyt/Comportement/*.py)
""",
)


TYPMOD2 = Attribute(
    value=(
        "EJ_HYME",
        "ELEMJOIN",
        "GRADVARI",
        "GRADSIGM",
        "INTERFAC",
        "INTSOLPI",
        "CABLE_GA",
        "PMF",
        "THM",
        "HHO",
        "HHO_GRAD",
        "JHMS",
        "XFEM_HM",
    ),
    comment="""
  TYPMOD2 : Complement au type de  modelisation utilise pour integrer les lois de comportement TYPMOD
            (TYPMOD2 est utilise dans NMCOMP et certaines routines LCxxxx)
           Les valeurs possibles de cet attribut sont :
  TYPMOD2= GRADVARI l'element utilise des comportements non locaux (TYPMOD(2)='GRADVARI' pour NMCOMP)
           ELEMJOIN l'element utilise des comportements d'elements de joints (CZM sur des modelisations *_JOINT)
           INTERFAC l'element utilise des comportements d'elements d'interface (CZM sur des modelisations *_INTERFACE)
           INTSOLPI l'element utilise des comportements d'elements d'interface sol-pieu (SP_* sur une modelisation 3D_SOLPIEUINT)
           PMF      l'element fait appel a des comportements 1D PMF
           THM      thermo-hydro-mechanic
           HHO      Hybrid High-Order elements
           HHO_GRAD = HHO + GRADVARI
           JHMS     Hydraulic joints
           XFEM_HM  XFEM with HM (and contact, sometimes ! )
""",
)

TYPMOD3 = Attribute(
    value=("SUSHI",),
    comment="""
  TYPMOD3 : Complement au type de  modelisation utilise pour integrer les lois de comportement TYPMOD
            (TYPMOD3 est utilise dans NMCOMP et certaines routines LCxxxx)
           Les valeurs possibles de cet attribut sont :
  TYPMOD3= SUSHI    finite volume
""",
)


XFEM = Attribute(
    value=("XH", "XH1", "XH2", "XH3", "XH4", "XHT", "XT"),
    comment="""
  XFEM = 'XH','XT','XHT'     : element XFEM de type Heaviside, cracktip, ou mixte
""",
)


XFEM_E = Attribute(
    value=("C", "H", "H2", "H3", "H4", "T"),
    comment="""
  XFEM_E = 'H','T','C','H2','H3','H4' : element XFEM maille "esclave" de type Heaviside, cracktip, ou mixte

""",
)


XFEM_M = Attribute(
    value=("C", "H", "H2", "H3", "H4", "T"),
    comment="""
  XFEM_M = 'H','T','C','H2','H3','H4' : element XFEM maille "maitre" de type Heaviside, cracktip, ou mixte
""",
)


XLAG = Attribute(
    value=("NOEUD",),
    comment="""
  XLAG = 'NOEUD'       : element XFEM ????
""",
)


# store all Attribute objects
ATTRS = objects_from_context(globals(), Attribute)
