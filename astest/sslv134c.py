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

"""
Create a 2d mesh of a plate with a crack.
crack on the left side: CDAI
crack on the right side: CDAE

This mesh can also be used as a 2d-axi mesh of a tube with an
axisymetrical circonferential crack, located at internal or external
surface (config CDAI or CDAE of the RSE-M).

Features of the mesh:
- only the half-plane x-y, y>=0 is modeled
- semi-circular radial mesh around the crack tip
- the area at the top of the mesh has a grid mesh
- the intermediate area has a free mesh

Optional features:
- mesh of a section of a tube (some group names are changed)
- entire model (entire plate/tube), without the symetry

Warning:
- Groups building is base on the ids! (not robust if modification)
"""

import math

import salome
from salome.geom import geomBuilder

import SMESH
from salome.smesh import smeshBuilder


# ----------------------------------------------------------------------
# parametres utilisateur
# ----------------------------------------------------------------------

# géométrie
# ---------

# géométrie : 'plaque' ou 'tube'
geom = "plaque"

# prise en compte de la symétrie : True or False (True par défaut)
isSym = True

t = 10.0  # thickness of the plate/pipe (mm)
rm_t = None  # Only for tube: ratio of the mean radius over the thickness
a_t = 0.2  # ratio of the length of the crack over the thickness (0<a/t<1)
config = "CDAI"  # config ('CDAI' or 'CDAE')
H = 10.0  # (demi)hauteur de la plaque/tube. Si 'None', on utilise la hauteur
# liées aux distances d'amortissement

# maillage
# --------

rtore = a / 10

# nombre d'éléments le long de la (demi-)circonference [minimum 2]
nb_circ = 16

# choix du type de distribution radiale : uniforme ou géométrique
type_distr_radial_seg = "uniforme"
# ~type_distr_radial_seg = 'geometrique'

# nombre de couches dans le tore [minimum 2]
# Rq : ne sert pas si distribution géométrique
nb_layers = 5

# nombre d'éléments dans l'épaisseur de la plaque/tube (vers le haut)
nb_elt_ep = 15

# nombre d'élements dans la hauteur (partie réglée)
nb_elt_L = 12

# distribution suivant la hauteur
# ~distr = 'uniforme'
distr = "linear"
# ~distr = 'power' # conseillé pour les tubes

# type of elements ('quadratic' or 'linear')
elt_type = "quadratic"
# ~elt_type = 'linear'

# export in a MED file?
export = False

# ----------------------------------------------------------------------
# Initialization
# ----------------------------------------------------------------------

if geom == "plaque":
    rm_t = 1.0  # valeur bidon
rm = rm_t * t
ri = rm - t / 2
re = rm + t / 2
a = a_t * t
if config == "CDAE":
    xa = re - a
if config == "CDAI":
    xa = ri + a

# distance amortissement des contraintes de flexions locales
d1 = 5 * math.sqrt(rm * t)

# distance amortissement de la propagation de l'ovalisation
d2 = 1.5 * math.sqrt(rm**3 / t)

# length of the plate/pipe
if not H:
    H = max(d1, d2)

# lenght of the zone libre
if geom == "plaque":
    Hlibre = t / 4.0
elif geom == "tube":
    Hlibre = t

# taille moyenne des mailles dans la zone libre
h_moy = t / nb_elt_ep

if export:
    dir_name = r"/home/D24752/data/"
    file_name = "{}-{}.mmed".format(geom, config)

# ----------------------------------------------------------------------
# Building the Geometry
# ----------------------------------------------------------------------

salome.salome_init()
geompy = geomBuilder.New()

# section globale
Face = geompy.MakeFaceHW(t, H, 1)
geompy.TranslateDXDYDZ(Face, rm, H / 2, 0)

# tore (disque en 2d)
Disk = geompy.MakeDiskR(rtore, 1)
geompy.TranslateDXDYDZ(Disk, xa, 0, 0)
vertex1 = geompy.MakeVertex(xa, 0, 0)

# ligne de démarquation entre zone libre en bas et zone réglée en haut
vertex2 = geompy.MakeVertex(ri, Hlibre, 0)
vertex3 = geompy.MakeVertex(re, Hlibre, 0)
Line = geompy.MakeLineTwoPnt(vertex2, vertex3)

Plaque = geompy.MakePartition(
    [Face], [Disk, Line, vertex1], [], [], geompy.ShapeType["FACE"], 0, [], 0
)

# création des groupes [a refaire mieux]
#  -> il faut à la rigueur maj les numero d'ids
tore = geompy.CreateGroup(Plaque, geompy.ShapeType["FACE"])
zone_libre = geompy.CreateGroup(Plaque, geompy.ShapeType["FACE"])
zone_regle = geompy.CreateGroup(Plaque, geompy.ShapeType["FACE"])
plaque = geompy.CreateGroup(Plaque, geompy.ShapeType["FACE"])

cercle = geompy.CreateGroup(Plaque, geompy.ShapeType["EDGE"])
lig = geompy.CreateGroup(Plaque, geompy.ShapeType["EDGE"])
lev = geompy.CreateGroup(Plaque, geompy.ShapeType["EDGE"])
peau_int = geompy.CreateGroup(Plaque, geompy.ShapeType["EDGE"])
peau_ext = geompy.CreateGroup(Plaque, geompy.ShapeType["EDGE"])
haut = geompy.CreateGroup(Plaque, geompy.ShapeType["EDGE"])

fond = geompy.CreateGroup(Plaque, geompy.ShapeType["VERTEX"])

geompy.UnionIDs(tore, [23])
geompy.UnionIDs(zone_libre, [2])
geompy.UnionIDs(zone_regle, [16])
geompy.UnionIDs(plaque, [2, 16, 23])
geompy.UnionIDs(cercle, [13])
l_ids_int = [15, 27]
l_ids_ext = [11, 25]
if config == "CDAE":
    geompy.UnionIDs(lig, l_ids_int)
    geompy.UnionIDs(lev, l_ids_ext)
elif config == "CDAI":
    geompy.UnionIDs(lig, l_ids_ext)
    geompy.UnionIDs(lev, l_ids_int)
geompy.UnionIDs(peau_int, [4, 18])
geompy.UnionIDs(peau_ext, [9, 22])
geompy.UnionIDs(haut, [20])
geompy.UnionIDs(fond, [26])

[Compound_1, Compound_2] = geompy.Propagate(zone_regle)

if geom == "plaque":
    geompy.TranslateDXDYDZ(Plaque, -t / 2.0, 0, 0)


geompy.addToStudy(Plaque, "Plaque")
geompy.addToStudyInFather(Plaque, plaque, "plaque")
geompy.addToStudyInFather(Plaque, tore, "tore")
geompy.addToStudyInFather(Plaque, zone_libre, "zone_libre")
geompy.addToStudyInFather(Plaque, zone_regle, "zone_regle")
geompy.addToStudyInFather(Plaque, cercle, "cercle")
geompy.addToStudyInFather(Plaque, lig, "lig")
geompy.addToStudyInFather(Plaque, lev, "lev_sup")
geompy.addToStudyInFather(Plaque, peau_int, "peau_int")
geompy.addToStudyInFather(Plaque, peau_ext, "peau_ext")
geompy.addToStudyInFather(Plaque, haut, "haut")
geompy.addToStudyInFather(Plaque, fond, "fond")
geompy.addToStudyInFather(Plaque, Compound_1, "Compound_1")
geompy.addToStudyInFather(Plaque, Compound_2, "Compound_2")

# ----------------------------------------------------------------------
# Building the Mesh
# ----------------------------------------------------------------------

smesh = smeshBuilder.New()

mesh = smesh.Mesh(Plaque)
mesh.Segment(cercle).NumberOfSegments(nb_circ)

# hypothèse pour le maillage radial
algo_radial = mesh.Quadrangle(algo=smeshBuilder.RADIAL_QUAD, geom=tore)
if type_distr_radial_seg == "uniforme":
    algo_radial.NumberOfLayers(nb_layers)
elif type_distr_radial_seg == "geometrique":
    # Progression geometrique
    LayersDistr = smesh.CreateHypothesis("LayerDistribution2D")
    GeomProgDistr = smesh.CreateHypothesis("GeometricProgression")
    GeomProgDistr.SetStartLength(rtore * math.pi / nb_circ)
    ratio = 1 / math.exp(math.pi / nb_circ)
    GeomProgDistr.SetCommonRatio(ratio)
    LayersDistr.SetLayerDistribution(GeomProgDistr)
    status = mesh.AddHypothesis(LayersDistr)


# hypothèse pour la partie libre
algo_libre = mesh.Triangle(algo=smeshBuilder.MG_CADSurf, geom=zone_libre)
algo_libre_para = algo_libre.Parameters()
algo_libre_para.SetPhySize(h_moy)
algo_libre_para.SetMinSize(rtore * math.pi / nb_circ)
algo_libre_para.SetMaxSize(h_moy)
# triangles (0), quadrangles dominant (1) ou quadrangles (2)
algo_libre_para.SetElementType(1)

# hypothèse pour la partie réglée
algo_regle = mesh.Quadrangle(algo=smeshBuilder.QUADRANGLE, geom=zone_regle)
Regular_1D_1 = mesh.Segment(geom=Compound_2)
Regular_1D_2 = mesh.Segment(geom=Compound_1)

H1 = Regular_1D_1.NumberOfSegments(nb_elt_ep)  # dans l'épaisseur
H2 = Regular_1D_2.NumberOfSegments(nb_elt_L)  # suivant la hauteur
if distr == "uniforme":
    pass
elif distr == "linear":
    # on donne la taille (en relatif) au départ (ligne 0)
    # et à l'arrivée (ligne 1)
    H2.SetTableFunction([0, 6.0, 1, 1.0])
elif distr == "power":
    # length of segments gradually changes depending on the Scale Factor
    # (ratio first / last segment)
    # Length of segments changes in geometric progression A = S**(1/(N-1))
    # For an edge of length L, length of the first segment is L*(1-A)/(1-A**N)
    H2.SetScaleFactor(6)


Sub_mesh_1 = algo_radial.GetSubMesh()
Sub_mesh_2 = algo_libre.GetSubMesh()
Sub_mesh_3 = algo_regle.GetSubMesh()

mesh.SetMeshOrder([[Sub_mesh_1, Sub_mesh_2, Sub_mesh_3]])
mesh.Compute()

Gr_plaque = mesh.GroupOnGeom(plaque, "PLAQUE", SMESH.FACE)
Gr_tore = mesh.GroupOnGeom(tore, "TORE", SMESH.FACE)
Gr_zone_libre = mesh.GroupOnGeom(zone_libre, "LIBRE", SMESH.FACE)
Gr_zone_regle = mesh.GroupOnGeom(zone_regle, "REGLE", SMESH.FACE)
Gr_lig = mesh.GroupOnGeom(lig, "LIG", SMESH.EDGE)
Gr_lev = mesh.GroupOnGeom(lev, "LEV_SUP", SMESH.EDGE)
if geom == "plaque":
    nameG = "GAUCHE"
    nameD = "DROITE"
elif geom == "tube":
    nameG = "PEAU_INT"
    nameD = "PEAU_EXT"
Gr_peau_int = mesh.GroupOnGeom(peau_int, nameG, SMESH.EDGE)
Gr_peau_ext = mesh.GroupOnGeom(peau_ext, nameD, SMESH.EDGE)
Gr_haut = mesh.GroupOnGeom(haut, "HAUT", SMESH.EDGE)
Gr_fond = mesh.GroupOnGeom(fond, "FRONT", SMESH.NODE)

# TRAITEMENT DU CAS DE LA PLAQUE/TUBE COMPLET
if not isSym:
    [
        plaque_sym,
        tore_sym,
        zone_libre_sym,
        zone_regle_sym,
        lig_sym,
        lev_sym,
        peau_int_sym,
        peau_ext_sym,
        haut_sym,
        fond_sym,
    ] = mesh.MirrorObject(
        mesh, SMESH.AxisStruct(0, 0, 0, 0, 1, 0), SMESH.SMESH_MeshEditor.PLANE, True, True
    )

    coincident_nodes = mesh.FindCoincidentNodesOnPart([mesh], 1e-05, [Gr_lev, lev_sym], 0)
    mesh.MergeNodes(coincident_nodes)

    coincident_nodes = mesh.FindCoincidentNodesOnPart([Gr_fond, fond_sym], 1e-05, [], 0)
    mesh.MergeNodes(coincident_nodes)

    equal_elements = mesh.FindEqualElements([mesh])
    mesh.MergeElements(equal_elements)

    # Ménage dans les groupes
    lev_sym.SetName("LEV_INF")
    haut_sym.SetName("BAS")
    new_grto = mesh.GetMesh().UnionListOfGroups([Gr_tore, tore_sym], "TORE")
    new_grpi = mesh.GetMesh().UnionListOfGroups([Gr_peau_int, peau_int_sym], nameG)
    new_grpe = mesh.GetMesh().UnionListOfGroups([Gr_peau_ext, peau_ext_sym], nameD)

    gr_to_remove = (
        lig_sym,
        Gr_plaque,
        Gr_zone_libre,
        Gr_zone_regle,
        plaque_sym,
        zone_libre_sym,
        zone_regle_sym,
        fond_sym,
        Gr_tore,
        tore_sym,
        Gr_peau_int,
        peau_int_sym,
        Gr_peau_ext,
        peau_ext_sym,
        Gr_lig,
    )
    for gr in gr_to_remove:
        mesh.RemoveGroup(gr)


# Set names of Mesh objects
smesh.SetName(algo_radial.GetAlgorithm(), "algo_radial")
smesh.SetName(algo_libre.GetAlgorithm(), "algo_libre")
smesh.SetName(algo_libre_para, "algo_libre_para")
smesh.SetName(Sub_mesh_1, "Sub-mesh_1")
smesh.SetName(Sub_mesh_2, "Sub-mesh_2")
smesh.SetName(mesh.GetMesh(), geom)

assert elt_type in ("linear", "quadratic")
if elt_type == "quadratic":
    mesh.ConvertToQuadratic(0)

# affichage dans l'arbre
if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()


# ----------------------------------------------------------------------
# Export the File
# ----------------------------------------------------------------------

# export in a MED file
if export:
    try:
        mesh.ExportMED(dir_name + file_name, 0, SMESH.MED_V2_2, 1, None, 1)
        pass
    except:
        print("Echec dans l'export du maillage. Vérifier les chemins")
