# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
# aslint: disable=C4008

import itertools as it
import os
import re
import shutil
import sys
from types import MappingProxyType
from typing import Any, ClassVar, Sequence, Literal, Mapping, Optional, Union

import numpy as np
from numpy.ma import mrecords as mr
from numpy.lib import recfunctions as rfn

from code_aster.Messages import UTMESS
from code_aster.Utilities import cut_long_lines, is_float, is_sequence
from code_aster.Utilities.misc import fmtF2PY

K_DEF = "K24"  # type par défaut des chaines de caractères

StrOrSeqStr = Union[str, Sequence[str]]
OrderType = Literal["CROISSANT", "DECROISSANT"]
FormatTable = Literal["TABLEAU", "ASTER", "NUMPY", "XMGRACE", "AGRAF", "TABLEAU_CROISE"]


class TableBase:
    """Classe pour partager les méthodes d'impression entre Table et Colonne
    (c'est surtout utile pour vérifier que l'extraction et les filtres sur les
    colonnes sont corrects).
    """

    DicForm: ClassVar[dict[str, Any]] = {
        "chead": None,  # entête avant la table
        "cfoot": "",  # lignes après la table
        "csep": " ",  # séparateur
        "ccom": "#",  # commentaire
        "ccpara": "",  # commentaire des labels
        "cdeb": "",  # début de ligne
        "cfin": "\n",  # fin de ligne
        "sepch": ";",  # séparateur entre deux lignes d'une cellule
        "formK": "%-12s",  # chaines
        "formR": "%12.5E",  # réels
        "formI": "%12d",  # entiers
    }

    data: np.ma.MaskedArray

    @property
    def para(self) -> list[str]:
        raise NotImplementedError

    @property
    def titr(self) -> str:
        raise NotImplementedError

    def Croise(self, **kwargs:  Any) -> "TableBase":
        raise NotImplementedError("Must be defined in a derived class")

    def Impr(
        self,
        FICHIER: Optional[str] = None,
        FORMAT: FormatTable = "TABLEAU",
        dform: Optional[dict[str, Any]] = None,
        **opts: Any,
    ):
        """Impresssion de la Table selon le format spécifié.
        FICHIER : nom du(des) fichier(s). Si None, on dirige vers stdout
        dform : dictionnaire de formats d'impression (format des réels,
           commentaires, saut de ligne...)
        opts  : selon FORMAT.
        """
        para = {
            "TABLEAU": {"mode": "a", "driver": self.ImprTableau},
            "ASTER": {"mode": "a", "driver": self.ImprTableau},
            "NUMPY": {"mode": "a", "driver": self.ImprNumpy},
            "XMGRACE": {"mode": "a", "driver": self.ImprGraph},
            "AGRAF": {"mode": "a", "driver": self.ImprTableau},
            "TABLEAU_CROISE": {"mode": "a", "driver": self.ImprTabCroise},
        }
        kwargs = {
            "FICHIER": FICHIER,
            "FORMAT": FORMAT,
            "dform": self.DicForm.copy(),
            "mode": para[FORMAT]["mode"],
        }
        if dform is not None and type(dform) is dict:
            kwargs["dform"].update(dform)
        # ajout des options
        kwargs.update(opts)
        if not kwargs.get("PAGINATION"):
            # call the associated driver
            para[FORMAT]["driver"](**kwargs)
        else:
            pagination_params = (
                list(kwargs["PAGINATION"])
                if is_sequence(kwargs["PAGINATION"])
                else [kwargs["PAGINATION"]]
            )
            del kwargs["PAGINATION"]
            npag = len(pagination_params)
            # paramètres hors ceux de la pagination
            lkeep = [p for p in self.para if p not in pagination_params]
            # création des listes des valeurs distinctes
            lvd = []
            for p in pagination_params:
                lvp = getattr(self, p).values()
                lvn = []
                for it in lvp:
                    if it is not None and it not in lvn:
                        lvn.append(it)
                lvn.sort()
                lvd.append(lvn)
            # création des n-uplets
            s = "[[" + ",".join(["x" + str(i) for i in range(npag)]) + "] "
            s += (
                " ".join(["for x" + str(i) + " in lvd[" + str(i) + "]" for i in range(npag)])
                + "]"
            )
            try:
                lnup = eval(s)
            except SyntaxError as s:
                UTMESS("F", "TABLE0_20")
            # pour chaque n-uplet, on imprime la sous-table
            for nup in lnup:
                tab = self
                for i in range(npag):
                    tab = tab & (getattr(tab, pagination_params[i]) == nup[i])
                    sl = ""
                    if tab.titr:
                        sl = "\n"
                    tab.titr += sl + pagination_params[i] + ": " + str(nup[i])
                tab[lkeep].Impr(**kwargs)

    def ImprTableau(self, **kwargs: Any) -> None:
        """Impression au format TABLEAU ou ASTER"""
        filepath = kwargs.get("FICHIER")
        filestream = None
        try:
            # file or stdout
            if filepath is not None:
                filestream = open(filepath, kwargs["mode"])
            else:
                filestream = sys.stdout
            # write
            filestream.write(self.ReprTable(**kwargs) + "\n")
            filestream.flush()
        finally:
            # close the file
            if filestream is not None and filepath is not None:
                filestream.close()

    def ImprNumpy(self, **kwargs: Any) -> None:
        numeric_fields = [
            name
            for name in self.data.dtype.names
            if np.issubdtype(self.data.dtype[name], np.integer)
            or np.issubdtype(self.data.dtype[name], np.floating)
        ]
        if len(numeric_fields) == 0:
            UTMESS("F", "TABLE0_5", valk=repr(numeric_fields))
        tab = rfn.structured_to_unstructured(self.data.filled(0))
        filename = kwargs.get("FICHIER")
        if filename is None:
            print(tab)
        else:
            np.save(filename, tab)
            if not filename.endswith(".npy"):
                shutil.move(f"{filename}.npy", filename)

    def ImprGraph(self, **kwargs: Any) -> None:
        """Impression au format XMGRACE : via le module Graph"""
        from .table_graph import AjoutParaCourbe, Graph

        args = kwargs.copy()
        if len(self.para) != 2:
            UTMESS("A", "TABLE0_21")
            return
        # suppression des lignes contenant une cellule vide
        tnv = getattr(self, self.para[0]).NON_VIDE() & getattr(self, self.para[1]).NON_VIDE()
        # objet Graph
        graph = Graph()
        dicC = {
            "Val": [
                list(getattr(tnv, tnv.para[0]).values()),
                list(getattr(tnv, tnv.para[1]).values()),
            ],
            "Lab": tnv.para,
        }
        AjoutParaCourbe(dicC, args)
        graph.AjoutCourbe(**dicC)

        # Surcharge des propriétés du graphique et des axes
        # (bloc quasiment identique dans impr_fonction_ops)
        if args.get("TITRE"):
            graph.Titre = args["TITRE"]
        if args.get("BORNE_X"):
            graph.Min_X = args["BORNE_X"][0]
            graph.Max_X = args["BORNE_X"][1]
        if args.get("BORNE_Y"):
            graph.Min_Y = args["BORNE_Y"][0]
            graph.Max_Y = args["BORNE_Y"][1]
        if args.get("LEGENDE_X"):
            graph.Legende_X = args["LEGENDE_X"]
        if args.get("LEGENDE_Y"):
            graph.Legende_Y = args["LEGENDE_Y"]
        if args.get("ECHELLE_X"):
            graph.Echelle_X = args["ECHELLE_X"]
        if args.get("ECHELLE_Y"):
            graph.Echelle_Y = args["ECHELLE_Y"]
        if args.get("GRILLE_X"):
            graph.Grille_X = args["GRILLE_X"]
        if args.get("GRILLE_Y"):
            graph.Grille_Y = args["GRILLE_Y"]
        try:
            graph.Trace(**args)
        except TypeError:
            UTMESS("A", "TABLE0_22")

    def ImprTabCroise(self, **kwargs: Any) -> None:
        """Impression au format TABLEAU_CROISE d'une table ayant 3 paramètres."""
        tabc = self.Croise()
        kwargs["FORMAT"] = "TABLEAU"
        tabc.Impr(**kwargs)

    def ReprTable(self, FORMAT="TABLEAU", dform=None, **kwargs: Any) -> str:
        """Représentation d'une Table ou d'une Colonne sous forme d'un tableau."""
        rows = self.rows
        para = self.para
        typ = self.type
        if not is_sequence(para):
            para = [self.para]
            typ = [self.type]
        if dform is None:
            dform = TableBase.DicForm.copy()
        if dform["chead"] is None:
            dform["chead"] = os.linesep.join(
                [dform["ccom"], dform["ccom"] + "-" * 80, dform["ccom"]]
            )
        # est-ce que l'attribut .type est renseigné ?
        typdef = typ != [None] * len(typ)
        is_aster_format = FORMAT == "ASTER"
        lspa = []
        if is_aster_format:
            # ['']+ pour ajouter un séparateur en début de ligne
            lspa.append("")
        # lmax : largeur max des colonnes = max(form{K,R,I},len(parametre))
        lmax = []
        # formats
        strfmt, strfmt_none = {}, {}
        for t, p in zip(typ, para):
            larg_max = max(
                [len(str(p))]
                + [len(FMT(dform, k, t) % 0) for k in ("formK", "formR", "formI")]
            )
            lspa.append(FMT(dform, "formK", t, larg_max, str(p)) % p)
            lmax.append(larg_max)
            assert t is not None, "Type de la colonne '%s' non défini" % p
            strfmt[p] = FMT(dform, "form" + t[0], t, larg_max)
            strfmt_none[p] = FMT(dform, "formK", t, larg_max)
        if typdef:
            stype = dform["csep"].join(
                [""]
                + [FMT(dform, "formK", typ[i], lmax[i]) % typ[i] for i in range(len(para))]
            )
        txt = []
        if is_aster_format:
            txt.append("#DEBUT_TABLE")
        if self.titr:
            if is_aster_format:
                txt.extend(["#TITRE " + lig for lig in self.titr.split("\n")])
            else:
                txt.extend([dform["ccom"] + lig for lig in self.titr.split("\n")])
        if " " not in dform["csep"]:
            lspa = [i.strip() for i in lspa]
        txt.append(dform["ccpara"] + dform["csep"].join(lspa))
        if is_aster_format and typdef:
            txt.append(stype)
        for r in rows:
            lig = []
            if is_aster_format:
                lig.append("")
            empty = True
            for t, p, lmax_i in zip(typ, para, lmax):
                val = r.get(p)
                if val is not None:
                    empty = False
                    lig.append(strfmt[p] % val)
                else:
                    s = strfmt_none[p] % "-"
                    # format AGRAF = TABLEAU + '\' devant les chaines de
                    # caractères !
                    if FORMAT == "AGRAF":
                        s = "\\" + s
                    lig.append(s)
            if not empty:
                lig2 = [dform["sepch"].join(ch.splitlines()) for ch in lig]
                txt.append(dform["csep"].join(lig2))
        if is_aster_format:
            txt.append("#FIN_TABLE")
        # ajout des debut et fin de ligne
        txt = [dform["cdeb"] + t + dform["cfin"] for t in txt]
        txt.insert(0, dform["chead"])
        txt.append(dform["cfoot"])
        return "".join(txt)

    def __repr__(self) -> str:
        return self.ReprTable()


class Table(TableBase):
    """Table using a masked numpy recarray under the hood.

    The attribute `data` contains the numpy masked record array.
    The object can be iterate over the rows and returns a read only dictionnary of the column
    values.
    A print of the table return a fancy display of the table.
    Each columns of the table can be accessed as an attribute which returns a `Colonne` object
    without any copy of the data. The colonne can be used to extract information using basic
    operators. 
    """

    def __init__(
        self,
        rows: Optional[Sequence[Mapping[str, Any]]] = None,
        para: Optional[Sequence[str]] = None,
        typ: Optional[Sequence[str]] = None,
        titr: str = "",
        nom: str = "",
    ) -> None:
        para = para or []
        typ = typ or []
        for name in para:
            if para.count(name) != 1:
                UTMESS("F", "TABLE0_23", valk=name)
        if len(typ) == len(para):
            typ = list(typ)
        else:
            typ = [None] * len(para)
        # to keep the compatibility with previous version the None types are saved
        # to be able to specify it when trying to get the types
        self._none_types: set[int] = set()
        self._name: str = nom
        self._title: str = titr
        self.data: np.ma.MaskedArray = self._create_recarray(para, typ, rows)
        self.referenceToDataStructure: list[Any] = []

    @property
    def para(self) -> list[str]:
        return list(self.data.dtype.names)

    @property
    def type(self) -> list[str]:
        return [
            self._get_aster_type(name) if i not in self._none_types else None
            for i, name in enumerate(self.data.dtype.names)
        ]

    @property
    def rows(self) -> tuple[MappingProxyType[str, Any]]:
        return tuple(self)

    @property
    def nom(self) -> str:
        return self._name

    @property
    def titr(self) -> str:
        return self._title

    @titr.setter
    def titr(self, title: str) -> None:
        self._title = title

    @classmethod
    def from_masked_array(
        cls, masked_array: np.ma.MaskedArray, nom: str = "", titr: str = ""
    ) -> "Table":
        table = cls(nom=nom, titr=titr)
        table.data = masked_array
        return table

    def getName(self) -> str:
        """Get the name of the Table."""
        return self._name

    def get_value(self, column: str, index: int) -> Optional[Any]:
        """Get the value at the given column and index.

        Returns None if the value is missing.
        """
        value = self.data[column][index]
        if np.ma.is_masked(value):
            return None
        return value

    def get_row_values(self, index: int) -> dict[str, Any]:
        """Get the values of the row in a dict format."""
        values = {}
        data = self.data[index]
        for name in self.data.dtype.names:
            value = data[name]
            if not np.ma.is_masked(value):
                values[name] = value
        return values

    def get_max_by_column(self, column: str) -> Optional[Any]:
        """Get the maximum value of the column.

        If all the values of the column are missing then it returns None.
        """
        max_value = self.data[column].max()
        return None if max_value is np.ma.masked else max_value.item()

    def get_min_by_column(self, column: str) -> Optional[Any]:
        """Get the minimum value of the column.

        If all the values of the column are missing then it returns None.
        """
        min_value = self.data[column].min()
        return None if min_value is np.ma.masked else min_value.item()

    def copy(self) -> "Table":
        """Get a copy of the Table."""
        return Table.from_masked_array(self.data.copy(), nom=self._name, titr=self._title)

    def add_para(self, para: StrOrSeqStr, typ: StrOrSeqStr) -> None:
        """Add a new column to the table."""
        if not is_sequence(para):
            para = [para]
        if not is_sequence(typ):
            typ = [typ]
        new_names = []
        new_typ = []
        names = set(self.data.dtype.names)
        for name, t in zip(para, typ):
            if name not in names:
                new_names.append(name)
                new_typ.append(t)
                names.add(name)
        dtypes = [
            (col_name, self._get_numpy_type(col_typ, col_name))
            for col_name, col_typ in zip(new_names, new_typ)
        ]
        combined_dtypes = list(self.data.dtype.descr) + dtypes
        new_data = np.ma.masked_all(len(self.data), dtype=np.dtype(combined_dtypes))
        for name in self.data.dtype.names:
            new_data[name] = self.data[name]
        self.data = new_data

    def append(self, obj: dict[str, Any]) -> None:
        """Ajoute une ligne (type dict) qui peut éventuellement définir un
        nouveau paramètre."""
        self.extend([obj])

    def extend(self, objs: list[dict[str, Any]]) -> None:
        """Ajoute plusieurs lignes (list of dict)."""
        para = list(self.data.dtype.names)
        typ = []
        for i, (n, t) in enumerate(zip(para, self.type)):
            if i not in self._none_types:
                typ.append(t)
            else:
                for j in range(len(objs)):
                    try:
                        t = typaster(objs[j][n])
                    except KeyError:
                        pass
                    else:
                        typ.append(t)
                        break
        new_para = list({p for obj in objs for p in obj} - set(para))
        new_typ: list[str] = []
        for p in new_para:
            for obj in objs:
                if p in obj:
                    new_typ.append(typaster(obj[p]))
                    break
        if len(new_typ) > 0:
            self.add_para(new_para, new_typ)
        new_data = self._create_recarray(para + new_para, typ + new_typ, objs)
        if len(self.data):
            self.data = np.ma.concatenate((self.data, new_data))
        else:
            self.data = new_data

    def SansColonneVide(self, l_para: Optional[list[str]] = None) -> None:
        """Retourne une copie de la table dans laquelle on a supprimé les colonnes
        vides (les lignes vides sont automatiquement supprimées).
        """
        names = l_para or list(self.data.dtype.names)
        keep = [name for name in names if not np.all(self.data[name].mask)]
        return Table.from_masked_array(rfn.repack_fields(self.data[keep]), titr=self.titr)

    def fromfunction(
        self,
        nom_para: str,
        funct: Any,
        l_para: Optional[StrOrSeqStr] = None,
        const: Optional[dict[str, Any]] = None,
    ):
        """Ajoute une colonne `nom_para` en évaluant la fonction `funct` sur
        la valeur des paramètres `l_para` (qui doivent exister dans la table).
        Si `l_para` n'est pas fourni, on prend `funct`.getVariables() (FORMULE Aster).
        On peut passer un dictionnaire de constantes dans `const`. Quand on
        utilise une FORMULE Aster, les constantes sont prises dans le contexte
        global.
        """
        para = self.data.dtype.names
        if not hasattr(funct, "__call__"):
            UTMESS("F", "TABLE0_25", valk=(funct.getName(), "__call__"))
        if nom_para in para:
            UTMESS("F", "TABLE0_24", valk=nom_para)
        if l_para is None:
            if not hasattr(funct, "getVariables"):
                UTMESS("F", "TABLE0_25", valk=(funct.getName(), "getVariables"))
            l_para = funct.getVariables()
        if not is_sequence(l_para):
            l_para = [l_para]
        else:
            l_para = list(l_para)
        not_found = ", ".join([p for p in l_para if p not in para])
        if not_found != "":
            UTMESS("F", "TABLE0_27", valk=not_found)
        if const is None:
            const = {}
        if not isinstance(const, dict):
            UTMESS("F", "TABLE0_28", valk=("const", "dict"))
        res = [
            None if is_masked else funct(*values, **const)
            for values, is_masked in zip(self.data.data[l_para], self.data.mask[l_para])
        ]
        self[nom_para] = res

    def sort(self, CLES: Optional[StrOrSeqStr] = None, ORDRE: OrderType = "CROISSANT"):
        """Tri de la table.
        CLES  : liste des clés de tri
        ORDRE : CROISSANT ou DECROISSANT
        """
        para = list(self.data.dtype.names)
        # vérification des arguments
        if CLES is None:
            # par défaut, on prend tous les paramètres
            CLES = para
        else:
            if not is_sequence(CLES):
                CLES = [CLES]
            else:
                CLES = list(CLES)
            not_found = ", ".join([p for p in CLES if p not in para])
            if not_found != "":
                UTMESS("F", "TABLE0_27", valk=not_found)
        if ORDRE not in ("CROISSANT", "DECROISSANT"):
            UTMESS("F", "TABLE0_29", valk=ORDRE)
        idx = self.data.argsort(order=CLES, fill_value=0.0)
        if ORDRE == "DECROISSANT":
            idx = idx[::-1]
        self.data = rfn.repack_fields(self.data[idx])

    def OrdreColonne(self, cols: StrOrSeqStr) -> None:
        """Sort the columns putting first the given columns.

        Ignore the columns not in the table.
        """
        if type(cols) not in (list, tuple):
            cols = [cols]
        current_names = set(self.data.dtype.names)
        new_names = [col for col in cols if col in current_names]
        new_names.extend(name for name in self.data.dtype.names if name not in cols)
        self.data = rfn.repack_fields(self.data[new_names])

    def difference(self, other):
        """Différence de deux tables: retire les lignes qui sont dans 'other'"""
        if self.data.dtype.names != other.data.dtype.names:
            UTMESS("A", "TABLE0_30")
            return Table()
        else:
            void_data_1 = self.data.filled().view((np.void, self.data.dtype.itemsize))
            void_data_2 = other.data.filled().view((np.void, other.data.dtype.itemsize))
            mask = ~np.isin(void_data_1, void_data_2)
            new_data = self.data[mask]
            return Table.from_masked_array(new_data, titr=self._title)

    def dict_CREA_TABLE(self):
        """Renvoie le dictionnaire des mots-clés à fournir à la commande CREA_TABLE
        pour produire une table_sdaster.
        """
        titr = cut_long_lines(self._title, 80)
        # il y a eu limite à 50 titres dans le fortran autant le limiter
        # maintenant
        dico = {"TITRE": ["%-80s" % lig for lig in titr.split("\n")][:1], "LISTE": []}
        # remplissage de chaque occurrence (pour chaque paramètre) du mot-clé
        # facteur LISTE
        for name in self.data.dtype.names:
            # nom du paramètre et type si K*
            d = {"PARA": name}
            typ = self._get_aster_type(name)
            if typ[0] == "K":
                mc = "LISTE_K"
                if typ not in ("K8", "K16", "K24"):
                    UTMESS("A", "TABLE0_32", valk=(name, K_DEF))
                    typ = K_DEF
                d["TYPE_K"] = typ
            elif typ == "I":
                mc = "LISTE_I"
            elif typ == "R":
                mc = "LISTE_R"
            else:
                UTMESS("F", "TABLE0_31", valk=name)
            # valeurs sans trou / avec trou
            values = self.data[name]
            if typ == "R":
                try:
                    check_nan(values)
                    check_inf(values)
                except ValueError as err:
                    UTMESS("F", "TABLE0_33", valk=(name, str(err)))
            if not np.any(values.mask):
                d[mc] = values.tolist()
            else:
                values_compressed = values.compressed()
                d["NUME_LIGN"] = [i+1 for i in range(len(values.mask)) if not values.mask[i]]
                d[mc] = values_compressed.tolist()
            if len(d[mc]) == 0:
                UTMESS("I", "TABLE0_34", valk=name)
            else:
                dico["LISTE"].append(d)
        if len(dico["LISTE"]) == 0:
            UTMESS("A", "TABLE0_35")
        return dico

    def values(self) -> dict[str, list[Any]]:
        """Renvoie la table sous la forme d'un dictionnaire de listes dont les
        clés sont les paramètres.
        """
        return {name: self.data[name].tolist() for name in self.data.dtype.names}

    def Array(self, Para: str, Champ: str) -> np.ndarray:
        new_array = np.zeros((len(self.data), 2), dtype=float)
        new_array[:, 0] = self.data[Para].astype(float).filled(np.nan)
        new_array[:, 1] = self.data[Champ].astype(float).filled(np.nan)
        return new_array

    def Croise(self) -> "Table":
        """Retourne un tableau croisé P3(P1,P2) à partir d'une table ayant
        trois paramètres (P1, P2, P3).
        """
        if len(self.data.dtype.names) != 3:
            UTMESS("A", "TABLE0_36")
            return Table()
        py, px, pz = self.data.dtype.names
        x = self.data[px].data
        x_mask = self.data[px].mask
        y = self.data[py].data
        y_mask = self.data[py].mask
        z = self.data[pz].data
        z_mask = self.data[pz].mask
        # names of the rows and columns
        rows = np.unique(x[~x_mask])
        columns = np.unique(y[~y_mask])
        # Group by (x, y) and pick the last one if duplicate
        nb_values = len(x)
        rev = np.arange(nb_values)[::-1]
        idx = np.lexsort((rev, y, x))
        x_ord, y_ord, z_ord = x[idx], y[idx], z[idx]
        z_mask_ord = z_mask[idx]
        # Remove duplicates
        pairs = np.empty(x_ord.shape[0], dtype=[("x", x_ord.dtype), ("y", y_ord.dtype)])
        pairs["x"] = x_ord
        pairs["y"] = y_ord
        _, first = np.unique(pairs, return_index=True)
        x_last, y_last, z_last = x_ord[first], y_ord[first], z_ord[first]
        z_mask_last = z_mask_ord[first]
        # Restructured into a 2d table
        idx_row = rows.searchsorted(x_last)
        idx_column = columns.searchsorted(y_last)
        nb_rows = len(rows)
        nb_columns = len(columns)
        z_data = np.zeros((nb_rows, nb_columns), dtype=z.dtype)
        z_data[idx_row, idx_column] = z_last
        z_mask_data = np.full((nb_rows, nb_columns), True, dtype=np.bool_)
        z_mask_data[idx_row, idx_column] = z_mask_last
        # Build the new data of the new table
        new_names = [f"{px}/{py}", *columns.astype(str)]
        arrays = [rows, *(z_data[:, i] for i in range(nb_columns))]
        masks = [
            np.zeros(nb_rows, dtype=bool),
            *(z_mask_data[:, i] for i in range(nb_columns)),
        ]
        new_data = np.ma.MaskedArray(
            data=np.rec.fromarrays(arrays, names=new_names),
            mask=np.rec.fromarrays(masks, names=new_names),
        )
        new_title = (
            self._title + ("\n" if self._title else "") + f"{pz} FONCTION DE {px} ET {py}"
        )
        return Table.from_masked_array(new_data, titr=new_title)

    def Renomme(self, pold, pnew):
        """Rename the old column name into the new column name.

        Raise a KeyError if the old one is missing or the new one exist already.
        """
        if pold not in self.data.dtype.names:
            raise KeyError("Paramètre %s inexistant dans cette table" % pold)
        elif pnew in self.data.dtype.names:
            raise KeyError("Le paramètre %s existe déjà dans la table" % pnew)
        else:
            self.data = rfn.rename_fields(self.data, {pold: pnew})

    def _create_recarray(
        self,
        para: Sequence[str],
        typ: Sequence[str],
        rows: Optional[Sequence[Mapping[str, Any]]] = None,
    ) -> np.ma.MaskedArray:
        nb_rows = len(rows) if rows else 0
        dtype = []
        for i, (column_name, column_type) in enumerate(zip(para, typ)):
            if column_type is None:
                self._none_types.add(i)
            dtype.append((column_name, self._get_numpy_type(column_type, column_name)))
        data = np.recarray((nb_rows,), dtype=dtype, titles=para)
        mask = np.zeros(data.shape[0], dtype=[(n, bool) for n in para])
        if rows is not None:
            for i, row in enumerate(rows):
                for col_name in para:
                    value = row[col_name] if col_name in row else None
                    if value is not None:
                        data[col_name][i] = row[col_name]
                    else:
                        mask[col_name][i] = True
        return np.ma.MaskedArray(data, mask=mask)

    def _get_numpy_type(self, column_type: str, column_name: Optional[str] = None) -> np.dtype:
        """Get the python type associated to the given Aster type.

        Args:
            column_name (str): Name of the column.
            column_type (str): Aster type of the column

        Returns:
            dtype: Numpy type
        """
        if column_type is None:
            return np.int
        elif column_type[0] == "K":
            if column_type in ("K8", "K16", "K24"):
                return np.dtype(column_type.replace("K", "U"))
            else:
                if column_name is not None:
                    UTMESS("A", "TABLE0_32", valk=(column_name, K_DEF))
                msg = f"The type '{column_type}' is not supported at the moment"
                raise ValueError(msg)
        elif column_type == "I":
            return np.int
        elif column_type == "R":
            return np.float
        else:
            if column_name is not None:
                UTMESS("F", "TABLE0_31", valk=column_name)
            else:
                msg = f"The type '{column_type}' is not supported at the moment"
                raise ValueError(msg)

    def _get_aster_type(self, column_name: str) -> str:
        dt = self.data[column_name].dtype
        if np.issubdtype(dt, np.floating):
            return "R"
        elif np.issubdtype(dt, np.int):
            return "I"
        elif np.issubdtype(dt, np.complex):
            return "C"
        elif dt.kind == "U":
            length = dt.itemsize // np.dtype("U1").itemsize
            return f"K{length}"
        else:
            return f"K24"

    def _to_masked_array(self, values: list[Any], name: str, aster_type: "str"):
        """
        Convert a sequence `values` (where None means “missing”) into a
        NumPy MaskedArray of the given scalar dtype, with all the None’s masked.
        Runs in pure C‐loops under the hood (via fromiter), so it’s very fast.
        """
        dtype = self._get_numpy_type(aster_type)
        n = len(values)
        mask = np.fromiter((v is None for v in values), dtype=bool, count=n)
        default = np.zeros((), dtype=dtype).item()
        data = np.fromiter(
            (v if v is not None else default for v in values), dtype=dtype, count=n
        )
        return np.ma.MaskedArray(data, mask=mask, dtype=[(name, dtype)])

    def __iter__(self):
        """Itère sur les lignes de la Table"""
        names = self.data.dtype.names
        infos = [
            (name, self.data.data[name], np.ma.getmaskarray(self.data[name])) for name in names
        ]
        n = len(self.data)
        for i in range(n):
            yield MappingProxyType({name: values[i] for name, values, m_values in infos if not m_values[i]})

    def __len__(self):
        """Retourne le nombre de ligne dans la Table/Colonne."""
        return len(self.data)

    def __getitem__(self, args):
        """Extrait la sous table composée des colonnes dont les paramètres sont dans args"""
        if not is_sequence(args):
            args = [args]
        else:
            args = list(args)
        keep_params = [name for name in self.data.dtype.names if name in args]
        if len(keep_params) > 0:
            mask = np.logical_and.reduce(
                [self.data[keep_params].mask[name] for name in keep_params]
            )
            return Table.from_masked_array(
                rfn.repack_fields(self.data[keep_params][~mask]), nom=self._name
            )
        return Table(titr=self._title, nom=self._name)

    def __setitem__(self, k_para, k_value):
        """Ajoute une colonne k_para dont les valeurs sont dans k_value"""
        if k_para in self.data.dtype.names:
            UTMESS("F", "TABLE0_24", valk=k_para)
        # check type coherence of the values
        last_type = typaster(k_value[0])
        for v in it.islice(k_value, 1, None):
            last_type = typaster(v, last_type)
        # transform the given values to a masked array
        values = self._to_masked_array(k_value, k_para, last_type)
        # combine the types
        combined_dtypes = list(self.data.dtype.descr) + list(values.dtype.descr)
        # find the new length of the table
        nb_values = len(values)
        nb_data = len(self.data)
        new_len = nb_data
        if nb_values > nb_data:
            new_len = nb_values
        # create an empty masked record array
        new_data = np.ma.masked_all(new_len, dtype=np.dtype(combined_dtypes))
        # copy old data
        new_data[list(self.data.dtype.names)][:nb_data] = self.data
        # copy new column
        new_data[k_para][:nb_values] = values
        self.data = new_data

    def __delitem__(self, args):
        """Supprime les colonnes correspondantes aux éléments de args"""
        if not is_sequence(args):
            args = [args]
        args = set(args)
        if len(args) > 0:
            unknown_params = args - set(self.data.dtype.names)
            if len(unknown_params) > 0:
                UTMESS("F", "TABLE0_27", valk=str(unknown_params))
            keep = [name for name in self.data.dtype.names if name not in args]
            self.data = rfn.repack_fields(self.data[keep])

    def __and__(self, other: "Table") -> "Table":
        """Intersection with the given table.

        Implement the operator &
        """
        if other.data.dtype.names != self.data.dtype.names:
            UTMESS("A", "TABLE0_30")
            return Table()
        else:
            void_data_1 = self.data.filled().view((np.void, self.data.dtype.itemsize))
            void_data_2 = other.data.filled().view((np.void, other.data.dtype.itemsize))
            mask = np.isin(void_data_1, void_data_2)
            new_data = self.data[mask]
            return Table.from_masked_array(new_data, titr=self._title)

    def __or__(self, other):
        """Union de deux tables (opérateur |)"""
        if self.data.dtype.names != other.data.dtype.names:
            UTMESS("A", "TABLE0_30")
            return Table()
        else:
            void_data_1 = self.data.filled().view((np.void, self.data.dtype.itemsize))
            void_data_2 = other.data.filled().view((np.void, other.data.dtype.itemsize))
            mask = ~np.isin(void_data_1, void_data_2)
            new_data = np.ma.concatenate((self.data, other.data[mask]), axis=0)
            return Table.from_masked_array(new_data, titr=self._title)

    def __getattr__(self, column: str) -> "Colonne":
        """Construit un objet intermediaire (couple table, colonne)"""
        if column not in self.data.dtype.names:
            return Colonne(self, "", None)
        else:
            return Colonne(self, column, self._get_aster_type(column))

    def __getstate__(self):
        """Return internal state.

        Returns:
            tuple: Internal state.
        """
        return (self.data, self._title, self._name, self._none_types)

    def __setstate__(self, state):
        """Restore internal state.

        Arguments:
            state (tuple): Internal state.
        """
        if isinstance(state, tuple) and len(state) == 4:
            self.data, self._title, self._name, self._none_types = state


class Colonne:
    """Classe intermédiaire pour mémoriser un couple (table, nom de colonne)
    et exprimer les critères d'extraction sous une forme naturelle en python
    en surchargeant les operateurs <, >, != et =.
    Alors on peut écrire la requete simple :
      soustable=t.a<10
    Ainsi que des requetes plus complexes :
      soustable=t.a<10 and t.b <4
    ou
      soustable=t.a<10 or t.b <4
    Les "alias" EQ, NE, LE, LT, GE, GT permettent à la macro IMPR_TABLE
    d'utiliser directement le mot-clé utilisateur CRIT_COMP défini dans le
    catalogue : getattr(Table,CRIT_COMP).
    """

    def __init__(self, table: Table, column: str, typ=None):
        self._table = table
        self._column = column
        self._type = typ
        self._title = ""

    @property
    def para(self) -> str:
        return self._column

    @property
    def titr(self) -> str:
        return self._title

    @property
    def data(self) -> np.ma.MaskedArray:
        return self._table.data

    def __le__(self, VALE: Any) -> bool:
        if is_sequence(VALE):
            crit = max(VALE)
        else:
            crit = VALE
        column = self._table.data[self._column]
        mask = (column.data <= crit) & ~column.mask
        return self._get_filtered_table(mask)

    def __lt__(self, VALE: Any) -> bool:
        if is_sequence(VALE):
            crit = max(VALE)
        else:
            crit = VALE
        column = self._table.data[self._column]
        mask = (column.data < crit) & ~column.mask
        return self._get_filtered_table(mask)

    def __ge__(self, VALE: Any) -> bool:
        if is_sequence(VALE):
            crit = min(VALE)
        else:
            crit = VALE
        column = self._table.data[self._column]
        mask = (column.data >= crit) & ~column.mask
        return self._get_filtered_table(mask)

    def __gt__(self, VALE: Any) -> bool:
        if is_sequence(VALE):
            crit = min(VALE)
        else:
            crit = VALE
        column = self._table.data[self._column]
        mask = (column.data > crit) & ~column.mask
        return self._get_filtered_table(mask)

    def __eq__(self, VALE: Any, CRITERE="RELATIF", PRECISION=0.0):
        mask_eq = self._get_equal_mask(VALE, CRITERE, PRECISION)
        return self._get_filtered_table(mask_eq)

    def __ne__(self, VALE, CRITERE="RELATIF", PRECISION=0.0):
        mask_eq = ~self._get_equal_mask(VALE, CRITERE, PRECISION)
        return self._get_filtered_table(mask_eq)

    def MAXI(self) -> Table:
        column = self._table.data[self._column]
        mask_eq = column == np.max(column)
        return self._get_filtered_table(mask_eq)

    def MINI(self) -> Table:
        column = self._table.data[self._column]
        mask_eq = column == np.min(column)
        return self._get_filtered_table(mask_eq)

    def MAXI_ABS(self) -> Table:
        column = self._table.data[self._column]
        max_value = np.max(np.abs(column))
        mask_eq = (column == max_value) | (column == -max_value)
        return self._get_filtered_table(mask_eq)

    def MINI_ABS(self) -> Table:
        column = self._table.data[self._column]
        min_value = np.min(np.abs(column))
        mask_eq = (column == min_value) | (column == -min_value)
        return self._get_filtered_table(mask_eq)

    def VIDE(self) -> Table:
        return self._get_filtered_table(np.ma.getmaskarray(self._table.data[self._column]))

    def NON_VIDE(self) -> Table:
        return self._get_filtered_table(~np.ma.getmaskarray(self._table.data[self._column]))

    def REGEXP(self, regexp: str) -> Table:
        column = self._table.data[self._column]
        if isinstance(regexp, str):
            pattern = re.compile(regexp)
            regex_search = np.vectorize(lambda x: pattern.search(x) is not None)
            mask = regex_search(column.data) & ~column.mask
        else:
            mask = np.zeros(len(column), dtype=bool)
        return self._get_filtered_table(mask)

    def values(self) -> list[Any]:
        """Get the values of the column."""
        data = self._table.data
        if self._column in data.dtype.names:
            return data[self._column].tolist()
        else:
            return [None] * len(data)

    def not_none_values(self):
        """Get the values that are not None."""
        return self._table.data[self._column].compressed().tolist()

    def _get_filtered_table(self, mask: np.ndarray) -> Table:
        new_data = self._table.data[mask]
        return Table.from_masked_array(new_data, titr=self._table.titr)

    def _get_equal_mask(self, VALE, CRITERE="RELATIF", PRECISION=0.0) -> np.ndarray:
        if not is_sequence(VALE):
            VALE = [VALE]
        column = self._table.data[self._column]
        if isinstance(VALE[0], str):
            stripped_vale = [v.strip() for v in VALE]
            column_str = column.astype(np.str)
            mask_eq = np.isin(column_str, stripped_vale)
        elif PRECISION == 0.0:
            mask_eq = np.isin(column, VALE)
        elif CRITERE == "ABSOLU":
            # can't use np.abs to keep the same behavior as the old version
            column_arr = column[:, None]
            mask_eq = np.any(
                ((column_arr - PRECISION) <= VALE) & ((column_arr + PRECISION) >= VALE), axis=1
            )
        else:
            sign = np.sign(column)
            sign[sign == 0] = -1.0
            signed_VALE = sign[:, None] * VALE
            signed_column = (sign * column)[:, None]
            mask_eq = np.any(
                ((1.0 - PRECISION) * signed_column <= signed_VALE)
                & ((1.0 + PRECISION) * signed_column >= signed_VALE),
                axis=1,
            )
        return mask_eq

    def __len__(self):
        return len(self._table.data)

    def __iter__(self):
        """Iterate over the elements of the column."""
        column = self._table.data[self._column]
        for mask, value in zip(column.mask, column.data):
            if mask:
                yield None
            else:
                yield value

    def __getitem__(self, i: Union[int, slice]):
        """Get the i-th value(s) of the column"""
        data = self._table.data[self._column]
        if isinstance(i, slice):
            return data[i].tolist()
        else:
            return None if data.mask[i] else data.data[i]

    def __getstate__(self):
        """Return internal state.

        Returns:
            tuple: Internal state.
        """
        return (self._table, self._column, self._type, self._title)

    def __setstate__(self, state):
        """Restore internal state.

        Arguments:
            state (tuple): Internal state.
        """
        if isinstance(state, tuple) and len(state) == 4:
            (self._table, self._column, self._type, self._title) = state

    # Aster operators equivalents
    LE = __le__
    LT = __lt__
    GE = __ge__
    GT = __gt__
    EQ = __eq__
    NE = __ne__


def merge(tab1: Table, tab2: Table, labels: Sequence[str]=(), restrict: bool=False, format_r: str=None):
    """Assemble les deux tables tb1 et tb2 selon une liste de labels communs.
    Si labels est vide:
     - les lignes de tb2 sont ajoutés à celles de tb1,
    sinon :
     - si on trouve les valeurs de tb2 sur les labels dans tb1 (et une seule fois),
       on surcharge tb1 avec les lignes de tb2 ;
     - sinon on ajoute la ligne de tb2 à la fin de tb1.
    Si format_r est fourni, on "arrondi" les réels selon ce format lors de la
    comparaison.
    """
    tab1_para = tab1.para
    tab2_para = tab2.para
    tab1_para_set = set(tab1_para)
    tab2_para_set = set(tab2_para)
    if not is_sequence(labels):
        labels = (labels,)
    labels_unique = []
    for label in labels:
        if label not in labels_unique:
            labels_unique.append(label)
    labels = labels_unique
    for key in labels:
        if key not in tab1_para_set:
            UTMESS("F", "TABLE0_27", valk=key)
        if key not in tab2_para_set:
            UTMESS("F", "TABLE0_27", valk=key)
    # ensemble des paramètres et des types
    title = "\n".join([tab1.titr, tab2.titr, "MERGE avec labels=%s" % repr(tuple(labels))])
    if len(labels) > 0:
        new_names = [name for name in tab2_para if name not in tab1_para_set]
        if len(new_names) == 0:
            new_data = tab1.data.copy()
        else:
            dtypes_to_add = [(name, tab2.data.dtype.fields[name][0]) for name in new_names]
            combined_dtypes = list(tab1.data.dtype.descr) + dtypes_to_add
            new_data = np.ma.masked_all(len(tab1.data), dtype=np.dtype(combined_dtypes))
            for name in tab1_para:
                new_data[name] = tab1.data[name]
        if format_r is None:
            data1 = tab1.data[labels]
            data2 = tab2.data[labels]
        else:
            data1 = _reformat_masked_recarray(tab1.data, format_r, labels)
            data2 = _reformat_masked_recarray(tab2.data, format_r, labels)
        common_dtype = _find_common_dtype(data1, data2)
        data1 = data1.astype(common_dtype)
        data2 = data2.astype(common_dtype)
        # TODO: check if the unique doesn't cause issues with the mask
        tab1_unique, tab1_index, tab1_counts = np.unique(data1, return_index=True, return_counts=True)
        tab1_index = tab1_index[tab1_counts == 1]
        tab1_unique = data1[tab1_index]
        _, tab2_index, tab2_counts = np.unique(data2, return_index=True, return_counts=True)
        tab2_index = tab2_index[tab2_counts == 1]
        tab2_unique = data2[tab2_index]
        (_, tab1_intersect_index, tab2_intersect_index,) = np.intersect1d(tab1_unique, tab2_unique, assume_unique=True, return_indices=True)
        replace_index_1 = tab1_index[tab1_intersect_index]
        replace_index_2 = tab2_index[tab2_intersect_index]
        for name in tab2_para:
            mask = ~(tab2.data.mask[name][replace_index_2])
            new_data[name][replace_index_1[mask]] = tab2.data.data[name][replace_index_2[mask]]
        if not restrict:
            matched = np.full(len(tab2), True, dtype=bool)
            matched[replace_index_2] = False
            new_data =rfn.stack_arrays(
                [new_data, tab2.data[matched]],
                usemask=True,
                autoconvert=True,
            )
    elif restrict:
        new_names = [name for name in tab2_para if name not in tab1_para_set]
        new_data = tab1.data.copy()
        if len(new_names) > 0:
            dtypes_to_add = [tab2.data.dtype.fields[name][0] for name in new_names]
            empty_columns = [np.ma.masked_all(len(tab1.data), dtype=dt) for dt in dtypes_to_add]
            new_data = rfn.append_fields(
                tab1.data,
                names=new_names,
                data=empty_columns,
                dtypes=dtypes_to_add,
                usemask=True,
            )
    else:
        new_data = rfn.stack_arrays(
            [tab1.data, tab2.data],
            usemask=True,
            asrecarray=False,
            autoconvert=True,
        )
    return Table.from_masked_array(new_data, titr=title)


def _find_common_dtype(data_1: np.ma.MaskedArray, data_2: np.ma.MaskedArray):
    descr = []
    dt1 = data_1.dtype
    dt2 = data_2.dtype
    for name in dt1.names:
        t1 = dt1.fields[name][0]
        t2 = dt2.fields[name][0]
        if t1.kind == 'U' and t2.kind == 'U':
            u_size = max(t1.itemsize // 4, t2.itemsize // 4)
            descr.append((name, f'U{u_size}'))
        else:
            descr.append((name, np.promote_types(t1, t2)))
    return np.dtype(descr)


def remove_twins(
    tab: Table, labels: Sequence[str] = (), format_r: Optional[str] = None
) -> None:
    """Remove lines if the values of `labels` have already been seen"""
    if format_r is None:
        _, idx = np.ma.unique(tab.data[list(labels)], return_index=True)
        tab.data = tab.data[idx]
    else:
        new_data = _reformat_masked_recarray(tab.data, format_r, labels)
        _, new_idx = np.unique(new_data, return_index=True)
        tab.data = tab.data[new_idx]


def _reformat_masked_recarray(
    data: np.ma.MaskedArray, format_r: str, labels: Sequence[str]
) -> np.ma.MaskedArray:
    # py_format_r = fmtF2PY(format_r)
    decimals = _extract_decimals(format_r)
    arrays = []
    for label in labels:
        column = data[label]
        if np.issubdtype(column.dtype, np.float):
            # new_column = np.char.mod(py_format_r, column)
            rounded = np.ma.round(column.astype(float), decimals=decimals)
            new_column = rounded.astype(column.dtype, copy=False)
            arrays.append(new_column)
        else:
            arrays.append(column)
    return mr.fromarrays(arrays, names=labels)


def _extract_decimals(fformat, default: int=5):
    """Extract number of decimals from a fortran format.

    Manage only the real fortran, e.g E12.5, 1PE13.6, D12.5...

    Returns the default value if it failed to find the number of decimals. 
    """
    mat_decimal = re.search(r"[eEdDfFgG]{1}[0-9]+\.([0-9]+)", fformat)
    if mat_decimal is not None:
        try:
            return int(mat_decimal.group(1))
        except (ValueError, TypeError):
            pass
    return default


def typaster(obj, prev=None, strict=False):
    """Retourne le type Aster ('R', 'I', K_DEF) correspondant à l'objet obj.
    Si prev est fourni, on vérifie que obj est du type prev.
    Si strict=False, on autorise que obj ne soit pas du type prev s'ils sont
    tous les deux numériques ; dans ce cas, on retourne le "type enveloppe" 'R'.
    """
    dtyp = {
        int: "I",
        float: "R",
        str: K_DEF,
        type(None): "I",
        np.str_: K_DEF,
        np.int: "I",
        np.int64: "I",
        np.int32: "I",
        np.int8: "I",
        np.float: "R",
        np.float64: "R",
        np.float32: "R",
        np.float16: "R",
    }
    if is_float(obj):
        obj = float(obj)
    if type(obj) in dtyp:
        typobj = dtyp[type(obj)]
        if prev is None or prev == typobj:
            return typobj
        elif prev[0] == typobj[0] == "K":
            if len(obj) <= int(prev[1:]):
                return prev
            else:
                raise TypeError(
                    "La longueur de la chaine %s est incompatible avec le type %s"
                    % (repr(obj), repr(prev))
                )
        elif strict:  # prev is not None et typobj != prev et strict
            raise TypeError("La valeur %s n'est pas de type %s" % (repr(obj), repr(prev)))
        elif prev in ("I", "R") and typobj in ("I", "R"):
            return "R"
        else:
            raise TypeError(
                "La valeur %s n'est pas compatible avec le type %s" % (repr(obj), repr(prev))
            )
    else:
        raise TypeError(
            "Une table ne peut contenir que des entiers, réels " "ou chaines de caractères."
        )


def check_nan(values: np.ndarray) -> None:
    """Raise ValueError exception if nan is found in values."""
    bool_val = np.isnan(values)
    if bool_val.any():
        indexes = np.argwhere(bool_val)[:, 0].astype(np.str)
        raise ValueError("NaN present at indexes %s" % ",".join(indexes))


def check_inf(values: np.ndarray) -> None:
    """Raise ValueError exception if inf is found in values."""
    bool_val = np.isinf(values)
    if bool_val.any():
        indexes = np.argwhere(bool_val)[:, 0].astype(np.str)
        raise ValueError("Inf present at indexes %s" % ",".join(indexes))


def FMT(dform, nform, typAster=None, larg=0, val=""):
    """Retourne un format d'impression Python à partir d'un type Aster ('R','I',
    'K8', 'K16'...). Si typAster is None, retourne dform[nform].
       larg : largeur minimale du format (val permet de ne pas ajouter des blancs
       si la chaine à afficher est plus longue que le format, on prend le partie
       de ne pas tronquer les chaines)
    """
    if typAster is None:
        fmt = dform[nform]
    elif typAster in ("I", "R"):
        if nform == "formK":
            # convertit %12.5E en %-12s
            fmt = re.sub(r"([0-9]+)[\.0-9]*[diueEfFgG]+", r"-\g<1>s", dform["form" + typAster])
        else:
            fmt = dform[nform]
    else:
        # typAster = Kn
        fmt = "%-" + typAster[1:] + "s"
    # on ajoute éventuellement des blancs pour atteindre la largeur demandée
    if larg != 0:
        fmt = " " * max(min(larg - len(val), larg - len(fmt % 0)), 0) + fmt
    return fmt
