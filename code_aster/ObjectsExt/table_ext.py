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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`Table` --- Representation of tables
**********************************************
"""

import aster
from libaster import Table, TableContainer, TableOfFunctions

from typing import Any

import numpy as np
try:
    import pandas as pd
except ImportError:
    HAS_PANDAS = False
else:
    HAS_PANDAS = True


from ..Objects.table_py import Table as TablePy
from ..Utilities import injector

from ..Objects.Serialization import InternalStateBuilder


class TableStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *Table*."""

    def restore(self, table):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            mesh (*DataStructure*): The *DataStructure* object to be restored.
        """
        super().restore(table)
        table.build()


@injector(Table)
class ExtendedTable:
    cata_sdj = "SD.sd_table.sd_table"
    internalStateBuilder = TableStateBuilder

    def __getitem__(self, key):
        """Retourne la valeur d'une cellule de la table.
        Exemple : TAB['INST', 1] retourne la 1ère valeur de la colonne 'INST'."""
        try:
            para, numlign = key
        except (TypeError, ValueError):
            raise RuntimeError("Table.__getitem__ takes exactly 2 arguments.")
        if para not in self.getParameters() or numlign > self.get_nrow():
            raise KeyError

        column = self.get_column(para)
        return column[numlign - 1]

    def get_column(self, para):
        """Retourne la colonne para"""
        exists, columnI, columnR, columnC, columnK = self.getValues(para)
        typ = self.getColumnType(para)
        if typ == "I":
            column = columnI
        elif typ == "R":
            column = columnR
        elif typ == "C":
            column = columnC
        else:
            column = columnK
        return [x if y else None for x, y in zip(column, exists)]

    def TITRE(self):
        """Retourne le titre d'une table Aster
        (Utile pour récupérer le titre et uniquement le titre d'une table dont
        on souhaite manipuler la dérivée).
        """
        titr = self.getTitle()
        if titr is None:
            titr = ""
        return titr

    def get_nrow(self):
        """Renvoie le nombre de lignes"""
        return self.getNumberOfLines()

    def get_nom_para(self):
        """Produit une liste des noms des colonnes"""
        return self.getParameters()

    def to_dataframe(self, filled_columns: dict[str, Any])->"pd.DataFrame":
        """Extract the table into a pandas DataFrame."""
        if not HAS_PANDAS:
            msg = (
                "Cannot create a DataFrame because the 'pandas' library is not installed.\n"
                "Please install it first."
            )
            raise ValueError(msg)
        names = self.getParameters()
        if not names:
            return pd.DataFrame()
        columns, masks = self._get_arrays(names)
        # fill the arrays with the given values for missing data
        for name, column, mask in zip(names, columns, masks):
            if name not in filled_columns:
                msg = (
                    f"Missing fill value for column '{name}'.\n"
                    f"Expected keys: {names}\n"
                    f"Provided keys: {list(filled_columns.keys())}"
                )
                raise ValueError(msg)
            fill_value = filled_columns[name]
            column[mask] = fill_value
        return pd.DataFrame({name: column for name, column in zip(names, columns)})

    def to_masked_recarray(self) -> np.ma.MaskedArray:
        """Extract the table into a masked record numpy array."""
        l_para = self.getParameters()
        if not l_para:
            return np.ma.MaskedArray()
        columns, masks = self._get_arrays(l_para)
        return np.ma.MaskedArray(
            np.core.records.fromarrays(columns, names=l_para),
            mask = np.core.records.fromarrays(masks, names=l_para)
        )

    def EXTR_TABLE(self, para=None):
        """Produit un objet TablePy à partir du contenu d'une table Aster.
        On peut limiter aux paramètres listés dans 'para'.
        """
        l_para = self.getParameters()
        if not l_para:
            return TablePy(titr=self.TITRE(), nom=self.getName())
        columns, masks = self._get_arrays(l_para)
        data = np.ma.MaskedArray(
            np.core.records.fromarrays(columns, names=l_para),
            mask = np.core.records.fromarrays(masks, names=l_para)
        )
        return TablePy.from_masked_array(data, titr=self.TITRE(), nom=self.getName())

    def _get_arrays(self, names: list[str]) -> tuple[list[np.ndarray], list[np.ndarray]]:
        """Get the arrays and their mask for each columns."""
        columns = []
        masks = []
        for name in names:
            exists, columnI, columnR, columnC, columnK = self.getValues(name)
            typ = self.getColumnType(name)
            if typ == "I":
                column = np.array(columnI, dtype=np.int)
            elif typ == "R":
                column = np.array(columnR, dtype=np.float)
            elif typ == "C":
                column = np.array(columnC, dtype=np.complex)
            else:
                column = np.array(columnK, dtype=typ.replace("K", "U"))
            columns.append(column)
            masks.append(~np.array(exists, dtype=bool))
        return columns, masks


@injector(TableOfFunctions)
class ExtendedTableOfFunctions:
    cata_sdj = "SD.sd_table_fonction.sd_table_fonction"


@injector(TableContainer)
class ExtendedTableContainer:
    cata_sdj = "SD.sd_table_container.sd_table_container"
