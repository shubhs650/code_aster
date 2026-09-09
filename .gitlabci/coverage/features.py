# coding=utf-8

"""
Build the list of combinations to be tested.
"""

import io
import os
import os.path as osp
import re
import sys
import tempfile
from pathlib import Path
from subprocess import Popen

from .show_keywords import ExtractKeywords

BASE = Path(__file__).parent


class Features:
    """This class extracts the possible combinations from code_aster source
    and installation directories.
    """

    def __init__(self, installdir, srcdir, logger):
        self._inst = installdir
        self._src = srcdir
        self._log = logger
        self._version = ""
        self._keywords = []  # cmde/mcf/mcs#value
        self._models = {}  # element: [phen#mod]
        self._calculs = {}  # element: [options]
        self._options = {}  # phen#mod: [options]

    def get_version(self):
        """Returns code_aster version from execution output.

        Returns:
            str: code_aster version.
        """
        return self._version

    def get_all(self):
        """Return the list of all features that should be tested.

        Returns:
            set(str): List of combinations.
        """
        self._get_keywords()
        self._get_options()
        res = set()
        res.update(self._options)
        res.update(self._keywords)
        return res

    def get_model(self, element):
        """ "Return the 'phen#mod' that uses the element.

        Returns:
            str: 'phen#mod'
        """
        found = self._models.get(element)
        if not found:
            if not (
                element.startswith("DDEPL_R")
                or element.startswith("DTEMP_R")
                or element.startswith("DPRES_C")
            ):
                print(f"ERROR: Element not found: {element}")
            return []
        return self._models[element]

    def _get_calculs(self):
        """Execute MAJ_CATA command and store the available elementary
        calculations.
        """
        self._log("INFO: executing 'code_aster.MAJ_CATA'...")
        run_aster = osp.join(self._inst, "bin", "run_aster")

        maj_cata = BASE / "maj_cata.py"
        export = tempfile.NamedTemporaryFile().name
        with open(export, "w") as fexp:
            fexp.write(f"F comm {maj_cata} D 1")
        tmpf = tempfile.NamedTemporaryFile().name
        with open(tmpf, "w") as fout:
            cmde = Popen([run_aster, export], stdout=fout, universal_newlines=True)
            cmde.wait()
        with open(tmpf, "r") as fout:
            output = fout.readlines()
        os.remove(tmpf)
        os.remove(export)

        re_vers = re.compile("Version +(?P<vers>[0-9]+\.[0-9]+\.[0-9]+)")
        calc = {}
        for line in output:
            mat = re_vers.search(line)
            if mat:
                self._version = mat.group("vers")
                self._log(f"INFO: code_aster version: {self._version}")
            if "&&CALCUL" not in line:
                continue
            line = line.strip().replace(" ", "")
            spl = line.split("/")
            assert len(spl) == 3, line
            _, elt, option = spl
            calc.setdefault(elt, [])
            calc[elt].append(option)
        self._calculs = calc

    def _get_models(self):
        """Build the list of available finite elements."""
        self._log("INFO: extracting finite elements...")
        sys.path.append(osp.join(self._src, "catalo"))
        from cataelem.elem import CataElem

        cel = CataElem()
        cel.build()

        models = {}
        for phen in cel.getPhenomenons():
            for modname, mod in phen.modelisations.items():
                for _, elt in mod.elements:
                    models.setdefault(elt.name, [])
                    models[elt.name].append("#".join([phen.code, modname]))
        self._models = models

    def _get_options(self):
        """Build the list of options available for each phen#mod."""
        if not self._calculs:
            self._get_calculs()
        if not self._models:
            self._get_models()
        self._log("INFO: grouping elements by phen#mod...")
        opts = {}
        flat = set()
        for elt, phenos in self._models.items():
            for phen in phenos:
                opts.setdefault(phen, [])
                for option in self._calculs[elt]:
                    if option not in opts[phen]:
                        opts[phen].append(option)
        for phen, options in opts.items():
            for option in options:
                flat.add("te:" + phen + "/" + option)
        self._options = flat

    def _get_keywords(self):
        """Extract the list of available keywords combinations."""
        self._log("INFO: extracting keywords...")
        orig = sys.stderr
        stream = io.StringIO()
        sys.stderr = stream
        words, errors = [], []
        sys.path.append(osp.join(self._inst, "lib", "aster"))
        try:
            words = ExtractKeywords().get()
            errors = stream.getvalue().splitlines()
        finally:
            stream.close()
            sys.stderr = orig
        for msg in errors:
            if "PRE_SEISME_NONL" in msg:
                continue
            sys.stderr.write(msg)
        kwd = set()
        for line in words:
            spl = line.split()
            key = "kw:" + spl.pop(0)
            if spl:
                key += "#" + spl[0]
            kwd.add(key)
        self._keywords = kwd
