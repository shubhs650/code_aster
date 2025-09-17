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

# person_in_charge: albert.alarcon at edf.fr

from ..CodeCommands import (
    CALC_CHAMP,
    CALC_TABLE,
    CREA_CHAMP,
    CREA_TABLE,
    DEFI_LIST_REEL,
    IMPR_RESU,
    IMPR_TABLE,
    MECA_STATIQUE,
    POST_RELEVE_T,
    STAT_NON_LINE,
)
from ..Messages import UTMESS
from ..Objects.table_py import Table


def combinaison_charge_ops(self, **args):
    """
    combinaison des calculs avec des chargements mécaniques et thermiques
    """

    args = _F(args)

    MODELE_MECA = args.get("MODELE_MECA")
    BLOC = args.get("BLOC")
    CARA_ELEM_MECA = args.get("CARA_ELEM_MECA")
    EXCIT_MECA = args.get("EXCIT_MECA")

    CHAM_MATER_MECA = args.get("CHAM_MATER_MECA")
    COMPORTEMENT = args.get("COMPORTEMENT")
    MODELE_THER = args.get("MODELE_THER")
    BLOC_THER = args.get("BLOC_THER")
    CARA_ELEM_THER = args.get("CARA_ELEM_THER")
    TABLE_COEF = args.get("TABLE_COEF")
    EXCIT_THER = args.get("EXCIT_THER")
    LIST_INST_THER = args.get("LIST_INST_THER")
    CHAM_RESU = args.get("CHAM_RESU")
    TABLE_RESU = args.get("TABLE_RESU")
    IMPRESSION = args.get("IMPRESSION")
    UNITE = args.get("UNITE")

    # si besoin de STAT_NON_LINE
    l_nonl = False
    l_ther = False
    l_membrane = False

    comport = []
    if COMPORTEMENT is not None:
        for ii in range(len(COMPORTEMENT)):
            comport.append(COMPORTEMENT[ii]["RELATION"])

    if "MULTIFIBRE" in comport or "CABLE" in comport:
        l_nonl = True

    if (MODELE_THER is not None) and (BLOC_THER is not None) and (CARA_ELEM_THER is not None):
        l_membrane = True

    # ------------------------------------------------
    # Calcul du nombre de cas simples et combinés
    # ------------------------------------------------
    nb_char = 0
    nb_char_meca = 0
    nb_char_ther = 0
    nb_combther1 = 0
    nb_combther2 = 0

    # lire la table des coefs
    tabcoef = TABLE_COEF.EXTR_TABLE().values()
    coefpara = TABLE_COEF.EXTR_TABLE().para

    # verifier les chargements MECA
    nb_char_meca = len(EXCIT_MECA)
    for ii in range(nb_char_meca):
        if EXCIT_MECA[ii]["NOM_CHAR"] not in coefpara:
            UTMESS("F", "COMBCH_1", valk=EXCIT_MECA[ii]["NOM_CHAR"])

    # verifier les chargements THER
    if EXCIT_THER is not None:
        l_ther = True
        nb_char_ther = len(EXCIT_THER)
        for ii in range(nb_char_ther):
            if EXCIT_THER[ii]["NOM_CHAR"] not in coefpara:
                UTMESS("F", "COMBCH_1", valk=EXCIT_THER[ii]["NOM_CHAR"])
    else:
        nb_char_ther = 0

    # nombre des lignes dans la table =  nb de combinaison
    nb_comb = len(tabcoef[coefpara[0]])

    list_combther1 = []
    nom_combther1 = []
    list_combther2 = []
    nom_combther2 = []
    if nb_char_ther > 0:
        for ii in range(nb_comb):
            flag = 0
            nom_ther = []
            for jj in range(nb_char_ther):
                if tabcoef[EXCIT_THER[jj]["NOM_CHAR"]][ii] != 0:
                    flag = flag + 1
                    nom_ther.append(EXCIT_THER[jj]["NOM_CHAR"])

            if flag == 2:
                nb_combther2 = nb_combther2 + 1
                list_combther2.append(ii)
                nom_combther2.append(nom_ther)
            elif flag == 1:
                nb_combther1 = nb_combther1 + 1
                list_combther1.append(ii)
                nom_combther1.append(nom_ther)
            elif flag > 3:
                UTMESS("F", "COMBCH_2", vali=ii + 1)

    combther1 = {}
    combther1["NUME"] = list_combther1
    combther1["NOM"] = nom_combther1
    combther2 = {}
    combther2["NUME"] = list_combther2
    combther2["NOM"] = nom_combther2

    nb_char = nb_char_meca + nb_char_ther * 2
    nb_tot = nb_char + nb_comb + nb_combther1 + nb_combther2 * 3

    # préparer pour la table comb finale
    list_nom = []
    list_nom.append(BLOC.userName)

    for ii in range(len(EXCIT_MECA)):
        list_nom.append(EXCIT_MECA[ii]["NOM_CHAR"])

    if nb_char_ther > 0:
        nomTHER = []
        for ii in range(len(EXCIT_THER)):
            nomTHER.append(EXCIT_THER[ii]["NOM_CHAR"])
            list_nom.append(EXCIT_THER[ii]["NOM_CHAR"] + "_MAX")
            list_nom.append(EXCIT_THER[ii]["NOM_CHAR"] + "_MIN")
            tabcoef[EXCIT_THER[ii]["NOM_CHAR"] + "_MAX"] = tabcoef[EXCIT_THER[ii]["NOM_CHAR"]]
            tabcoef[EXCIT_THER[ii]["NOM_CHAR"] + "_MIN"] = tabcoef[EXCIT_THER[ii]["NOM_CHAR"]]

    list_nb = [*range(nb_char + 1)]

    # mettre à jour la table des coef en ajoutant les nouveaux colonnes _MAX/MIN/SXYZ
    coefpara = list(tabcoef.keys())

    #
    # ---------------------------------------------------------------------
    # Impression des tables des chargements - combinaisons si demandé
    # ---------------------------------------------------------------------
    num_unite1 = -1
    num_unite2 = -1
    num_unite3 = -1

    for ii in range(len(TABLE_RESU)):
        if TABLE_RESU[ii]["OPTION"] == "COEF_COMB":
            nomco1 = TABLE_RESU[ii]["TABLE"]
            if TABLE_RESU[ii]["UNITE"] is not None:
                num_unite1 = TABLE_RESU[ii]["UNITE"]
        elif TABLE_RESU[ii]["OPTION"] == "CALC_COMB":
            nomco2 = TABLE_RESU[ii]["TABLE"]
            if TABLE_RESU[ii]["UNITE"] is not None:
                num_unite2 = TABLE_RESU[ii]["UNITE"]
        elif TABLE_RESU[ii]["OPTION"] == "EXTREMA":
            nomco3 = TABLE_RESU[ii]["TABLE"]
            if TABLE_RESU[ii]["UNITE"] is not None:
                num_unite3 = TABLE_RESU[ii]["UNITE"]

    numcomb = []
    for ii in range(nb_tot):
        numcomb.append("COMB_" + str(ii + 1))

    #
    # Création de la liste imbriquée col qui contient les coeff pour tous les cas
    # col[ii][:] => charge_ii    /  col[:][jj] => COMB_jj
    col = [None] * (nb_char + 1)
    for ii in range(nb_char + 1):
        col[ii] = [None] * nb_tot
    #
    # 1ère colonne pour BLOC
    for jj in range(nb_tot):
        col[0][jj] = 1.0

    # Chargements unitaires
    for jj in range(nb_char):
        for ii in range(1, nb_char + 1):
            if (ii - 1) == jj:
                col[ii][jj] = 1.0
            else:
                col[ii][jj] = 0.0

    # Chargements combinés
    flag1 = 0
    flag2 = 0
    jj = nb_char
    i_comb = 0
    while jj < nb_tot:
        for ii in range(1, nb_char + 1):
            col[ii][jj] = tabcoef[list_nom[ii]][i_comb]

        if nb_char_ther == 0:
            jj = jj + 1
        elif i_comb in combther1["NUME"]:
            # un chargement thermique dans la ligne
            for ii in range(1, nb_char + 1):
                col[ii][jj + 1] = tabcoef[list_nom[ii]][i_comb]

            ind0 = combther1["NUME"].index(i_comb)
            ind1 = list_nom.index(combther1["NOM"][ind0][0] + "_MAX")

            col[ind1][jj + 1] = 0.0
            col[ind1 + 1][jj] = 0.0

            jj = jj + 2
            flag1 = flag1 + 1

        elif i_comb in combther2["NUME"]:
            # deux chargements thermiques dans la ligne
            for ii in range(1, nb_char + 1):
                col[ii][jj + 1] = tabcoef[list_nom[ii]][i_comb]
                col[ii][jj + 2] = tabcoef[list_nom[ii]][i_comb]
                col[ii][jj + 3] = tabcoef[list_nom[ii]][i_comb]

            ind0 = combther2["NUME"].index(i_comb)
            ind1 = list_nom.index(combther2["NOM"][ind0][0] + "_MAX")
            ind2 = list_nom.index(combther2["NOM"][ind0][1] + "_MAX")
            # ther1_max
            col[ind1][jj + 2] = 0.0
            col[ind1][jj + 3] = 0.0
            # ther1_min
            col[ind1 + 1][jj] = 0.0
            col[ind1 + 1][jj + 1] = 0.0
            # ther2_max
            col[ind2][jj + 1] = 0.0
            col[ind2][jj + 3] = 0.0
            # ther2_min
            col[ind2 + 1][jj] = 0.0
            col[ind2 + 1][jj + 2] = 0.0

            jj = jj + 4
            flag2 = flag2 + 1
        else:
            # 0 chargement thermique dans la ligne
            jj = jj + 1

        i_comb = i_comb + 1

    # Création et impression de la table au format code aster
    coeff = {}
    coeff["LISTE"] = []
    coeff["LISTE"].append(_F(LISTE_K=numcomb, TYPE_K="K8", PARA="INTITULE"))
    for ii in range(nb_char + 1):
        coeff["LISTE"].append(_F(LISTE_R=col[ii], PARA=list_nom[ii]))

    tabcomb = CREA_TABLE(**coeff)
    self.register_result(tabcomb, nomco1)

    # -----------------------------------------------------------------------
    #                          CALCULs individuels
    # -----------------------------------------------------------------------
    calc = [None] * nb_char
    UTMESS("I", "COMBCH_4", valk=("CALCULS INDIVIDUELS - MECANIQUE",))

    ##### Chargement mécaniques linéaires ou non-linéaires
    for charge in range(nb_char_meca):
        if l_nonl:
            comport = {}
            comport["COMPORTEMENT"] = []
            for ii in range(len(COMPORTEMENT)):
                if COMPORTEMENT[ii]["RELATION"] in ["CABLE"]:
                    if COMPORTEMENT[ii]["GROUP_MA"] is not None:
                        comport["COMPORTEMENT"].append(
                            _F(
                                GROUP_MA=COMPORTEMENT[ii]["GROUP_MA"],
                                RELATION=COMPORTEMENT[ii]["RELATION"],
                                DEFORMATION="GROT_GDEP",
                            )
                        )
                    else:
                        comport["COMPORTEMENT"].append(
                            _F(
                                TOUT="OUI",
                                RELATION=COMPORTEMENT[ii]["RELATION"],
                                DEFORMATION="GROT_GDEP",
                            )
                        )
                else:
                    if COMPORTEMENT[ii]["GROUP_MA"] is not None:
                        comport["COMPORTEMENT"].append(
                            _F(
                                GROUP_MA=COMPORTEMENT[ii]["GROUP_MA"],
                                RELATION=COMPORTEMENT[ii]["RELATION"],
                            )
                        )
                    else:
                        comport["COMPORTEMENT"].append(
                            _F(TOUT="OUI", RELATION=COMPORTEMENT[ii]["RELATION"])
                        )

            instmeca = DEFI_LIST_REEL(VALE=(0.0, 1.0))
            calc[charge] = STAT_NON_LINE(
                MODELE=MODELE_MECA,
                CHAM_MATER=CHAM_MATER_MECA,
                CARA_ELEM=CARA_ELEM_MECA,
                EXCIT=(_F(CHARGE=BLOC), _F(CHARGE=EXCIT_MECA[charge]["CHAR_MECA"])),
                INCREMENT=_F(LIST_INST=instmeca, NUME_INST_FIN=1),
                CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=100),
                **comport,
            )
        else:
            calc[charge] = MECA_STATIQUE(
                MODELE=MODELE_MECA,
                CHAM_MATER=CHAM_MATER_MECA,
                CARA_ELEM=CARA_ELEM_MECA,
                EXCIT=(_F(CHARGE=BLOC), _F(CHARGE=EXCIT_MECA[charge]["CHAR_MECA"])),
            )

    nb = charge
    ##### Chargement thermiques

    if l_ther:
        UTMESS("I", "COMBCH_4", valk=("CALCULS INDIVIDUELS - THERMIQUE",))
        UTMESS("A", "COMBCH_3", valk=("MODELE_THER", "CARA_ELEM_THER", "BLOC_THER"))

    if l_membrane:
        model = MODELE_THER
        carael = CARA_ELEM_THER
        blocfin = BLOC_THER
    else:
        model = MODELE_MECA
        carael = CARA_ELEM_MECA
        blocfin = BLOC

    for charge in range(nb_char_ther):
        if l_nonl:
            comport = {}
            comport["COMPORTEMENT"] = []
            for ii in range(len(COMPORTEMENT)):
                if COMPORTEMENT[ii]["RELATION"] == "CABLE":
                    if COMPORTEMENT[ii]["GROUP_MA"] is not None:
                        comport["COMPORTEMENT"].append(
                            _F(
                                GROUP_MA=COMPORTEMENT[ii]["GROUP_MA"],
                                RELATION="CABLE",
                                DEFORMATION="GROT_GDEP",
                            )
                        )
                    else:
                        comport["COMPORTEMENT"].append(
                            _F(TOUT="OUI", RELATION="CABLE", DEFORMATION="GROT_GDEP")
                        )
                else:
                    if COMPORTEMENT[ii]["GROUP_MA"] is not None:
                        comport["COMPORTEMENT"].append(
                            _F(
                                GROUP_MA=COMPORTEMENT[ii]["GROUP_MA"],
                                RELATION=COMPORTEMENT[ii]["RELATION"],
                            )
                        )
                    else:
                        comport["COMPORTEMENT"].append(
                            _F(TOUT="OUI", RELATION=COMPORTEMENT[ii]["RELATION"])
                        )

            calc[nb + charge + 1] = STAT_NON_LINE(
                MODELE=model,
                CHAM_MATER=EXCIT_THER[charge]["CHAM_MATER_THER"],
                CARA_ELEM=carael,
                EXCIT=(_F(CHARGE=blocfin),),
                INCREMENT=_F(LIST_INST=LIST_INST_THER),
                CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=100),
                **comport,
            )
        else:
            calc[nb + charge + 1] = MECA_STATIQUE(
                MODELE=model,
                CHAM_MATER=EXCIT_THER[charge]["CHAM_MATER_THER"],
                CARA_ELEM=carael,
                LIST_INST=LIST_INST_THER,
                EXCIT=(_F(CHARGE=blocfin),),
            )

    UTMESS("I", "COMBCH_4", valk=("  FIN DE CALCULS INDIVIDUELS   ",))

    # --------------------------------------------------------------------------------------
    #                         Calculs des efforts et des contraintes pour chaque cas
    #                              les déplacements sont calculés par défaut
    # --------------------------------------------------------------------------------------

    for ii in range(nb_char - nb_char_ther):
        calc[ii] = CALC_CHAMP(
            reuse=calc[ii],
            RESULTAT=calc[ii],
            CONTRAINTE=("EFGE_ELNO", "EFGE_NOEU", "SIEF_ELNO"),
            FORCE=("FORC_NODA",),
        )

    type_champ = {
        "DEPL": "NOEU_DEPL_R",
        "EFGE_NOEU": "NOEU_SIEF_R",
        "SIEF_ELNO": "ELNO_SIEF_R",
        "FORC_NODA": "NOEU_DEPL_R",
    }
    nom_champ = {"DEPL": None, "EFGE_NOEU": None, "SIEF_ELNO": None, "FORC_NODA": None}
    nom_table = {"DEPL": None, "EFGE_NOEU": None, "SIEF_ELNO": None, "FORC_NODA": None}

    list_cham = []
    list_cmp = []
    for ii in range(len(CHAM_RESU)):
        list_cham.append(CHAM_RESU[ii]["NOM_CHAM"])
        if CHAM_RESU[ii]["NOM_CMP"] is not None:
            list_cmp.append(CHAM_RESU[ii]["NOM_CMP"])
        else:
            list_cmp.append("TOUT")

    list_sortie = []
    for cc in range(len(TABLE_RESU)):
        if TABLE_RESU[cc]["OPTION"] == "EXTREMA":
            if TABLE_RESU[cc]["CRIT_COMP"] is not None:
                if "TOUT" in TABLE_RESU[cc]["CRIT_COMP"]:
                    list_sortie = ["MAXI", "MINI", "MAXI_ABS", "MINI_ABS"]
                else:
                    list_sortie = TABLE_RESU[cc]["CRIT_COMP"]
            else:
                list_sortie = ["MAXI", "MINI", "MAXI_ABS", "MINI_ABS"]

    dict_list_sortie = {
        "MAXI": "MAX",
        "MINI": "MIN",
        "MAXI_ABS": "MAXI_ABS",
        "MINI_ABS": "MINI_ABS",
    }
    imprresu = {}
    imprresu["RESU"] = []
    comblist = []
    extrlist = []
    # for champ in list_cham :
    for kk in range(len(list_cham)):
        champ = list_cham[kk]
        nom_champ[champ] = [None] * nb_tot

        # Pour les chargements mécaniques
        for ii in range(nb_char_meca):
            # Création des champs de résultats  #
            nom_champ[champ][ii] = CREA_CHAMP(
                OPERATION="EXTR",
                TYPE_CHAM=type_champ[champ],
                RESULTAT=calc[ii],
                NUME_ORDRE=1,
                NOM_CHAM=champ,
            )
            medname = "COMB%03d_" % (ii + 1)
            medname = medname + champ
            imprresu["RESU"].append(_F(CHAM_GD=nom_champ[champ][ii], NOM_CHAM_MED=medname))

        # Pour les chargements thermiques
        for ii in range(nb_char_ther):
            # Création des champs de résultats
            nom_champ[champ][nb_char_meca + 2 * ii] = CREA_CHAMP(
                OPERATION="EXTR",
                TYPE_CHAM=type_champ[champ],
                RESULTAT=calc[nb_char_meca + ii],
                TOUT_ORDRE="OUI",
                TYPE_MAXI="MAXI",
                NOM_CHAM=champ,
            )
            nom_champ[champ][nb_char_meca + 2 * ii + 1] = CREA_CHAMP(
                OPERATION="EXTR",
                TYPE_CHAM=type_champ[champ],
                RESULTAT=calc[nb_char_meca + ii],
                TOUT_ORDRE="OUI",
                TYPE_MAXI="MINI",
                NOM_CHAM=champ,
            )
            medname = "COMB%03d_" % (nb_char_meca + 2 * ii + 1)
            medname = medname + champ
            imprresu["RESU"].append(_F(CHAM_GD=nom_champ[champ][ii], NOM_CHAM_MED=medname))
            medname = "COMB%03d_" % (nb_char_meca + 2 * ii + 1 + 1)
            medname = medname + champ
            imprresu["RESU"].append(_F(CHAM_GD=nom_champ[champ][ii], NOM_CHAM_MED=medname))

        for ii in range(nb_char, nb_tot):
            combinaison = {}
            combinaison["ASSE"] = []
            # Création des champs de résultats par combinaison linéaire des chargements unitaires
            for jj in range(nb_char):
                combinaison["ASSE"].append(
                    _F(
                        CUMUL="OUI",
                        TOUT="OUI",
                        CHAM_GD=nom_champ[champ][jj],
                        COEF_R=col[jj + 1][ii],
                    )
                )

            nom_champ[champ][ii] = CREA_CHAMP(
                MODELE=MODELE_MECA, OPERATION="ASSE", TYPE_CHAM=type_champ[champ], **combinaison
            )

            # preparer Impression du fichier résultat RMED
            medname = "COMB%03d_" % (ii + 1)
            medname = medname + champ
            imprresu["RESU"].append(_F(CHAM_GD=nom_champ[champ][ii], NOM_CHAM_MED=medname))

        # ---------------------------------------------------------------------------------------------------
        #                         Création et impression de la table de résultats pour tous les chargements
        #                         Unité, liste de champs, liste de composantes définies
        #                         dans le dictionnaire table_resu
        # ---------------------------------------------------------------------------------------------------

        nb_comp = len(list_cmp[kk])
        nom_table[champ] = [None] * nb_comp

        for comp in range(nb_comp):
            actcmp = {}
            actcmp["ACTION"] = []
            for ii in range(nb_tot):
                if "TOUT" in list_cmp[kk]:
                    actcmp["ACTION"].append(
                        _F(CHAM_GD=nom_champ[champ][ii], INTITULE=numcomb[ii], OPERATION="EXTREMA")
                    )
                else:
                    actcmp["ACTION"].append(
                        _F(
                            CHAM_GD=nom_champ[champ][ii],
                            INTITULE=numcomb[ii],
                            NOM_CMP=list_cmp[kk][comp],
                            OPERATION="EXTREMA",
                        )
                    )

            nom_table[champ][comp] = POST_RELEVE_T(**actcmp)
            nom_table[champ][comp] = CALC_TABLE(
                reuse=nom_table[champ][comp],
                TABLE=nom_table[champ][comp],
                ACTION=_F(OPERATION="SUPPRIME", NOM_PARA="CHAM_GD"),
            )

            comblist = comblist + list(nom_table[champ][comp].EXTR_TABLE().rows)

            # ---------------------------------------------------------------------------------------------------
            #                         Création et impression de la table des résultats enveloppe
            #                         Unité, liste de champs, liste de composantes et liste d'extremas définies
            #                         dans le dictionnaire table_enve
            # ---------------------------------------------------------------------------------------------------
            tab = [[None] * len(list_sortie) for i in range(nb_comp)]
            res = [None] * nb_comp
            for extr in range(len(list_sortie)):
                tab[comp][extr] = CALC_TABLE(
                    TABLE=nom_table[champ][comp],
                    ACTION=(
                        _F(
                            OPERATION="FILTRE",
                            NOM_PARA="EXTREMA",
                            VALE_K=dict_list_sortie[list_sortie[extr]],
                        ),
                        _F(OPERATION="FILTRE", NOM_PARA="VALE", CRIT_COMP=list_sortie[extr]),
                        _F(
                            OPERATION="AJOUT_COLONNE",
                            NOM_PARA="TYPE_VALE",
                            VALE_COLONNE=list_sortie[extr],
                        ),
                        _F(OPERATION="FILTRE", NOM_PARA="TYPE_VALE", CRIT_COMP="NON_VIDE"),
                    ),
                )

                if extr == 0:
                    res[comp] = tab[comp][extr]
                else:
                    res[comp] = CALC_TABLE(
                        reuse=res[comp],
                        TABLE=res[comp],
                        ACTION=_F(OPERATION="COMB", TABLE=tab[comp][extr], RESTREINT="NON"),
                    )

            res[comp] = CALC_TABLE(
                reuse=res[comp],
                TABLE=res[comp],
                ACTION=_F(OPERATION="SUPPRIME", NOM_PARA="CHAM_GD"),
            )
            extrlist = extrlist + list(res[comp].EXTR_TABLE().rows)

    combpara = nom_table[champ][comp].EXTR_TABLE().para
    combtype = nom_table[champ][comp].EXTR_TABLE().type
    extrpara = res[comp].EXTR_TABLE().para
    extrtype = res[comp].EXTR_TABLE().type

    # ------------------------------------------------
    #        Impression des tables et MED
    # ------------------------------------------------

    if num_unite1 > 0:
        IMPR_TABLE(TABLE=tabcomb, FORMAT="TABLEAU", UNITE=num_unite1)

    table2 = Table(comblist, combpara, combtype)
    motcles = table2.dict_CREA_TABLE()
    tabcalc = CREA_TABLE(TYPE_TABLE="TABLE", **motcles)
    self.register_result(tabcalc, nomco2)

    if num_unite2 > 0:
        IMPR_TABLE(
            TABLE=tabcalc,
            FORMAT="TABLEAU",
            UNITE=num_unite2,
            NOM_PARA=("INTITULE", "EXTREMA", "NOEUD", "CMP", "VALE"),
        )

    table3 = Table(extrlist, extrpara, extrtype)
    motcles = table3.dict_CREA_TABLE()
    tabextr = CREA_TABLE(TYPE_TABLE="TABLE", **motcles)
    self.register_result(tabextr, nomco3)

    if num_unite3 > 0:
        IMPR_TABLE(
            TABLE=tabextr,
            FORMAT="TABLEAU",
            UNITE=num_unite3,
            NOM_PARA=("INTITULE", "NOEUD", "CMP", "VALE", "TYPE_VALE"),
        )

    if IMPRESSION == "OUI":
        IMPR_RESU(UNITE=UNITE, **imprresu)

    return
