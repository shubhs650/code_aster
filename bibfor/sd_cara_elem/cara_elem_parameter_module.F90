! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
! as lint: disable=W1502
!       W1502: More than 19 continuation lines
!
module cara_elem_parameter_module
!
!
!   Pour l'opérateur AFFE_CARA_ELEM
!
! --------------------------------------------------------------------------------------------------
!
!&<
    implicit none
!
!   Pour l'utilisation des entiers codés :
!       - pas de zéro dans les entiers codé
!       - il ne faut pas dépasser 30
!       - Aucune vérification n'est faite
!
!   Utilisation
!       Icode : est l'entier codé
!       liste : [ i1, i2, i3 , ... ]   Pas nécessairement dans l'ordre
!       Codage
!           Icode = ACE_CodeEntier( liste )
!           Icode = ACE_CodeEntier( liste , Icode) Ajout des éléments de liste dans Icode
!
!       Décodage
!           ACE_IsInCode( Icode , vali )
!               si True    ==> vali fait partie de Icode
!               si False   ==> vali ne fait pas partie de Icode
!
!   Forme des différentes section
    integer(kind=8), parameter :: ACE_SECTION_GENERALE      = 0
    integer(kind=8), parameter :: ACE_SECTION_RECTANGLE     = 1
    integer(kind=8), parameter :: ACE_SECTION_CERCLE        = 2
!   Variation des sections le long de l'axe
    integer(kind=8), parameter :: ACE_SECTION_CONSTANTE     = 0
    integer(kind=8), parameter :: ACE_SECTION_AFFINE        = 1
    integer(kind=8), parameter :: ACE_SECTION_HOMOTHETIQUE  = 2
!
    integer(kind=8), parameter :: ACE_NOTHING        =  0
!
!   Nombre total de mots clefs facteur
    integer(kind=8), parameter :: ACE_NB_MCLEF = 17
    character(len=16), parameter :: ACE_MCLEF(ACE_NB_MCLEF) = [ &
        'POUTRE          ','COQUE           ','DISCRET         ','ORIENTATION     ', &
        'CABLE           ','BARRE           ','MASSIF          ', &
        'POUTRE_FLUI     ','RIGI_PARASOL    ','GRILLE          ','RIGI_MISS_3D    ', &
        'DISCRET_2D      ','MEMBRANE        ','MASS_AJOU       ','MULTIFIBRE      ', &
        'MASS_REP        ','RIGI_GRILLE     ']
!
    integer(kind=8), parameter :: ACE_POUTRE         =  1
    integer(kind=8), parameter :: ACE_COQUE          =  2
    integer(kind=8), parameter :: ACE_DISCRET        =  3
    integer(kind=8), parameter :: ACE_ORIENTATION    =  4
    integer(kind=8), parameter :: ACE_CABLE          =  5
    integer(kind=8), parameter :: ACE_BARRE          =  6
    integer(kind=8), parameter :: ACE_MASSIF         =  7
    integer(kind=8), parameter :: ACE_POUTRE_FLUI    =  8
    integer(kind=8), parameter :: ACE_RIGI_PARASOL   =  9
    integer(kind=8), parameter :: ACE_GRILLE         = 10
    integer(kind=8), parameter :: ACE_RIGI_MISS_3D   = 11
    integer(kind=8), parameter :: ACE_DISCRET_2D     = 12
    integer(kind=8), parameter :: ACE_MEMBRANE       = 13
    integer(kind=8), parameter :: ACE_MASS_AJOU      = 14
    integer(kind=8), parameter :: ACE_MULTIFIBRE     = 15
    integer(kind=8), parameter :: ACE_MASS_REP       = 16
    integer(kind=8), parameter :: ACE_RIGI_GRILLE    = 17
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ACE_NB_GRMA_MA = 3
    character(len=16), parameter :: ACE_GRMA_MA(ACE_NB_GRMA_MA) =[ &
        'GROUP_MA        ','GROUP_MA_POI1   ','GROUP_MA_SEG2   ' ]
    character(len=8),  parameter :: ACE_GRMA_TY(ACE_NB_GRMA_MA) =[ &
        'GROUP_MA',        'GROUP_MA',        'GROUP_MA' ]
    integer(kind=8), parameter :: ACE_GR_MAI  = 1
    integer(kind=8), parameter :: ACE_GR_PO1  = 2
    integer(kind=8), parameter :: ACE_GR_SE2  = 3
!
    integer(kind=8), parameter :: MCLEF_GRP_MA(ACE_NB_GRMA_MA*ACE_NB_MCLEF) = [ &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_GR_PO1,  ACE_GR_SE2 , &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_PO1, ACE_GR_SE2,  ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_GR_PO1,  ACE_NOTHING, &
        ACE_GR_MAI, ACE_NOTHING, ACE_NOTHING, &
        ACE_GR_MAI, ACE_GR_PO1,  ACE_NOTHING, &
        ACE_GR_MAI, ACE_GR_SE2,  ACE_NOTHING]
!
! --------------------------------------------------------------------------------------------------
!   Toutes les cartes créées par AFFE_CARA_ELEM
    integer(kind=8), parameter :: ACE_NB_CARTE     = 7
    integer(kind=8), parameter :: ACE_NB_CARTE_CMP = 3
    character(len=10), parameter :: ACE_CARTE(ACE_NB_CARTE*ACE_NB_CARTE_CMP) = [ &
        '.CARDINFO ', 'CINFDI_R  ', 'DISCRET   ', &
        '.CARDISCK ', 'CADISK_R  ', 'DISCRET   ', &
        '.CARDISCM ', 'CADISM_R  ', 'DISCRET   ', &
        '.CARDISCA ', 'CADISA_R  ', 'DISCRET   ', &
        '.CVENTCXF ', 'VENTCX_F  ', 'STRX_LINEI', &
        '.CARCABLE ', 'CACABL_R  ', 'CABLE     ', &
        '.CARGENBA ', 'CAGNBA_R  ', 'BARRE     ' ]
!
    integer(kind=8), parameter :: ACE_CAR_DINFO    =  1
    integer(kind=8), parameter :: ACE_CAR_DISCK    =  2
    integer(kind=8), parameter :: ACE_CAR_DISCM    =  3
    integer(kind=8), parameter :: ACE_CAR_DISCA    =  4
    integer(kind=8), parameter :: ACE_CAR_CVCXF    =  5
    integer(kind=8), parameter :: ACE_CAR_CABLE    =  6
    integer(kind=8), parameter :: ACE_CAR_BARRE    =  7
!
! Les cartes restantes
!           '.CARARCPO', 'CAARPO_R ', 'POUTRE    ',
!           '.CARGENPO', 'CAGNPO_R ', 'POUTRE    ',
!           '.CARGEOPO', 'CAGEPO_R ', 'POUTRE    ',
!           '.CARMASSI', 'CAMA_R   ', 'MASSIF    ',
!           '.CARCOQUE', 'CACOQU_R ', 'COQUE     ',
!           '.CARCOQUF', 'CACOQU_F ', 'COQUE     ',
!           '.CARPOUFL', 'CAPOUF_R ', 'POUTREFLUI',
!           '.CARORIEN', 'CAORIE_R ', 'ORIENT    ',
!
! --------------------------------------------------------------------------------------------------
!
!       ACE_NB_ELEMENT              : Nombre de type d'éléments différents
!       ACE_NM_ELEMENT              : Nom des types d'élément (POUTRE DISCRET COQUE ...)
!                                           el = ACE_NM_ELEMENT(ii)
!       ACE_NU_(el)                 : Numéro du type de l'élément (POUTRE DISCRET COQUE ...)
!       ACE_NB_(el)                 : Nombre de support dans le type
!       ACE_EL_(el)[ACE_NB_(el)]    : Liste des supports dans le type
!
!   Nombre de type d'éléments
    integer(kind=8), parameter :: ACE_NB_ELEMENT  =  10
    character(len=16),parameter :: ACE_NM_ELEMENT(ACE_NB_ELEMENT) =[ &
        'POUTRE          ', 'DISCRET         ', 'COQUE           ', 'CABLE           ', &
        'BARRE           ', 'GRILLE          ', 'MEMBRANE        ', 'MASSIF          ', &
        'MASSIF THM      ', 'MASSIF HHO      ']
!
    integer(kind=8), parameter :: ACE_NU_POUTRE    =  1
    integer(kind=8), parameter :: ACE_NU_DISCRET   =  2
    integer(kind=8), parameter :: ACE_NU_COQUE     =  3
    integer(kind=8), parameter :: ACE_NU_CABLE     =  4
    integer(kind=8), parameter :: ACE_NU_BARRE     =  5
    integer(kind=8), parameter :: ACE_NU_GRILLE    =  6
    integer(kind=8), parameter :: ACE_NU_MEMBRANE  =  7
    integer(kind=8), parameter :: ACE_NU_MASSIF    =  8
    integer(kind=8), parameter :: ACE_NU_THHMM     =  9
    integer(kind=8), parameter :: ACE_NU_HHO       = 10
!
    integer(kind=8), parameter :: ACE_NB_POUTRE    = 14
    character(len=16), parameter :: ACE_EL_POUTRE(ACE_NB_POUTRE) =[ &
        'MECA_POU_D_T    ', 'MECA_POU_D_E    ', 'MECA_POU_D_T_GD ', 'MEFS_POU_D_T    ', &
        'MECA_POU_D_TG   ', 'MECA_POHO_HEXA8 ', 'MECA_POHO_HEXA20', 'MET3SEG3        ', &
        'MET6SEG3        ', 'MET3SEG4        ', 'MECA_POU_D_EM   ', 'MECA_POU_D_TGM  ', &
        'MECA_POU_D_SQUE ', 'MECA_S8P2       ']
!
    integer(kind=8), parameter :: ACE_MECA_POU_D_T          =  1
    integer(kind=8), parameter :: ACE_MECA_POU_D_E          =  2
    integer(kind=8), parameter :: ACE_MECA_POU_D_T_GD       =  3
    integer(kind=8), parameter :: ACE_MEFS_POU_D_T          =  4
    integer(kind=8), parameter :: ACE_MECA_POU_D_TG         =  5
    integer(kind=8), parameter :: ACE_MECA_POHO_HEXA8       =  6
    integer(kind=8), parameter :: ACE_MECA_POHO_HEXA20      =  7
    integer(kind=8), parameter :: ACE_MET3SEG3              =  8
    integer(kind=8), parameter :: ACE_MET6SEG3              =  9
    integer(kind=8), parameter :: ACE_MET3SEG4              = 10
    integer(kind=8), parameter :: ACE_MECA_POU_D_EM         = 11
    integer(kind=8), parameter :: ACE_MECA_POU_D_TGM        = 12
    integer(kind=8), parameter :: ACE_MECA_POU_D_SQUE       = 13
    integer(kind=8), parameter :: ACE_MECA_S8P2             = 14
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_DISCRET   = 8
    character(len=16), parameter :: ACE_EL_DISCRET(ACE_NB_DISCRET) =[ &
        'MECA_DIS_T_N    ', 'MECA_DIS_T_L    ', 'MECA_DIS_TR_N   ', 'MECA_DIS_TR_L   ', &
        'MECA_2D_DIS_T_N ', 'MECA_2D_DIS_T_L ', 'MECA_2D_DIS_TR_N', 'MECA_2D_DIS_TR_L' ]
!
    integer(kind=8), parameter :: ACE_MECA_DIS_T_N      = 1
    integer(kind=8), parameter :: ACE_MECA_DIS_T_L      = 2
    integer(kind=8), parameter :: ACE_MECA_DIS_TR_N     = 3
    integer(kind=8), parameter :: ACE_MECA_DIS_TR_L     = 4
    integer(kind=8), parameter :: ACE_MECA_2D_DIS_T_N   = 5
    integer(kind=8), parameter :: ACE_MECA_2D_DIS_T_L   = 6
    integer(kind=8), parameter :: ACE_MECA_2D_DIS_TR_N  = 7
    integer(kind=8), parameter :: ACE_MECA_2D_DIS_TR_L  = 8
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_COQUE = 29
    character(len=16), parameter :: ACE_EL_COQUE(ACE_NB_COQUE) = [ &
                                    'THCOTR3         ', 'THCOTR6         ', &
                                    'THCOQU4         ', 'THCOQU8         ', &
                                    'THCOTR7         ', 'THCOQU9         ', &
                                    'MEDKTR3         ', 'MEDSTR3         ', &
                                    'MET3TR3         ', 'MEDKQU4         ', &
                                    'MEDSQU4         ', 'MEQ4QU4         ', &
                                    'MECXSE3         ', 'MEDKTG3         ', &
                                    'MEDKQG4         ', 'MEQ4GG4         ', &
                                    'MET3GG3         ', 'THCASE3         ', &
                                    'THCPSE3         ', 'MEC3QU9H        ', &
                                    'MEC3TR7H        ', 'MEBODKT         ', &
                                    'MEBODST         ', 'MEBOQ4G         ', &
                                    'MEBOCQ3         ', 'THCOSE3         ', &
                                    'THCOSE2         ', 'PLAQ_MITC_QUAD9 ', &
                                    'PLAQ_MITC_TRIA6 ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_CABLE     = 2
    character(len=16),parameter :: ACE_EL_CABLE(ACE_NB_CABLE) =[ &
        'MECABL2         ', 'MEPOULI         ' ]
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_BARRE     = 7
    character(len=16),parameter :: ACE_EL_BARRE(ACE_NB_BARRE) =[ &
        'MECA_BARRE      ', 'MECA_2D_BARRE   ', 'MECGSEG3        ', &
        'MEBA2DSE2       ', 'MEBA2DSE3       ', 'MEBA3DSE2       ', 'MEBA3DSE3       ' ]
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_GRILLE    = 6
    character(len=16),parameter :: ACE_EL_GRILLE(ACE_NB_GRILLE) =[ &
        'MEGCQU4         ', 'MEGMTR3         ', 'MEGMQU4         ','MEGMTR6         ', &
        'MEGMQU8         ', 'MEGCTR3         ' ]
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_MEMBRANE = 6
    character(len=16),parameter :: ACE_EL_MEMBRANE(ACE_NB_MEMBRANE) =[ &
        'MEMBTR3         ', 'MEMBTR6         ', 'MEMBQU4         ', 'MEMBQU8         ', &
        'MEMBTR7         ', 'MEMBQU9         ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_MASSIF = 61
    character(len=16), parameter :: ACE_EL_MASSIF(ACE_NB_MASSIF) = [ &
        'EIAXQS8         ', 'EIAXQU8         ', 'EIPLQS8         ', 'EIPLQU8         ', &
        'MEAXQU4         ', 'MEAXQU8         ', 'MEAXQU9         ', 'MEAXTR3         ', &
        'MEAXTR6         ', 'MECA_HEXA20     ', 'MECA_HEXA27     ', 'MECA_HEXA8      ', &
        'MECA_HEXS20     ', 'MECA_HEXS8      ', 'MECA_PENTA15    ', 'MECA_PENTA18    ', &
        'MECA_PENTA6     ', 'MECA_PYRAM13    ', 'MECA_PYRAM5     ', 'MECA_TETRA10    ', &
        'MECA_TETRA4     ', 'MECA_TETRS10    ', 'MECPQU4         ', 'MECPQU8         ', &
        'MECPQU9         ', 'MECPTR3         ', 'MECPTR6         ', 'MEDPQU4         ', &
        'MEDPQU8         ', 'MEDPQU9         ', 'MEDPTR3         ', 'MEDPTR6         ', &
        'MEEI_HEXA20     ', 'MEEI_HEXS20     ', 'MEEI_PENTA15    ', 'MEEI_PENTS15    ', &
        'MET3SEG3        ', 'MET3SEG4        ', 'MET6SEG3        ', 'THAXQU4         ', &
        'THAXQU8         ', 'THAXQU9         ', 'THAXTR3         ', 'THAXTR6         ', &
        'THER_HEXA20     ', 'THER_HEXA27     ', 'THER_HEXA8      ', 'THER_HEXA8_D    ', &
        'THER_PENTA15    ', 'THER_PENTA6     ', 'THER_PENTA6_D   ', 'THER_PYRAM13    ', &
        'THER_PYRAM5     ', 'THER_TETRA10    ', 'THER_TETRA4     ', 'THER_TETRA4_D   ', &
        'THPLQU4         ', 'THPLQU8         ', 'THPLQU9         ', 'THPLTR3         ', &
        'THPLTR6         ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_THHMM     = 75
    character(len=16),parameter :: ACE_EL_THHMM(ACE_NB_THHMM) = [ &
        'HM_DPQ8S        ', 'HM_AXIS_QU8S    ', 'HM_DPTR6S       ', 'HM_AXIS_TR6S    ', &
        'HM_HEXA20S      ', 'HM_PENTA15S     ', 'HM_TETRA10S     ', 'THM_DPQ8S       ', &
        'THM_AXIS_QU8S   ', 'THM_DPTR6S      ', 'THM_AXIS_TR6S   ', 'THM_HEXA20S     ', &
        'THM_PENTA15S    ', 'THM_TETRA10S    ', 'H_DPQ8S         ', 'H_DPTR6S        ', &
        'H_HEXA20S       ', 'H_PENTA15S      ', 'H_TETRA10S      ', 'THHM_DPQ8S      ', &
        'THHM_AXIS_QU8S  ', 'THHM_DPTR6S     ', 'THHM_AXIS_TR6S  ', 'THHM_HEXA20S    ', &
        'THHM_PENTA15S   ', 'THHM_TETRA10S   ', 'HHM_DPQ8S       ', 'HHM_AXIS_QU8S   ', &
        'HHM_DPTR6S      ', 'HHM_AXIS_TR6S   ', 'HHM_HEXA20S     ', 'HHM_PENTA15S    ', &
        'HHM_TETRA10S    ', 'THH_DPQ8S       ', 'THH_AXIS_QU8S   ', 'THH_DPTR6S      ', &
        'THH_AXIS_TR6S   ', 'THH_HEXA20S     ', 'THH_PENTA15S    ', 'THH_TETRA10S    ', &
        'HH_DPQ8S        ', 'HH_AXIS_QU8S    ', 'HH_DPTR6S       ', 'HH_AXIS_TR6S    ', &
        'HH_HEXA20S      ', 'HH_PENTA15S     ', 'HH_TETRA10S     ', 'THH2M_DPQ8S     ', &
        'THH2M_AXIS_QU8S ', 'THH2M_DPTR6S    ', 'THH2M_AXIS_TR6S ', 'THH2M_HEXA20S   ', &
        'THH2M_PENTA15S  ', 'THH2M_TETRA10S  ', 'HH2M_DPQ8S      ', 'HH2M_AXIS_QU8S  ', &
        'HH2M_DPTR6S     ', 'HH2M_AXIS_TR6S  ', 'HH2M_HEXA20S    ', 'HH2M_PENTA15S   ', &
        'HH2M_TETRA10S   ', 'THH2_DPQ8S      ', 'THH2_AXIS_QU8S  ', 'THH2_DPTR6S     ', &
        'THH2_AXIS_TR6S  ', 'THH2_HEXA20S    ', 'THH2_PENTA15S   ', 'THH2_TETRA10S   ', &
        'HH2_DPQ8S       ', 'HH2_AXIS_QU8S   ', 'HH2_DPTR6S      ', 'HH2_AXIS_TR6S   ', &
        'HH2_HEXA20S     ', 'HH2_PENTA15S    ', 'HH2_TETRA10S    ' ]
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_HHO     = 36
    character(len=16),parameter :: ACE_EL_HHO(ACE_NB_HHO) = [ &
        'MECA3DH27_HHO111', 'MECA3DH27_HHO222', 'MECA3DP19_HHO111', 'MECA3DP19_HHO222', &
        'MECA3DP21_HHO111', 'MECA3DP21_HHO222', 'MECA3DT15_HHO111', 'MECA3DT15_HHO222', &
        'MECA_DPQ9_HHO111', 'MECA_DPQ9_HHO222', 'MECA_DPT7_HHO111', 'MECA_DPT7_HHO222', &
        'THER2DQ9_HHO000 ', 'THER2DQ9_HHO111 ', 'THER2DQ9_HHO222 ', 'THER2DT7_HHO000 ', &
        'THER2DT7_HHO111 ', 'THER2DT7_HHO222 ', 'THER3DH27_HHO000', 'THER3DH27_HHO111', &
        'THER3DH27_HHO222', 'THER3DP19_HHO000', 'THER3DP19_HHO111', 'THER3DP19_HHO222', &
        'THER3DP21_HHO000', 'THER3DP21_HHO111', 'THER3DP21_HHO222', 'THER3DT15_HHO000', &
        'THER3DT15_HHO111', 'THER3DT15_HHO222', 'THERAXQ9_HHO000 ', 'THERAXQ9_HHO111 ', &
        'THERAXQ9_HHO222 ', 'THERAXT7_HHO000 ', 'THERAXT7_HHO111 ', 'THERAXT7_HHO222 ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_TYPE_ELEM = ACE_NB_POUTRE + ACE_NB_DISCRET  + &
                                  ACE_NB_COQUE  + ACE_NB_CABLE  + ACE_NB_BARRE       + &
                                  ACE_NB_MASSIF + ACE_NB_GRILLE + ACE_NB_MEMBRANE    + &
                                  ACE_NB_THHMM  + ACE_NB_HHO

    type elem_supp_carac
        character(len=16) :: catanom(ACE_NB_TYPE_ELEM)
        integer(kind=8)   :: acenum(ACE_NB_TYPE_ELEM)
        integer(kind=8)   :: catanum(ACE_NB_TYPE_ELEM)
        integer(kind=8)   :: aceind(ACE_NB_ELEMENT,2)
        integer(kind=8)   :: MaxCataNum
    end type elem_supp_carac

    type (elem_supp_carac) :: elem_supp

contains
!
! --------------------------------------------------------------------------------------------------
!   Retourne l'entier codé d'une liste d'entier []
!   Si dg_ est présent, on y ajoute les éléments de liste, et on retourne le tout
!       Pas d'ordre dans la liste d'entier, il peut y avoir des répétitions
    integer(kind=8) function ACE_CodeEntier(liste,dg_) result(dg)
        integer(kind=8),intent(in)              :: liste(:)
        integer(kind=8),optional, intent(in)    :: dg_
!
        integer(kind=8) :: ii, code
!
        dg = 0
        if ( present(dg_) ) dg = dg_
        do ii = lbound(liste,1), ubound(liste,1)
            code = lshift(1,liste(ii))
            dg   = ior(dg,code)
        enddo
    end function ACE_CodeEntier
!
! Indique si l'entier icode fait partie de l'entier codé irgcmp
    logical function ACE_IsInCode(icode, irgcmp)
        integer(kind=8), intent(in) :: icode, irgcmp
!
        integer(kind=8) :: code
!
        code = lshift(1,irgcmp)
        ACE_IsInCode = iand(icode,code) .eq. code
    end function ACE_IsInCode
!
!
! --------------------------------------------------------------------------------------------------
!
!   Codage avec un entier du type de section
!
!   Pour l'instant c'est "pas simple" on garde l'ancienne numération
!   Dans les composantes des PhysicalQuantity on trouve des TSEC et TVAR
!       Le codage est en "dur", on les remplace progressivement par
!           TSEC = ACE_SECTION_[ GENERALE | RECTANGLE | CERCLE ]
!           TVAR = ACE_SECTION_[ CONSTANTE | AFFINE | HOMOTHETIQUE ]
!       Comme c'est progressif on garde le même codage
!
!   Quand tout sera "bien fait", on pourra faire
!       code = FOR + VAR et ne garder qu'une composante FORVAR
!       FOR = Forme de la section
!                   GENERALE    RECTANGLE   CERCLE
!                   10          20          30
!       VAR = Variation des sections le long de l'axe
!                   CONSTANTE   AFFINE  HOMOTHETIQUE
!                   1           2       3
!       On gardera, en changeant ce qu'il faut
!           ACE_CodeFormeVaria
!           ACE_DeCodeNumFormeVaria
!           ACE_DeCodeNomFormeVaria
!
!
    integer(kind=8) function ACE_CodeFormeVaria(TF, TV) result(code)
!
#include "asterfort/assert.h"
!
        integer(kind=8), intent(in) :: TF, TV
!
        code = 0
!
        if (TF .eq. ACE_SECTION_GENERALE) then
            code = code + 10
        else if (TF .eq. ACE_SECTION_RECTANGLE) then
            code = code + 20
        else if (TF .eq. ACE_SECTION_CERCLE) then
            code = code + 30
        else
            write(*,*) 'ACE_CodeFormeVaria TF', TF
            ASSERT( .FALSE. )
        end if
!
        if (TV .eq. ACE_SECTION_CONSTANTE) then
            code = code + 1
        else if (TV .eq. ACE_SECTION_AFFINE) then
            code = code + 2
        else if (TV .eq. ACE_SECTION_HOMOTHETIQUE) then
            code = code + 3
        else
            write(*,*) 'ACE_CodeFormeVaria TV', TV
            ASSERT( .FALSE. )
        end if
!
    end function ACE_CodeFormeVaria
!
!
    subroutine ACE_DeCodeNumFormeVaria(TFV,resultat)
!
#include "asterfort/assert.h"
!
        integer(kind=8), intent(in) :: TFV
        integer(kind=8), intent(out) :: resultat(2)
!
        integer(kind=8) :: code
!
        code = TFV
        if ( (10 .lt. code).and.(code .lt. 20) ) then
            resultat(1) = ACE_SECTION_GENERALE
            code = code - 10
        else if ( (20 .lt. code).and.(code .lt. 30) ) then
            resultat(1) = ACE_SECTION_RECTANGLE
            code = code - 20
        else if ( (30 .lt. code).and.(code .lt. 40)) then
            resultat(1) = ACE_SECTION_CERCLE
            code = code - 30
        else
            write(*,*) 'ACE_DeCodeNumFormeVaria TFV ', TFV
            ASSERT( .FALSE. )
        end if
!
        if ( code .eq. 1) then
            resultat(2) = ACE_SECTION_CONSTANTE
        else if ( code .eq. 2) then
            resultat(2) = ACE_SECTION_AFFINE
        else if ( code .eq. 3) then
            resultat(2) = ACE_SECTION_HOMOTHETIQUE
        else
            write(*,*) 'ACE_DeCodeNumFormeVaria TFV, code ', TFV, code
            ASSERT( .FALSE. )
        end if
!
    end subroutine ACE_DeCodeNumFormeVaria
!
!
    subroutine ACE_DeCodeNomFormeVaria(TFV,resultat)
!
#include "asterfort/assert.h"
!
        integer(kind=8), intent(in) :: TFV
        character(len=*), intent(out) :: resultat(2)
!
        integer(kind=8) :: code
!
        code = TFV
        if ( (10 .lt. code).and.(code .lt. 20) ) then
            resultat(1) = 'GENERALE'
            code = code - 10
        else if ( (20 .lt. code).and.(code .lt. 30) ) then
            resultat(1) = 'RECTANGLE'
            code = code - 20
        else if ( (30 .lt. code).and.(code .lt. 40)) then
            resultat(1) = 'CERCLE'
            code = code - 30
        else
            write(*,*) 'ACE_DeCodeNomFormeVaria TFV ', TFV
            ASSERT( .FALSE. )
        end if
!
        if ( code .eq. 1) then
            resultat(2) = 'CONSTANTE'
        else if ( code .eq. 2) then
            resultat(2) = 'AFFINE'
        else if ( code .eq. 3) then
            resultat(2) = 'HOMOTHETIQUE'
        else
            write(*,*) 'ACE_DeCodeNomFormeVaria TFV, code ', TFV, code
            ASSERT( .FALSE. )
        end if
!
    end subroutine ACE_DeCodeNomFormeVaria
!
! --------------------------------------------------------------------------------------------------
!   Initialisation des éléments pouvant être affectés :
!       elem_supp%catanom   ACE_EL_[POUTRE|DISCRET|COQUE|CABLE|BARRE|GRILLE|MEMBRANE|MASSIF|THHMM]
!       elem_supp%acenum    ACE_NU_[ même liste] Classification de cara_elem
!       elem_supp%catanum   Leur numéro dans le catalogue de code_aster
!       elem_supp%aceind    Index de début des elem_supp%[catanom|catanum|acenum] l'éléments
!                               Exemple pour les DISCRET, ils vont de :
!                                       elem_supp%aceind(ACE_NU_DISCRET,1)
!                                       elem_supp%aceind(ACE_NU_DISCRET,2)
!
!
!   Il y a ACE_NB_ELEMENT type d'élément traité par AFFE_CARA_ELEM
!       Chaque ACE_NB_ELEMENT à plusieurs modélisations possibles
!           Pour les POUTRES : ACE_NU_POUTRE Numéro du type d'élément
!                              ACE_NB_POUTRE Nb modèles de poutres possibles
!                              ACE_EL_POUTRE(ii) : ii ème modèle
!                              ACE_EL_POUTRE(ACE_MECA_POU_D_T) = 'MECA_POU_D_T ' Nom catalogue Aster
!
!
!
    Subroutine ACE_Init_elem_affe()
!
#include "asterfort/assert.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
!
        integer(kind=8)   :: nbtel,ii
!
        nbtel = 0
        do ii = nbtel+1, nbtel+ACE_NB_POUTRE
            elem_supp%catanom(ii) = ACE_EL_POUTRE(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_POUTRE
        enddo
        elem_supp%aceind(ACE_NU_POUTRE,1) = nbtel+1
        elem_supp%aceind(ACE_NU_POUTRE,2) = nbtel+ACE_NB_POUTRE
        nbtel = nbtel + ACE_NB_POUTRE
!
        do ii = nbtel+1, nbtel+ACE_NB_DISCRET
            elem_supp%catanom(ii) = ACE_EL_DISCRET(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_DISCRET
        enddo
        elem_supp%aceind(ACE_NU_DISCRET,1) = nbtel+1
        elem_supp%aceind(ACE_NU_DISCRET,2) = nbtel+ACE_NB_DISCRET
        nbtel = nbtel + ACE_NB_DISCRET
!
        do ii = nbtel+1, nbtel+ACE_NB_COQUE
            elem_supp%catanom(ii) = ACE_EL_COQUE(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_COQUE
        enddo
        elem_supp%aceind(ACE_NU_COQUE,1) = nbtel+1
        elem_supp%aceind(ACE_NU_COQUE,2) = nbtel+ACE_NB_COQUE
        nbtel = nbtel + ACE_NB_COQUE
!
        do ii = nbtel+1, nbtel+ACE_NB_CABLE
            elem_supp%catanom(ii) = ACE_EL_CABLE(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_CABLE
        enddo
        elem_supp%aceind(ACE_NU_CABLE,1) = nbtel+1
        elem_supp%aceind(ACE_NU_CABLE,2) = nbtel+ACE_NB_CABLE
        nbtel = nbtel + ACE_NB_CABLE
!
        do ii = nbtel+1, nbtel+ACE_NB_BARRE
            elem_supp%catanom(ii) = ACE_EL_BARRE(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_BARRE
        enddo
        elem_supp%aceind(ACE_NU_BARRE,1) = nbtel+1
        elem_supp%aceind(ACE_NU_BARRE,2) = nbtel+ACE_NB_BARRE
        nbtel = nbtel + ACE_NB_BARRE
!
        do ii = nbtel+1, nbtel+ACE_NB_MASSIF
            elem_supp%catanom(ii) = ACE_EL_MASSIF(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_MASSIF
        enddo
        elem_supp%aceind(ACE_NU_MASSIF,1) = nbtel+1
        elem_supp%aceind(ACE_NU_MASSIF,2) = nbtel+ACE_NB_MASSIF
        nbtel = nbtel + ACE_NB_MASSIF
!
        do ii = nbtel+1, nbtel+ACE_NB_GRILLE
            elem_supp%catanom(ii) = ACE_EL_GRILLE(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_GRILLE
        enddo
        elem_supp%aceind(ACE_NU_GRILLE,1) = nbtel+1
        elem_supp%aceind(ACE_NU_GRILLE,2) = nbtel+ACE_NB_GRILLE
        nbtel = nbtel + ACE_NB_GRILLE
!
        do ii = nbtel+1, nbtel+ACE_NB_MEMBRANE
            elem_supp%catanom(ii) = ACE_EL_MEMBRANE(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_MEMBRANE
        enddo
        elem_supp%aceind(ACE_NU_MEMBRANE,1) = nbtel+1
        elem_supp%aceind(ACE_NU_MEMBRANE,2) = nbtel+ACE_NB_MEMBRANE
        nbtel = nbtel + ACE_NB_MEMBRANE
!
        do ii = nbtel+1, nbtel+ACE_NB_THHMM
            elem_supp%catanom(ii) = ACE_EL_THHMM(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_THHMM
        enddo
        elem_supp%aceind(ACE_NU_THHMM,1) = nbtel+1
        elem_supp%aceind(ACE_NU_THHMM,2) = nbtel+ACE_NB_THHMM
        nbtel = nbtel + ACE_NB_THHMM
!
        do ii = nbtel+1, nbtel+ACE_NB_HHO
            elem_supp%catanom(ii) = ACE_EL_HHO(ii-nbtel)
            elem_supp%acenum(ii) = ACE_NU_HHO
        enddo
        elem_supp%aceind(ACE_NU_HHO,1) = nbtel+1
        elem_supp%aceind(ACE_NU_HHO,2) = nbtel+ACE_NB_HHO
        nbtel = nbtel + ACE_NB_HHO
!
        ASSERT( nbtel .eq. ACE_NB_TYPE_ELEM )
!       Récupération des numéros des types éléments
        do ii = 1, ACE_NB_TYPE_ELEM
            call jenonu(jexnom('&CATA.TE.NOMTE', elem_supp%catanom(ii)), elem_supp%catanum(ii))
        enddo
        elem_supp%MaxCataNum = maxval(elem_supp%catanum(:))
    end subroutine ACE_Init_elem_affe

!&>
end module cara_elem_parameter_module
