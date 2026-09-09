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
subroutine racotu(iprno, lonlis, klisno, noepou, mesh, &
                  ligrel, model, caraElem, numddl, &
                  lisrel, coorig)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/afretu.h"
#include "asterfort/assvec.h"
#include "asterfort/calcul.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/normev.h"
#include "asterfort/raorfi.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
#include "asterfort/vemare.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8) :: lonlis, iprno(*)
    character(len=8) :: klisno(lonlis), noepou, mesh, caraElem, model
    character(len=14) :: numddl
    character(len=19) :: ligrel, lisrel
    real(kind=8) :: coorig(3)
!
! --------------------------------------------------------------------------------------------------
!
!     RACCORD COQUE_TUYAU PAR DES RELATIONS LINEAIRES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 3
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    integer(kind=8), parameter :: nbmode = 3, nbcmp = 6*(nbmode-1)
    integer(kind=8) :: numno1
    character(len=8) :: nocmp(nbcmp), nomddl(4)
    character(len=24) :: valech
    real(kind=8) :: coef(4), eg1(3), eg2(3), eg3(3)
    real(kind=8) :: rayon, coori1(3), gp1(3)
    integer(kind=8) :: imod, info, ifm, idch1
    integer(kind=8) :: iwi1wo1, k
    integer(kind=8) :: nbcoef, idec, ival, nbec, ino, i
    real(kind=8), pointer :: vale(:) => null()
    character(len=24), parameter :: numeModeField = '&&RACOTU.NUME_MODE'
    character(len=24), parameter :: mapPipeAxis = '&&RACOTU.CAXE_TUY'
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, info)

! - Initialisations
    lchin = ' '
    lpain = ' '
    lchout = ' '
    lpaout = ' '

! - CALCUL DU RAYON DU MAILLAGE COQUE A L'AIDE DU PREMIER N
    call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=vale)
    numno1 = char8_to_int(klisno(1))
    coori1(1) = vale(3*(numno1-1)+1)
    coori1(2) = vale(3*(numno1-1)+2)
    coori1(3) = vale(3*(numno1-1)+3)
    gp1 = coorig-coori1
    call normev(gp1, rayon)
!
!     CREATION D'UNE CARTE CONTENANT LE POINT P ORIGINE DE PHI
!
    call raorfi(mesh, ligrel, noepou, caraElem, coorig, &
                eg1, eg2, eg3, '&&RACOTU', rayon)
!
! --- DETERMINATION DE 3 LISTES  DE VECTEURS PAR ELEMENT PRENANT
! --- LEURS VALEURS AUX NOEUDS DES ELEMENTS. CES VALEURS SONT :
!     VECT_COEF_UM
! --- 1/PI* SOMME/S_ELEMENT(COS(M.PHI).P(1,J).NI)DS
! --- 1/PI* SOMME/S_ELEMENT(SIN(M.PHI).P(1,J).NI)DS
!     VECT_COEF_VM
! --- 1/PI* SOMME/S_ELEMENT(COS(M.PHI).P(2,J).NI)DS
! --- 1/PI* SOMME/S_ELEMENT(SIN(M.PHI).P(2,J).NI)DS
!     VECT_COEF_WM
! --- 1/PI* SOMME/S_ELEMENT(COS(M.PHI).P(3,J).NI)DS
! --- 1/PI* SOMME/S_ELEMENT(SIN(M.PHI).P(3,J).NI)DS
!
! --- OU P EST LA MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
! --- (E1,E2,E3) DEFINI SUR LE BORD ORIENTE DE LA COQUE
!     ------------------------------
! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = mesh//'.COORDO'
    lpain(2) = 'PORIGIN'
    lchin(2) = '&&RAPOCO.CAORIGE'
    lpain(3) = 'PORIGFI'
    lchin(3) = '&&RACOTU.CAORIFI'
    lpain(4) = 'PNUMMOD'
    lchin(4) = numeModeField(1:19)
    nbFieldIn = 4

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem, caorienZ_=mapPipeAxis)

! - Add output fields
    lpaout(1) = 'PVECTU1'
    lchout(1) = '&&RAPOTU.COEF_UM'
    lpaout(2) = 'PVECTU2'
    lchout(2) = '&&RAPOTU.COEF_VM'
    lpaout(3) = 'PVECTU3'
    lchout(3) = '&&RAPOTU.COEF_WM'

! - CREATION DES .RERR DES VECTEURS EN SORTIE DE CALCUL
    call vemare('V', '&&RAPOTU', model)

! - CREATION DU .RELR
    call jedetr('&&RAPOTU           .RELR')
    call reajre('&&RAPOTU', ' ', 'V')
!
!  #################################################################
!
!
!     RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE DDL WO
!
!
!  #################################################################
!

    imod = 0
    call mecact('V', numeModeField, 'LIGREL', ligrel, 'NUMMOD', &
                ncmp=1, nomcmp='NUM', si=imod)
    call calcul('S', 'CARA_SECT_POUT5', ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')
    call jedetr('&&RAPOTU           .RELR')
    call reajre('&&RAPOTU', lchout(3), 'V')
    call assvec('V', 'CH_DEPL_3', 1, '&&RAPOTU           .RELR', [1.d0], numddl)
    valech = 'CH_DEPL_3          .VALE'
    nbcoef = 1
    idec = 0
    nomddl(1) = 'WO'
    coef(1) = -2.d0
!
    call afretu(iprno, lonlis, klisno, noepou, mesh, &
                valech, nbcoef, idec, coef, nomddl, &
                lisrel)
!
!
!  #################################################################
!
!   CAS PARTICULIER, POUR GERER LES RELATIONS SUR WI1 et WO1
!
!    Compte tenu de la relation particulière entre WI1 et VO1 d'une part
!    et WO1 et VI1 d'autre part, la relation cinematique ne peut pas porter que sur
!    la projection de WO1 ou WI1 sur cos(phi) ou sin(phi). Il faut aussi prendre
!    en compte les projections de VO1 et VI1 sur cos(phi) et sin(phi).
!
!    Il faut donc construire le champ WI1+VO1 d'une part, et WO1-VI1 d'autre part
!
!  #################################################################
!
!   On calcul les projections sur U1, V1 et W1 - imod = 1
    imod = 1
    if (info .eq. 2) then
        write (ifm, *) 'RELATIONS SUR LE MODE ', imod
    end if
    call mecact('V', numeModeField, 'LIGREL', ligrel, 'NUMMOD', &
                ncmp=1, nomcmp='NUM', si=imod)
    call calcul('S', 'CARA_SECT_POUT5', ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')
!
    call jedetr('&&RAPOTU           .RELR')
    call reajre('&&RAPOTU', lchout(3), 'V')
    call assvec('V', 'CH_DEPL_3', 1, '&&RAPOTU           .RELR', [1.d0], numddl)
    valech = 'CH_DEPL_3          .VALE'
!
!   On recupere les champs WO1 et WI1, stockes dans 'CH_DEPL_3          .VALE'
!
    call jeveuo(valech, 'L', idch1)
    call dismoi('NB_EC', 'DEPL_R', 'GRANDEUR', repi=nbec)
    if (nbec .gt. 11) then
        call utmess('F', 'MODELISA_94')
    end if
!
!   On alloue un vecteur spécifique qui recueuilleura WI1+VO1 et WO1-VI1
!   Et on recopie WI1 et WO1 dedans
!
    call wkvect('&&RACOTU.WI1WO1', 'V V R', lonlis*6, iwi1wo1)
    do i = 1, lonlis
        ino = char8_to_int(klisno(i))
!           ADRESSE DE LA PREMIERE CMP DU NOEUD INO DANS LES CHAMNO
        ival = iprno((ino-1)*(nbec+2)+1)
        do k = 1, 6
            zr(iwi1wo1+6*(i-1)+(k-1)) = zr(idch1+ival-1+(k-1))
        end do
    end do
!
!   On construit et on assemble les champs VO1 et VI1
!
    call jedetr('&&RAPOTU           .RELR')
    call reajre('&&RAPOTU', lchout(2), 'V')
    call assvec('V', 'CH_DEPL_2', 1, '&&RAPOTU           .RELR', [1.d0], numddl)
    valech = 'CH_DEPL_2          .VALE'
!   On recupere ces champs, et on les combine avec les champs WI1 et WO1
    call jeveuo(valech, 'L', idch1)
    do i = 1, lonlis
        ino = char8_to_int(klisno(i))
!           ADRESSE DE LA PREMIERE CMP DU NOEUD INO DANS LES CHAMNO
        ival = iprno((ino-1)*(nbec+2)+1)
!       Construction de wi1 + vo1
        zr(iwi1wo1+6*(i-1)+0) = zr(iwi1wo1+6*(i-1)+0)+zr(idch1+ival-1+3)
        zr(iwi1wo1+6*(i-1)+1) = zr(iwi1wo1+6*(i-1)+1)+zr(idch1+ival-1+4)
        zr(iwi1wo1+6*(i-1)+2) = zr(iwi1wo1+6*(i-1)+2)+zr(idch1+ival-1+5)
!       Construction de wo1-vi1
        zr(iwi1wo1+6*(i-1)+3) = zr(iwi1wo1+6*(i-1)+3)-zr(idch1+ival-1+0)
        zr(iwi1wo1+6*(i-1)+4) = zr(iwi1wo1+6*(i-1)+4)-zr(idch1+ival-1+1)
        zr(iwi1wo1+6*(i-1)+5) = zr(iwi1wo1+6*(i-1)+5)-zr(idch1+ival-1+2)
!
    end do
!
!   On "detourne" l'utilisation de idec pour traiter le cas particulier WI1
    idec = -1
!
    nbcoef = 1
    nomddl(1) = 'WI1'
    coef(1) = -1.d0
!
    valech = '&&RACOTU.WI1WO1         '
    call afretu(iprno, lonlis, klisno, noepou, mesh, &
                valech, nbcoef, idec, coef, nomddl, &
                lisrel)
!
!   On "detourne" l'utilisation de idec pour traiter le cas particulier WO1
    idec = -2
!
    nbcoef = 1
    nomddl(1) = 'WO1'
    coef(1) = -1.d0
    call afretu(iprno, lonlis, klisno, noepou, mesh, &
                valech, nbcoef, idec, coef, nomddl, &
                lisrel)

!
!  #################################################################
!
!
!   RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE
!   DDL UIM, UOM, VIM, VOM, WIM, WOM, M VARIANT DE 2 A NBMODE
!
!
!  #################################################################
!
    nocmp(1) = 'UI2'
    nocmp(2) = 'VI2'
    nocmp(3) = 'WI2'
    nocmp(4) = 'UO2'
    nocmp(5) = 'VO2'
    nocmp(6) = 'WO2'
    nocmp(7) = 'UI3'
    nocmp(8) = 'VI3'
    nocmp(9) = 'WI3'
    nocmp(10) = 'UO3'
    nocmp(11) = 'VO3'
    nocmp(12) = 'WO3'
!
!     IDEC=0 SIGNIFIE QUE ON UTILISE LES TERMES EN COS(M*PHI)
!     IDEC=3 SIGNIFIE QUE ON UTILISE LES TERMES EN SIN(M*PHI)
!
    do imod = 2, nbmode
        if (info .eq. 2) then
            write (ifm, *) 'RELATIONS SUR LE MODE ', imod
        end if
        call mecact('V', numeModeField, 'LIGREL', ligrel, 'NUMMOD', &
                    ncmp=1, nomcmp='NUM', si=imod)
        call calcul('S', 'CARA_SECT_POUT5', ligrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    'V', 'OUI')
        call jedetr('&&RAPOTU           .RELR')
        call reajre('&&RAPOTU', lchout(1), 'V')
        call assvec('V', 'CH_DEPL_1', 1, '&&RAPOTU           .RELR', [1.d0], numddl)
        valech = 'CH_DEPL_1          .VALE'
!
!        RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE DDL UIM
!
        idec = 0
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+1)
        coef(1) = -1.d0
        call afretu(iprno, lonlis, klisno, noepou, mesh, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!        RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE DDL UOM
!
        idec = 3
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+4)
        coef(1) = -1.d0
        call afretu(iprno, lonlis, klisno, noepou, mesh, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!        RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE DDL VOM
!
        call jedetr('&&RAPOTU           .RELR')
        call reajre('&&RAPOTU', lchout(2), 'V')
        call assvec('V', 'CH_DEPL_2', 1, '&&RAPOTU           .RELR', [1.d0], numddl)
        valech = 'CH_DEPL_2          .VALE'
        idec = 0
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+5)
        coef(1) = -1.d0
        call afretu(iprno, lonlis, klisno, noepou, mesh, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!        RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE DDL VIM
!        IDEC=3 SIGNIFIE QUE ON UTILISE LES TERMES EN SIN(M*PHI)
!
        idec = 3
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+2)
        coef(1) = -1.d0
        call afretu(iprno, lonlis, klisno, noepou, mesh, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!        RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE DDL WIM
!
        call jedetr('&&RAPOTU           .RELR')
        call reajre('&&RAPOTU', lchout(3), 'V')
        call assvec('V', 'CH_DEPL_3', 1, '&&RAPOTU           .RELR', [1.d0], numddl)
        valech = 'CH_DEPL_3          .VALE'
        idec = 0
!
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+3)
        coef(1) = -1.d0
!
        call afretu(iprno, lonlis, klisno, noepou, mesh, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!        RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE DDL WOM
!
        idec = 3
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+6)
        coef(1) = -1.d0
!
        call afretu(iprno, lonlis, klisno, noepou, mesh, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
    end do
!
! --- DESTRUCTION DES OBJETS DE TRAVAIL
!
    call jedetr('&&RACOTU.WI1WO1         ')
    call jedetr('&&RAPOTU           .RERR')
    call jedetr('&&RAPOTU           .RERR')
    call detrsd('CARTE', numeModeField)
    call detrsd('RESUELEM', '&&RAPOTU.COEF_UM')
    call detrsd('RESUELEM', '&&RAPOTU.COEF_VM')
    call detrsd('RESUELEM', '&&RAPOTU.COEF_WM')
    call detrsd('CHAMP_GD', mapPipeAxis)
    call detrsd('CHAMP_GD', '&&RACOTU.CAORIFI')
    call detrsd('CHAMP_GD', 'CH_DEPL_1')
    call detrsd('CHAMP_GD', 'CH_DEPL_2')
    call detrsd('CHAMP_GD', 'CH_DEPL_3')
!
    call jedema()
end subroutine
