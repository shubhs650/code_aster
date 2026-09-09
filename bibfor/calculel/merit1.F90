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
subroutine merit1(modelZ, caraElemZ, materCodeZ, &
                  loadNameZ, &
                  timeMap, matrElem, resuElemPref, &
                  indxMatrElem, jvBase)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ, materCodeZ
    character(len=*), intent(in) :: loadNameZ
    character(len=24), intent(in) :: timeMap
    character(len=19), intent(in) :: matrElem, resuElemPref
    integer(kind=8), intent(in) :: indxMatrElem
    character(len=1), intent(in) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
! Pseudo-Thermic for fluid
!
! CALCUL DES MATRICES ELEMENTAIRES DE RIGIDITE THERMIQUE
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREES:
!
!     LES NOMS QUI SUIVENT SONT LES PREFIXES UTILISATEUR K8:
!        MODELE : NOM DU MODELE
!        NCHAR  : NOMBRE DE CHARGES
!        LCHAR  : LISTE DES CHARGES
!        MATE   : CHAMP DE MATERIAUX
!        CARA   : CHAMP DE CARAC_ELEM
!        MATEL  : NOM DU MATR_ELEM (N RESUELEM) PRODUIT
!        PREFCH : PREFIXE DES NOMS DES RESUELEM STOCKES DANS MATEL
!        NUMERO : NUMERO D'ORDRE A PARTIR DUQUEL ON NOMME LES RESUELEM
!        TIME   : CHAMPS DE TEMPSR
!
!     SORTIES:
!        MATEL  : EST REMPLI.
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOut)
    integer(kind=8), parameter :: numeHarm = 0
    character(len=8) :: model, caraElem, loadName
    character(len=16) :: option
    character(len=24) :: chgeom, chharm
    character(len=24) :: modelLigrel, loadLigrel
    integer(kind=8) :: ilires, iret, nbFieldIn
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    model = modelz
    caraElem = caraElemZ
    loadName = loadNameZ
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! - Get geometry field
    call megeom(modelZ, chgeom)

! - Create field for Fourier
    call meharm(modelZ, numeHarm, chharm)
!
    call jeexin(matrElem//'.RERR', iret)
    if (iret .gt. 0) then
        call jedetr(matrElem//'.RERR')
        call jedetr(matrElem//'.RELR')
    end if
    call memare('V', matrElem, modelZ, 'RIGI_THER')

! - Set output field
    lpaout(1) = 'PMATTTR'
    lchout(1) = resuElemPref(1:8)//'.ME000'

! - Compute volumic terms
    ilires = 0
    if (model .ne. ' ') then
! ----- Add input fields
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PMATERC'
        lchin(2) = materCodeZ(1:24)
        lpain(3) = 'PINSTR'
        lchin(3) = timeMap
        lpain(4) = 'PHARMON'
        lchin(4) = chharm
        nbFieldIn = 4

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElemZ)

        option = 'RIGI_THER'
        ilires = ilires+1
        call codent(ilires+indxMatrElem, 'D0', lchout(1) (12:14))
        call calcul('S', option, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')
        call reajre(matrElem, lchout(1), jvBase)
    end if

! - Compute load terms
    if (loadName .ne. " ") then
        call exisd('CHAMP_GD', loadName//'.CHTH.CMULT', iret)
        if (iret .ne. 0) then
            lpain(1) = 'PDDLMUR'
            lchin(1) = loadName//'.CHTH.CMULT'
            nbFieldIn = 1
            ilires = ilires+1
            call codent(ilires+indxMatrElem, 'D0', lchout(1) (12:14))
            loadLigrel = loadName//'.CHTH.LIGRE'
            option = 'THER_DDLM_R'
            call calcul('S', option, loadLigrel, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, jvBase, &
                        'OUI')
            call reajre(matrElem, lchout(1), jvBase)
        end if
    end if
!
    call jedema()
end subroutine
