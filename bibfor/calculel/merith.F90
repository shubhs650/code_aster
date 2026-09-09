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
subroutine merith(modelZ, loadNameZ, materCodeZ, caraElemZ, &
                  timeMapZ, matrElemZ, jvBaseZ)
!
    implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/memare.h"
#include "asterfort/merit1.h"
#include "asterfort/merit2.h"
#include "asterfort/reajre.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: modelZ, loadNameZ, materCodeZ, caraElemZ, jvBaseZ, timeMapZ
    character(len=*), intent(inout) :: matrElemZ
!
! --------------------------------------------------------------------------------------------------
!
! Pseudo-Thermic for fluid
!
! CALCUL DES MATRICES ELEMENTAIRES DE RIGIDITE THERMIQUE
!
! --------------------------------------------------------------------------------------------------
!
!
!     LES NOMS QUI SUIVENT SONT LES PREFIXES UTILISATEUR K8:
!        MODELZ : NOM DU MODELE
!        NCHAR  : NOMBRE DE CHARGES
!        LCHAR  : LISTE DES CHARGES
!        MATE   : CHAMP DE MATERIAUX
!        CARAZ  : CHAMP DE CARAC_ELEM
!        TIMEZ  : CHAMPS DE TEMPSR
!        MATELZ : NOM  DU  MATELE (N RESUELEM) PRODUIT
!                  ( ISO      , 'RIGIDI_TH'             )
!                  ( CAL_TI   , 'DDLMUR_THER'           )
!                  ( ISO_FACE , 'RIGIDI_TH_COEFHR/F'    )
!        NH     : NUMERO DE L'HARMONIQUE DE FOURIER(SI PAS FOURIER NH=0)
!
!     SORTIES:
!        MATELZ   : LE MATELE EST REMPLI.
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iResuElem, iret, nbResuElem1, nbResuElem2, indxMatrElem
    character(len=8) :: model, caraElem, loadName
    character(len=19) :: matrElem
    character(len=19), parameter :: matrElem1 = '&MERITH1', matrElem2 = '&MERITH2'
    character(len=1) :: jvBase
    character(len=24) :: timeMap, resuElemPref, materCode
    character(len=24), pointer :: resuElem1(:) => null(), resuElem2(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    model = modelZ
    materCode = materCodeZ
    caraElem = caraElemZ
    timeMap = timeMapZ
    jvBase = jvBaseZ
    matrElem = matrElemZ
    loadName = loadNameZ
!
!     -- RIGIDITE CORRESPONDANT AUX ELEMENTS ISO ET AUX ELEMENTS CAL_TI:
    resuElemPref = '&MERITH1'
    indxMatrElem = 0
    call merit1(model, caraElem, materCode, &
                loadName, &
                timeMap, matrElem1, resuElemPref, &
                indxMatrElem, jvBase)
    call jeexin(matrElem1(1:19)//'.RELR', iret)
    nbResuElem1 = 0
    if (iret .ne. 0) then
        call jelira(matrElem1(1:19)//'.RELR', 'LONUTI', nbResuElem1)
        call jeveuo(matrElem1(1:19)//'.RELR', 'L', vk24=resuElem1)
    end if
!
!     -- RIGIDITE CORRESPONDANT AUX ELEMENTS D'ECHANGE:
    resuElemPref = '&MERITH2'
    indxMatrElem = nbResuElem1
    call merit2(model, caraElem, &
                loadName, &
                timeMap, matrElem2, resuElemPref, &
                indxMatrElem, jvBase)
    call jeexin(matrElem2(1:19)//'.RELR', iret)
    nbResuElem2 = 0
    if (iret .ne. 0) then
        call jelira(matrElem2(1:19)//'.RELR', 'LONUTI', nbResuElem2)
        call jeveuo(matrElem2(1:19)//'.RELR', 'L', vk24=resuElem2)
    end if

! - Create MATR_ELEM
    call memare(jvBase, matrElem, model, 'RIGI_THER')

! - Add RESU_ELEM
    do iResuElem = 1, nbResuElem1
        call reajre(matrElem, resuElem1(iResuElem), jvBase)
    end do
    do iResuElem = 1, nbResuElem2
        call reajre(matrElem, resuElem2(iResuElem), jvBase)
    end do

! - Clean
    call jedetr(matrElem1(1:19)//'.RELR')
    call jedetr(matrElem1(1:19)//'.RERR')
    call jedetr(matrElem2(1:19)//'.RELR')
    call jedetr(matrElem2(1:19)//'.RERR')

    matrElemZ = matrElem
    call jedema()
end subroutine
