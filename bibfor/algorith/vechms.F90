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
subroutine vechms(model, materField, materCode, caraElem, varplu, listLoad, &
                  partps, vectElem)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/exixfe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/lisico.h"
#include "asterfort/lislch.h"
#include "asterfort/lislco.h"
#include "asterfort/lisnbg.h"
#include "asterfort/lisnnb.h"
#include "asterfort/lisnol.h"
#include "asterfort/vechmp.h"
#include "asterfort/vechmx.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: materField, caraElem, materCode
    real(kind=8), intent(in) :: partps(3)
    character(len=19), intent(in) :: listLoad, varplu
    character(len=19), intent(in) :: vectElem
!
! --------------------------------------------------------------------------------------------------
!
! CALCUL DES VECTEURS ELEMENTAIRES DES CHARGEMENTS MECANIQUES
! DE NEUMANN STANDARD (VOIR DEFINITION DANS LISDEF)
!
! CALCUL EFFECTIF - BOUCLE SUR LES CHARGES
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMO   : NOM DU MODELE
! IN  LISCHA : SD LISTE DES CHARGES
! IN  PARTPS : TABLEAU DONNANT T+, DELTAT ET THETA (POUR LE THM)
! IN  CARELE : CARACTERISTIQUES DES POUTRES ET COQUES
! IN  MATE   : MATERIAU CODE
! IN  VARPLU : VARIABLES DE COMMANDE A L'INSTANT T+
! IN  nbFieldInMax   : NOMBRE MAXI DE CHAMPS D'ENTREE
! IN  LPAIN  : LISTE DES PARAMETRES IN
! IN  LCHIN  : LISTE DES CHAMPS IN
! IN  LASTIN : NOMBRE EFFECTIF DE CHAMPS IN
! OUT VECELE : VECT_ELEM RESULTAT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100
    character(len=8) :: lpain(nbFieldInMax)
    character(len=19) :: lchin(nbFieldInMax)
    integer(kind=8) :: iLoad, nbLoad, nbFieldIn
    character(len=24), parameter :: listLoadIndxJv = '&&NOMLIS'
    integer(kind=8) :: genrec, ier
    aster_logical :: lneum, lXFEM
    integer(kind=8) :: nbLoadIndx, nbneum
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - INITIALISATIONS
    call exixfe(model, ier)
    lXFEM = ier .ne. 0
    call detrsd('VECT_ELEM', vectElem)

! - Get parameters for loads
    call lisnnb(listLoad, nbLoad)
    nbneum = lisnbg(listLoad, 'NEUM_MECA')

    if (nbneum .gt. 0) then
! ----- Set input fields
        call vechmp(model, materField, materCode, caraElem, &
                    varplu, lXFEM, partps, &
                    nbFieldInMax, lpain, lchin, nbFieldIn)

! ----- LISTE DES INDEX DES CHARGES
        call lisnol(listLoad, 'NEUM_MECA', listLoadIndxJv, nbLoadIndx)
        ASSERT(nbLoadIndx .gt. 0)

        do iLoad = 1, nbLoad
            call lislco(listLoad, iLoad, genrec)
            lneum = lisico('NEUM_MECA', genrec)
            if (lneum) then
                call vechmx(model, listLoad, iLoad, nbLoadIndx, listLoadIndxJv, &
                            nbFieldInMax, lpain, lchin, nbFieldIn, vectElem)
            end if
        end do
    end if
!
    call jedetr(listLoadIndxJv)
!
    call jedema()
end subroutine
