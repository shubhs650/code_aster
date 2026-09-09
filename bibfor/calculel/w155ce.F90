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
subroutine w155ce(resultOut, resultIn, nbStore, listStore)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterc/getfac.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/exlima.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/rsexch.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: resultOut, resultIn
    integer(kind=8), intent(in) :: nbStore, listStore(nbStore)
!
! --------------------------------------------------------------------------------------------------
!
!     COMMANDE :  POST_CHAMP / COQU_EXCENT
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyw = 'COQU_EXCENT'
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iret, iStore, numeStore, ibid, nocc, iocc
    character(len=8) :: model, caraElem, mplan, modelSave
    character(len=4) :: physQuanScal
    character(len=16) :: fieldName
    character(len=19) :: fieldIn, fieldOut, ligrel
    integer(kind=8) :: vali(2), iexi
    aster_logical :: lreel, lnoeu, ldetli, lvide
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)

!
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "
!
!
!     -- 1. : Y-A-T-IL QUELQUE CHOSE A FAIRE ?
!     ----------------------------------------
    call getfac('COQU_EXCENT', nocc)
    if (nocc .eq. 0) then
        goto 30
    end if
    ASSERT(nocc .lt. 10)
!
!
    modelSave = ' '
    ldetli = ASTER_FALSE
    lvide = ASTER_TRUE
    do iocc = 1, nocc
!
!     -- 2.  : NOMSYM, MPLAN :
!     --------------------------------------------------

        call getvtx(factorKeyw, 'NOM_CHAM', iocc=iocc, scal=fieldName, nbret=ibid)
        ASSERT(fieldName .eq. 'EFGE_ELNO' .or. fieldName .eq. 'EFGE_ELGA')
        call getvtx(factorKeyw, 'MODI_PLAN', iocc=iocc, scal=mplan, nbret=ibid)
        ASSERT(mplan .eq. 'OUI')
        lnoeu = fieldName .eq. 'EFGE_ELNO'
!
!
!     -- 3. : BOUCLE SUR LES NUMERO D ORDRE
!     --------------------------------------------------
        do iStore = 1, nbStore
            numeStore = listStore(iStore)
            call rsexch(' ', resultIn, fieldName, numeStore, fieldIn, iret)
            if (iret .eq. 0) then
!
!         -- 3.1 : MODELE, CARELE, LIGREL :
                call rslesd(resultIn, numeStore, model_=model, cara_elem_=caraElem)
                if (model .ne. modelSave) then
                    if (ldetli) then
                        call detrsd('LIGREL', ligrel)
                    end if
                    call exlima(' ', 1, 'G', model, ligrel)
                    modelSave = model
                    ldetli = ASTER_FALSE
                    if (ligrel(1:8) .ne. model) then
                        ldetli = ASTER_TRUE
                    end if
                end if
!
                call rsexch(' ', resultOut, fieldName, numeStore, fieldOut, iret)
                ASSERT(iret .eq. 100)
!
                call jelira(fieldIn//'.CELV', 'TYPE', cval=physQuanScal)
                if (physQuanScal .eq. 'R') then
                    lreel = ASTER_TRUE
                else if (physQuanScal .eq. 'C') then
                    lreel = ASTER_FALSE
                else
                    ASSERT(ASTER_FALSE)
                end if

! ------------- Add input fields
                nbFieldIn = 1
                if (lnoeu) then
                    if (lreel) then
                        lpain(nbFieldIn) = 'PEFFONR'
                    else
                        lpain(nbFieldIn) = 'PEFFONC'
                    end if
                else
                    if (lreel) then
                        lpain(nbFieldIn) = 'PEFFOGR'
                    else
                        lpain(nbFieldIn) = 'PEFFOGC'
                    end if
                end if
                lchin(nbFieldIn) = fieldIn

! ------------- Add fields for orientation
                call setOrieFields(nbFieldInMax, lpain, lchin, &
                                   nbFieldIn, caraElem)

! ------------- Set output fields
                if (lnoeu) then
                    if (lreel) then
                        lpaout(1) = 'PEFFOENR'
                    else
                        lpaout(1) = 'PEFFOENC'
                    end if
                else
                    if (lreel) then
                        lpaout(1) = 'PEFFOEGR'
                    else
                        lpaout(1) = 'PEFFOEGC'
                    end if
                end if
                lchout(1) = fieldOut

! ------------- Compute
                call calcul('C', 'EFGE_EXCENT', ligrel, &
                            nbFieldIn, lchin, lpain, &
                            nbFieldOut, lchout, lpaout, &
                            'G', 'OUI')
!
                call jeexin(lchout(1)//'.CELV', iexi)
                if (iexi .eq. 0) then
                    vali(1) = iocc
                    vali(2) = numeStore
                    call utmess('A', 'CALCULEL2_19', ni=2, vali=vali)
                else
                    ldetli = ASTER_FALSE
                    lvide = ASTER_FALSE
                    call rsnoch(resultOut, fieldName, numeStore)
                end if
            end if
        end do
    end do
!
    if (ldetli) then
        call detrsd('LIGREL', ligrel)
    end if
    if (lvide) then
        call utmess('F', 'CALCULEL2_20')
    end if
!
30  continue
    call jedema()
end subroutine
