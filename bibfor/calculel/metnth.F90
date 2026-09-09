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
subroutine metnth(model, loadNameJv, loadInfoJv, &
                  caraElem, materCodeZ, &
                  timeMap, tempPrev, matrElem)
!
    use coorSyst_module, only: setOrieFields
    use loadTherCompute_module
    use loadTherCompute_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/load_list_info.h"
#include "asterfort/megeom.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: loadNameJv, loadInfoJv
    character(len=8), intent(in) :: caraElem
    character(len=*), intent(in) :: materCodeZ
    character(len=24), intent(in) :: timeMap
    character(len=24), intent(in) :: tempPrev
    character(len=24), intent(in) :: matrElem
!
! --------------------------------------------------------------------------------------------------
!
! Thermic - Matrix
!
! Elementary matrix for convection (volumic and surfacic terms)
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: l_stat = ASTER_TRUE
    real(kind=8) :: theta
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    character(len=16), parameter :: option = 'RIGI_THER_CONV'
    character(len=24) :: chgeom
    character(len=24) :: chvite, modelLigrel, loadField, resuElem
    integer(kind=8) :: iret, iconv
    integer(kind=8) :: nbLoad, iLoad, nbFieldIn
    aster_logical :: noLoadInList
    character(len=8) :: loadName
    character(len=24), pointer :: listLoadName(:) => null()
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    character(len=8), pointer :: loadFieldVale(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    ASSERT(model .ne. ' ')
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "

! - Stationnary !
    ASSERT(l_stat)
    theta = r8vide()

! - Get loads
    call load_list_info(noLoadInList, nbLoad, listLoadName, listLoadInfo, &
                        loadNameJv, loadInfoJv)

! - Geometry field
    call megeom(model, chgeom)

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = materCodeZ
    lpain(3) = 'PINSTR'
    lchin(3) = timeMap(1:19)
    lpain(4) = 'PTEMPEI'
    lchin(4) = tempPrev(1:19)
    nbFieldIn = 4

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Generate new RESU_ELEM name
    call jedetr(matrElem(1:19)//'.RELR')
    call memare('V', matrElem, model, option)
    resuElem = matrElem(1:8)//'.ME000'

! - Set output field
    lpaout(1) = 'PMATTTR'
    lchout(1) = resuElem(1:19)

!
    chvite = '????'
    iconv = 0
    do iLoad = 1, nbLoad
        loadName = listLoadName(iLoad) (1:8)
        loadField = loadName(1:8)//'.CHTH.CONVE'
        call exisd('CHAMP_GD', loadField, iret)
        if (iret .ne. 0) then
            iconv = iconv+1
            if (iconv .gt. 1) then
                call utmess('F', 'CHARGES8_5')
            end if

            call memare('V', matrElem, model, option)

! --------- Get speed field
            call jeveuo(loadField(1:19)//'.VALE', 'L', vk8=loadFieldVale)
            chvite = loadFieldVale(1)
            lpain(5) = 'PVITESR'
            lchin(5) = chvite(1:19)
            nbFieldIn = 5

! --------- Compute
            call codent(iLoad, 'D0', lchout(1) (12:14))
            call calcul('S', option, modelLigrel, &
                        6, lchin, lpain, &
                        1, lchout, lpaout, &
                        'V', 'OUI')
            call reajre(matrElem, lchout(1), 'V')
!
        end if
    end do
!
    call jedema()
!
end subroutine
