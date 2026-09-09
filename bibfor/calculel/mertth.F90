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
subroutine mertth(model, loadNameJv, loadInfoJv, &
                  caraElem, materCode, &
                  timeMapMatr, timeMapMove, &
                  tempPrev, tempIter, &
                  matrElem)
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
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/load_list_info.h"
#include "asterfort/megeom.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/setStructFields.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: loadNameJv, loadInfoJv
    character(len=8), intent(in) :: caraElem
    character(len=24), intent(in) :: materCode
    character(len=24), intent(in) :: timeMapMatr, timeMapMove
    character(len=24), intent(in) :: tempPrev, tempIter
    character(len=19), intent(inout) :: matrElem
!
! --------------------------------------------------------------------------------------------------
!
! Thermic - Matrix
!
! Elementary matrix for transport (volumic and surfacic terms)
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  caraElem         : name of elementary characteristics (field)
! In  loadNameJv       : name of object for list of loads name
! In  loadInfoJv       : name of object for list of loads info
! In  timeMapMatr      : time (<CARTE>)
! In  timeMapMove      : modified time (<CARTE>) for THER_NON_LINE_MO
! In  tempPrev         : previous temperature
! In  tempIter         : temperature field at current Newton iteration
! IO  matrElem         : name of matrElem result
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: l_stat = ASTER_TRUE, lMove = ASTER_TRUE
    character(len=1), parameter :: jvBase = "V"
    real(kind=8) :: theta
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    character(len=16), parameter :: option = 'RIGI_THER_TRANS'
    character(len=24) :: modelLigrel, loadLigrel
    character(len=24) :: chgeom, resuElem
    integer(kind=8) :: iret, nbFieldIn
    integer(kind=8) :: nbLoad, iLoad, loadNume
    aster_logical :: noLoadInList
    character(len=13) :: loadPreObject
    character(len=8) :: loadName
    character(len=24), pointer :: listLoadName(:) => null()
    integer(kind=8), pointer :: listLoadInfo(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

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

! - Allocate result
    call jeexin(matrElem(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        matrElem = '&&METRIG'
        call memare('V', matrElem, model, 'RIGI_THER')
    else
        call jedetr(matrElem(1:19)//'.RELR')
    end if

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = materCode(1:19)
    lpain(3) = 'PTEMPER'
    lchin(3) = tempPrev(1:19)
    lpain(4) = 'PTEMPEI'
    lchin(4) = tempIter(1:19)
    nbFieldIn = 4

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Generate new RESU_ELEM name
    resuElem = matrElem(1:8)//'.ME001'

! - Set output field
    lpaout(1) = 'PMATTTR'
    lchout(1) = resuElem(1:19)

! - Compute "volumic" term
    call calcul('S', option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')
    call reajre(matrElem, lchout(1), jvBase)

! - Add load terms
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "
    do iLoad = 1, nbLoad
        loadName = listLoadName(iLoad) (1:8)
        loadNume = listLoadInfo(nbLoad+iLoad+1)
        loadPreObject = loadName(1:8)//'.CHTH'
        loadLigrel = loadPreObject(1:13)//'.LIGRE'

! ----- Standard input fields
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom(1:19)
        lpain(2) = 'PTEMPEI'
        lchin(2) = tempIter(1:19)
        lpain(3) = 'PDEPLAR'
        lchin(3) = '&&DEPPLU'
        nbFieldIn = 3

! ----- Set output field
        lpaout(1) = 'PMATTTR'

        if (loadNume .gt. 0) then
            call compLoadMatr(l_stat, theta, &
                              model, timeMapMatr, &
                              loadNume, &
                              loadPreObject, loadLigrel, &
                              nbFieldIn, lpain, lchin, &
                              jvBase, resuElem, matrElem, &
                              lMove, timeMapMove)

        end if
    end do

!
end subroutine
