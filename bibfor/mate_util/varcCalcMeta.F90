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
subroutine varcCalcMeta(modelZ, &
                        nbFieldIn, nbFieldOut, &
                        lpain, lchin, &
                        lpaout, lchout, &
                        vectElemZ)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reajre.h"
!
    character(len=*), intent(in) :: modelZ
    integer(kind=8), intent(in) :: nbFieldIn, nbFieldOut
    character(len=8), intent(in) :: lpain(*), lpaout(*)
    character(len=19), intent(in) :: lchin(*)
    character(len=19), intent(inout) :: lchout(*)
    character(len=*), intent(in) :: vectElemZ
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! Call to CALCUL special for metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  nbFieldIn        : effective number of input fields
! In  nbFieldOut       : effective number of output fields
! In  lpain            : list of input parameters
! In  lchin            : list of input fields
! In  lpaout           : list of output parameters
! IO  lchout           : list of output fields
! In  vectElem         : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'CHAR_MECA_META_Z'
    character(len=1), parameter :: jvBase = "V"
    character(len=8) :: newnom
    character(len=19) :: resuElem, modelLigrel
    integer(kind=8) :: nbResuElem
    character(len=24), pointer :: relr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)

! - Get last resuElem
    call jelira(vectElemZ(1:8)//'           .RELR', 'LONUTI', nbResuElem)
    call jeveuo(vectElemZ(1:8)//'           .RELR', 'L', vk24=relr)
    resuElem = relr(nbResuElem) (1:19)
    newnom = resuElem(10:16)

! - Generate new resuElem
    call gcnco2(newnom)
    resuElem(10:16) = newnom(2:8)
    call corich('E', resuElem, ichin_=-1)
    lchout(1) = resuElem

! - Compute
    call calcul('C', option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')
    call reajre(vectElemZ, resuElem, jvBase)
!
end subroutine
