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
subroutine w190ca(model, caraElem, chmar1, chefge, chamfer, chefge0, chmar2)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/calcul.h"
#include "asterfort/exlim3.h"
#include "asterfort/setStructFields.h"
!
    character(len=8), intent(in) :: model, caraElem
    character(len=19), intent(in) :: chmar1, chmar2, chefge, chamfer, chefge0
!
! --------------------------------------------------------------------------------------------------
!
!  VERI_FERRAILLAGE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    character(len=16), parameter:: option = 'MARG_ELEM'
    character(len=19) :: ligrel
!
! --------------------------------------------------------------------------------------------------
!
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

    call exlim3('AFFE', 'G', model, ligrel)

! - Set input fields
    nbFieldIn = 1
    lpain(nbFieldIn) = 'PVFER1'
    lchin(nbFieldIn) = chmar1
    nbFieldIn = nbFieldIn+1
    lpain(nbFieldIn) = 'PEFFORR'
    lchin(nbFieldIn) = chefge
    nbFieldIn = nbFieldIn+1
    lpain(nbFieldIn) = 'PEFFOR0'
    lchin(nbFieldIn) = chefge0
    nbFieldIn = nbFieldIn+1
    lpain(nbFieldIn) = 'PVFER0'
    lchin(nbFieldIn) = chamfer

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Set output field
    lpaout(1) = 'PVFER2'
    lchout(1) = chmar2

! - Compute
    call calcul('S', option, ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'G', 'OUI')
!
end subroutine
