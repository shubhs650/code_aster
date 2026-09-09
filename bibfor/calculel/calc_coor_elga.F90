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
subroutine calc_coor_elga(modelZ, ligrel, chgeom, chgaus, &
                          caraElemZ_)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/xajcin.h"
!
    character(len=*), intent(in) :: modelZ
    character(len=19), intent(in) :: ligrel
    character(len=19), intent(in) :: chgeom, chgaus
    character(len=*), optional, intent(in) :: caraElemZ_
!
! --------------------------------------------------------------------------------------------------
!
! Compute <CARTE> with informations on Gauss points
!
! --------------------------------------------------------------------------------------------------
!
! In  model      : model
! In  ligrel     : list of elements where computing
! In  chgeom     : name of <CARTE> for geometry
! In  chgaus     : name of <CARTE> with informations on Gauss points
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    character(len=16), parameter :: option = 'COOR_ELGA'
    integer(kind=8) :: nfiss
!
! --------------------------------------------------------------------------------------------------
!

! - Initializations
    call dismoi('NB_FISS_XFEM', modelZ, 'MODELE', repi=nfiss)
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    nbFieldIn = 1

! - Add fields for orientation
    if (present(caraElemZ_)) then
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElemZ_)
    end if

! - Add fields for XFEM
    if (nfiss .gt. 0) then
        call xajcin(modelZ, option, nbFieldInMax, lchin, lpain, nbFieldIn)
    end if

! - Add output field and compute
    lpaout(1) = 'PCOORPG'
    lchout(1) = chgaus
    call calcul('S', option, ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')
!
end subroutine
