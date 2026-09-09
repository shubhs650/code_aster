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

subroutine spmats(imater, matint, matpou)
!
! --------------------------------------------------------------------------------------------------
!
!  Retourne les noms des materiaux des éléments 3D_INTERF_POU
!
! --------------------------------------------------------------------------------------------------
!   out
!       matint     : nom du materiau de l'interface
!       matpou     : nom du materiau de la poutre
! --------------------------------------------------------------------------------------------------
!
!
    use Behaviour_type

    implicit none

#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/rcmats.h"
#include "asterfort/rccome.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), intent(in) :: imater
    character(len=*), intent(out) :: matint, matpou
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: nomats(2)
    integer(kind=8) :: nbmats, kmat, icodre
    character(len=8) :: nomi
!
! --------------------------------------------------------------------------------------------------

! - Get all material names
    call rcmats(imater, nbmats, nomats)
    ASSERT(nbmats .eq. 2)

! - Sort materials using their phenomenon name
    do kmat = 1, nbmats
        nomi = nomats(kmat)
        call rccome(nomi, 'INTERF_POU_', icodre)
        if (icodre .eq. 0) then
            matint = nomi
        else
            matpou = nomi
        end if
    end do

end subroutine
