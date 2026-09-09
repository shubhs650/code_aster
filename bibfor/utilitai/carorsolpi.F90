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

subroutine carorsolpi(nutyel, ntyele, IsIntsolpi, ino1, ino2)
!
!
! --------------------------------------------------------------------------------------------------
!
!                           DETECTE S'IL S'AGIT D'UN ELEMENT INTERF_POU
!
!   OUT
!       IsIntsolpi : TRUE s'il s'agit d'un élément INTERF_POU
!       ino1, ino2 : indices des noeuds de la poutre
!
! --------------------------------------------------------------------------------------------------
!
    use cara_elem_parameter_module
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/teattr.h"
!

    integer(kind=8) :: nutyel, ntyele(*), ino1, ino2
    logical :: IsIntsolpi
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jj, iret2
    character(len=16) :: nunomel
    character(len=8) :: typmod2
!
! --------------------------------------------------------------------------------------------------
!
    do jj = 1, ACE_NB_TYPE_ELEM
        if (nutyel .eq. ntyele(jj)) then
            call jenuno(jexnum('&CATA.TE.NOMTE', nutyel), nunomel)
            call teattr('C', 'TYPMOD2', typmod2, iret2, typel=nunomel)
            if (typmod2 .eq. 'INTSOLPI') then
                IsIntsolpi = .TRUE.
                ino1 = 25
                ino2 = 23
                exit
            end if
        end if
    end do
!
end subroutine
