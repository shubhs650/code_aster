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
subroutine kitPrepBehaviour(compor, nvi_tot, comporFlua, comporPlas)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    integer(kind=8), intent(in) :: nvi_tot
    character(len=16), intent(out) :: comporFlua(COMPOR_SIZE)
    character(len=16), intent(out) :: comporPlas(COMPOR_SIZE)
!
! --------------------------------------------------------------------------------------------------
!
! KIT_DDI
!
! Prepare fields for behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  compor          : behaviour
! In  nvi_tot         : total number of internal variables
! Out comporFlua      : behaviour for creep
! Out comporPlas      : behaviour for plasticity
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nume_plas, nume_flua
    integer(kind=8) :: nvi_flua, nvi_plas
!
! --------------------------------------------------------------------------------------------------
!
    comporFlua = 'VIDE'
    comporPlas = 'VIDE'
    read (compor(CREEP_NVAR), '(I16)') nvi_flua
    read (compor(PLAS_NVAR), '(I16)') nvi_plas
    read (compor(PLAS_NUME), '(I16)') nume_plas
    read (compor(CREEP_NUME), '(I16)') nume_flua
    ASSERT(nvi_tot .eq. (nvi_flua+nvi_plas))

! - Prepare COMPOR <CARTE> for creeping
    comporFlua(RELA_NAME) = compor(CREEP_NAME)
    write (comporFlua(NVAR), '(I16)') nvi_flua
    comporFlua(DEFO) = compor(DEFO)
    write (comporFlua(NUME), '(I16)') nume_flua

! - Prepare COMPOR <CARTE> for plasticity
    comporPlas(RELA_NAME) = compor(PLAS_NAME)
    write (comporPlas(NVAR), '(I16)') nvi_plas
    comporPlas(DEFO) = compor(DEFO)
    write (comporPlas(NUME), '(I16)') nume_plas
!
end subroutine
