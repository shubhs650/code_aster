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
subroutine varinonu_wrap(ligrel, compor, list_elem, list_vari, list_cmp)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/varinonu.h"
    character(len=*) :: ligrel, compor, list_elem, list_vari, list_cmp
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_elem, nb_cmp
    integer(kind=8), pointer :: v_list_elem(:) => null()
    character(len=16), pointer :: v_list_vari(:) => null()
    character(len=8), pointer :: v_list_cmp(:) => null()
!
    call jeveuo(list_elem, 'L', vi=v_list_elem)
    call jelira(list_elem, 'LONMAX', nb_elem)
    call jeveuo(list_vari, 'L', vk16=v_list_vari)
    call jelira(list_vari, 'LONMAX', nb_cmp)
    call jeveuo(list_cmp, 'E', vk8=v_list_cmp)
!
    call varinonu(ligrel, compor, nb_elem, v_list_elem, nb_cmp, v_list_vari, v_list_cmp)
!
end subroutine
