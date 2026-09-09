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

subroutine te0507(nomopt, nomte)
!
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_matrix_module
    use HHO_size_module
    use HHO_statcond_module
    use HHO_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/tecach.h"
#include "asterfort/readVector.h"
#include "asterfort/writeVector.h"
!
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - Static condensation
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16), intent(in) :: nomte, nomopt
!
! --- Local variables
!
    integer(kind=8) :: cbs, fbs, total_dofs, iad, iret
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
!
    real(kind=8), dimension(MSIZE_TDOFS_VEC) :: rhs_elem, rhs_cond, rhs_decond
    type(HHO_matrix) :: lhs_elem, lhs_cond, lhs_decond
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
    call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
    call lhs_elem%initialize(total_dofs, total_dofs)
    call lhs_elem%read("PMAELS1", ASTER_TRUE)
!
    call tecach("OON", "PVEELE1", "L", iret, iad=iad)
    if (iret == 0) then
        call readVector("PVEELE1", total_dofs, rhs_elem)
    else
        rhs_elem = 0.d0
    end if
!
    call hhoCondStatic(cbs, lhs_elem, rhs_elem, &
                       lhs_cond, rhs_cond, lhs_decond, rhs_decond)
!
    call lhs_cond%write("PMATUUR", ASTER_TRUE)
    call writeVector("PVECTUR", total_dofs, rhs_cond)
!
    call lhs_decond%write("PMATUND", ASTER_FALSE)
    call writeVector("PVECTUD", total_dofs, rhs_decond)
!
    call lhs_elem%free()
    call lhs_cond%free()
    call lhs_decond%free()
!
end subroutine
