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

subroutine te0574(nomopt, nomte)
!
    use HHO_init_module, only: hhoInfoInitFace
    use HHO_matrix_module
    use HHO_size_module
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
    integer(kind=8) :: fbs, iad, iret
    type(HHO_Data) :: hhoData
    type(HHO_Face) :: hhoFace
!
    real(kind=8), dimension(MSIZE_FDOFS_VEC) :: rhs_elem
    type(HHO_matrix) :: lhs_elem
!
    call hhoInfoInitFace(hhoFace, hhoData)
!
    call hhoMecaFaceDofs(hhoFace, hhoData, fbs)
!
    call lhs_elem%initialize(fbs, fbs, 0.d0)
    call tecach("OON", "PMAELS1", "L", iret, iad=iad)
    if (iret == 0) then
        call lhs_elem%read("PMAELS1", ASTER_TRUE)
    end if
!
    call tecach("OON", "PVEELE1", "L", iret, iad=iad)
    if (iret == 0) then
        call readVector("PVEELE1", fbs, rhs_elem)
    else
        rhs_elem = 0.d0
    end if
!
    call lhs_elem%write("PMATUUR", ASTER_TRUE)
    call writeVector("PVECTUR", fbs, rhs_elem)
!
    call lhs_elem%free()
!
end subroutine
