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
subroutine te0550(option, nomte)
!
! --------------------------------------------------------------------------------------------------
!     CALCUL DES FORCES ELEMENTAIRES LINEIQUES POUR LES ELEMENTS BARRE
! --------------------------------------------------------------------------------------------------
!
! option : nom de l'option à calculer
!       CHAR_MECA_PESA_R    : charges de pesanteur
! --------------------------------------------------------------------------------------------------
!
    implicit none
    character(len=*) :: option, nomte
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/teattr.h"
#include "asterfort/te0550_implement.h"
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
! --------------------------------------------------------------------------------------------------
    character(len=8) :: attrib
    integer(kind=8):: nno, npg, ndim_sp
    integer(kind=8):: jv_poids, jv_vff, jv_dxi_ff, jv_geom, jv_sect, jv_materc, jv_pesa, jv_fext
! --------------------------------------------------------------------------------------------------

! - Get parameters of element
    call elrefe_info(fami=fami, nno=nno, npg=npg, jpoids=jv_poids, jvf=jv_vff, jdfde=jv_dxi_ff)
    call teattr('S', 'DIM_COOR_MODELI', attrib)
    read (attrib, '(I8)') ndim_sp

! - Get arguments of the option
    call jevech('PCAGNBA', 'L', jv_sect)
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PPESANR', 'L', jv_pesa)
    call jevech('PMATERC', 'L', jv_materc)
    call jevech('PVECTUR', 'E', jv_fext)

! - Option computation
    call te0550_implement(option, fami, nno, npg, ndim_sp, &
                          zr(jv_poids), zr(jv_vff), zr(jv_dxi_ff), zr(jv_sect), zr(jv_geom), &
                          zr(jv_pesa), zi(jv_materc), zr(jv_fext))

end subroutine
