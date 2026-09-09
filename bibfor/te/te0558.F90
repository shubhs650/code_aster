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
subroutine te0558(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/nmbamb.h"
#include "asterfort/ngforc.h"
#include "asterfort/teattr.h"
#include "asterfort/elrefe_info.h"
    character(len=16) :: option, nomte
!-----------------------------------------------------------------------
! FORC_NODA: Nodal forces
! REFE_FORC_NODA: to be implemented
!
! Elements: MECA_BARRE and MECA_BARRE_2D
! ----------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
! --------------------------------------------------------------------------------------------------
    character(len=8):: attrib
    integer(kind=8) :: nno, npg, ndim_sp, ndim_fe, nddl, neps
    integer(kind=8) :: jv_geom, jv_cont, jv_vectu
    integer(kind=8) :: jv_poids, jv_dfde
    real(kind=8):: aire
    real(kind=8), allocatable:: b(:, :, :), w(:, :), ni2ldc(:, :)
! --------------------------------------------------------------------------------------------------

    ! Get parameters of element
    call elrefe_info(fami=fami, ndim=ndim_fe, nno=nno, npg=npg, jpoids=jv_poids, jdfde=jv_dfde)
    call teattr('S', 'DIM_COOR_MODELI', attrib)
    read (attrib, '(I8)') ndim_sp

    ! Parametres de l'option
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PSIEFR', 'L', jv_cont)
    call jevech('PVECTUR', 'E', jv_vectu)

!   Cross section area is useless since Bt * N*(1/S) * (w*S) = Bt * (N*1) * (w*1)
    aire = 1.d0

    ! Kinematics matrix B
    call nmbamb(ndim_sp, nno, npg, zr(jv_geom), aire, &
                zr(jv_dfde), zr(jv_poids), nddl, neps, b, w, ni2ldc)

    ! Force nodale
    call ngforc(w, b, ni2ldc, zr(jv_cont), zr(jv_vectu))

    deallocate (b, w, ni2ldc)

end subroutine
