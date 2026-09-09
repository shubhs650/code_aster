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

subroutine te0549(option, nomte)

    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/ngmatr.h"
#include "asterfort/nmbamb.h"
#include "asterfort/rcvalb.h"
#include "asterfort/teattr.h"
#include "jeveux.h"

    character(len=16), intent(in) :: option, nomte

! --------------------------------------------------------------------------------------------------
! Elements: BARRE / 2D_BARRE
! Options: RIGI_MECA
! --------------------------------------------------------------------------------------------------
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: ndim_sp, nno, npg, nddl, neps, g
    integer(kind=8):: jv_poids, jv_dfde, jv_materc, jv_geom, jv_sect, jv_matuu
    integer(kind=8) :: iok(1)
    real(kind=8) :: aire
    character(len=8) :: attrib
    real(kind=8), allocatable :: b(:, :, :), w(:, :), ni2ldc(:, :), dsidep(:, :, :)
! --------------------------------------------------------------------------------------------------

    ! Get parameters of element
    call elrefe_info(fami=fami, nno=nno, npg=npg, jpoids=jv_poids, jdfde=jv_dfde)
    call teattr('S', 'DIM_COOR_MODELI', attrib)
    read (attrib, '(I8)') ndim_sp

    ! Get option arguments
    call jevech('PCAGNBA', 'L', jv_sect)
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_materc)
    call jevech('PMATUUR', 'E', jv_matuu)

    ! Cross section area
    aire = zr(jv_sect)

    ! Kinematics
    call nmbamb(ndim_sp, nno, npg, zr(jv_geom), aire, &
                zr(jv_dfde), zr(jv_poids), nddl, neps, b, w, ni2ldc)

    ! Young modulus and local matrices
    allocate (dsidep(1, 1, npg))
    do g = 1, npg
        call rcvalb(fami, g, 1, '+', zi(jv_materc), ' ', 'ELAS', 0, ' ', [0.d0], &
                    1, 'E', dsidep(1, 1, g), iok, 1)
    end do

    ! Compute matrix for the element
    call ngmatr(nddl, neps, npg, w, b, dsidep, matsym=ASTER_TRUE, matuu=zr(jv_matuu))

    ! Memory management
    deallocate (b, w, ni2ldc, dsidep)

end subroutine
