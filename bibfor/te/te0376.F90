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
subroutine te0376(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nmspfn.h"
#include "asterfort/matrot.h"
#include "asterfort/spmats.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Element: 3D_INTERF_POU
!
! Options: FORC_NODA
!
! --------------------------------------------------------------------------------------------------
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ndim = 3
    character(len=4), parameter :: fami = "RIGI"
    integer(kind=8) :: nno, nno_s, nno_p, npg, nddl
    integer(kind=8) :: ipoids, ivf, ivfs, ivfp, icoopg
    integer(kind=8) :: iorie, igeom, imater
    integer(kind=8) :: icontm, ivect
    real(kind=8) :: pgl(3, 3)
    character(len=8) :: matint, matpou
!
! --------------------------------------------------------------------------------------------------
!
! - Get general pointer info on element
    call elrefe_info(fami=fami, nno=nno, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf)

    nno_s = 8
    nno_p = 2
    nddl = 3*nno_s+6*nno_p
    ivfs = ivf
    ivfp = icoopg

    ASSERT(nno .eq. 27)
    ASSERT(npg .eq. 2 .or. npg .eq. 3 .or. npg .eq. 4)

! - Get input fields
    call jevech('PCAORIE', 'L', iorie)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imater)
    call jevech('PSIEFR', 'L', icontm)

! - Get multiple materials
    call spmats(imater, matint, matpou)

! - Get output fields
    call jevech('PVECTUR', 'E', ivect)
!
! - Get orientation
!
    call matrot(zr(iorie), pgl)

! - Main calculation
    call nmspfn(ndim, nno, nddl, &
                nno_p, nno_s, npg, &
                zr(ipoids), zr(ivfs), zr(ivfp), &
                pgl, zr(igeom), 3, zi(imater), matpou, &
                zr(icontm), zr(ivect))
!
end subroutine
