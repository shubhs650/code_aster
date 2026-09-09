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
subroutine te0016(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tefrep.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: MECA_3D, MECA_3D_SI
!           MECA_3D_INCO_UPO, MECA_3D_INCO_UPG, MECA_INCO_3D_UP
!           MECA_3D_GRAD_INCO, MECA_3D_GRAD_VARI,
!           MECA_3D_DIL (FORMULATIONS DIL AND DIL_INCO)
!
! Options: CHAR_MECA_FR3D3D
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ndim = 3
    integer(kind=8) :: jvWeight, jvShape, jvDShape
    integer(kind=8) :: jvGeom, jvForc, jvVect
    integer(kind=8) :: nno, npg
    integer(kind=8) :: kpg, iNode, iDof
    integer(kind=8) :: kdec, ldec
    real(kind=8) :: jacWeight, fx, fy, fz
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, &
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShape)

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call tefrep(option, 'PFR3D3D', jvForc)

! - Get output fields
    call jevech('PVECTUR', 'E', jvVect)
    do iDof = 1, ndim*nno
        zr(jvVect+iDof-1) = 0.d0
    end do

! - Loop on Gauss points
    do kpg = 1, npg
        ldec = (kpg-1)*nno

! ----- Compute jacobian with Gauss weight
        call dfdm3d(nno, kpg, jvWeight, jvDShape, zr(jvGeom), jacWeight)

! ----- Compute force at Gauss point from node value
        fx = 0.d0
        fy = 0.d0
        fz = 0.d0
        do iNode = 1, nno
            kdec = ndim*(iNode-1)
            fx = fx+zr(jvShape-1+ldec+iNode)*zr(jvForc+kdec)
            fy = fy+zr(jvShape-1+ldec+iNode)*zr(jvForc+kdec+1)
            fz = fz+zr(jvShape-1+ldec+iNode)*zr(jvForc+kdec+2)
        end do
        do iNode = 1, nno
            kdec = ndim*(iNode-1)
            zr(jvVect+kdec) = zr(jvVect+kdec)+jacWeight*fx*zr(jvShape+ldec+iNode-1)
            zr(jvVect+kdec+1) = zr(jvVect+kdec+1)+jacWeight*fy*zr(jvShape+ldec+iNode-1)
            zr(jvVect+kdec+2) = zr(jvVect+kdec+2)+jacWeight*fz*zr(jvShape+ldec+iNode-1)
        end do
    end do
!
end subroutine
