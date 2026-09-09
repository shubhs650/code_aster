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
subroutine te0093(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/tefrep.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: MECA_AXIS, MECA_AXIS_SI
!           MECA_AXIS_INCO_UPO, MECA_AXIS_INCO_UPG, MECA_AXIS_INCO_UP
!           MECA_AXIS_GRAD_INCO, MECA_AXIS_GRAD_VARI
!           MECA_C_PLAN, MECA_C_PLAN_SI
!           MECA_D_PLAN, MECA_D_PLAN_SI
!           MECA_D_PLAN_INCO_UPO, MECA_D_PLAN_INCO_UPG, MECA_D_PLAN_INCO_UP
!           MECA_D_PLAN_GRAD_INCO, MECA_D_PLAN_GRAD_VARI
!           MECA_D_PLAN_DIL (FORMULATIONS DIL AND DIL_INCO)
!
! Options: CHAR_MECA_FR2D2D
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ndimSpace = 2
    integer(kind=8) :: jvWeight, jvShape, jvDShape
    integer(kind=8) :: jvGeom, jvForc, jvVect
    integer(kind=8) :: nno, npg
    integer(kind=8) :: kpg, iNode, iDof
    integer(kind=8) :: jdec, kdec
    real(kind=8) :: jacWeight, r, fx, fy
    aster_logical :: lAxis
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, &
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShape)
    lAxis = lteatt('AXIS', 'OUI')

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call tefrep(option, 'PFR2D2D', jvForc)

! - Get output fields
    call jevech('PVECTUR', 'E', jvVect)
    do iDof = 1, ndimSpace*nno
        zr(jvVect+iDof-1) = 0.d0
    end do

! - Loop on Gauss points
    do kpg = 1, npg
        kdec = (kpg-1)*nno

! ----- Compute jacobian
        call dfdm2d(nno, kpg, jvWeight, jvDShape, zr(jvGeom), jacWeight)
        if (lAxis) then
            r = 0.d0
            do iNode = 1, nno
                r = r+zr(jvGeom+ndimSpace*(iNode-1))*zr(jvShape+kdec+iNode-1)
            end do
            jacWeight = jacWeight*r
        end if

! ----- Compute force at Gauss point from node value
        fx = 0.d0
        fy = 0.d0
        do iNode = 1, nno
            jdec = (iNode-1)*ndimSpace
            fx = fx+zr(jvShape+kdec+iNode-1)*zr(jvForc+jdec)
            fy = fy+zr(jvShape+kdec+iNode-1)*zr(jvForc+jdec+1)
        end do

! ----- Compute force
        do iNode = 1, nno
            zr(jvVect+ndimSpace*(iNode-1)) = zr(jvVect+ndimSpace*(iNode-1))+ &
                                             jacWeight*fx*zr(jvShape+kdec+iNode-1)
            zr(jvVect+ndimSpace*(iNode-1)+1) = zr(jvVect+ndimSpace*(iNode-1)+1)+ &
                                               jacWeight*fy*zr(jvShape+kdec+iNode-1)
        end do
    end do
!
end subroutine
