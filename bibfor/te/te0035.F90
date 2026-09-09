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
subroutine te0035(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dxbsig.h"
#include "asterfort/dxefgi_fonc.h"
#include "asterfort/dxefgi.h"
#include "asterfort/dxefgt.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT/DKTG/DST/Q4G/Q4GG
!
! Options: CHAR_MECA_TEMP_R
!          CHAR_MECA_EPSI_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ncomp = 6
    integer(kind=8) :: nno, npg, ivf
    integer(kind=8) :: i, jvGeom, jvVect, jvEpsi, kpg
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: epsiR(32)
    real(kind=8) :: bEfge(24), efge(32)
    character(len=8) :: epsiF(6)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno, npg=npg, jvf=ivf)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of geometry
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! - CALCUL DES EFFORTS GENERALISES D'ORIGINE THERMIQUE AUX POINTS D'INTEGRATION
    if (option .eq. 'CHAR_MECA_TEMP_R') then
        call dxefgt(plateCara, plateOrie, efge)
    else if (option .eq. 'CHAR_MECA_EPSI_R') then
        call jevech('PEPSINR', 'L', jvEpsi)
        do kpg = 1, npg
            epsiR(ncomp*(kpg-1)+1) = zr(jvEpsi+ncomp*(kpg-1)+1-1)
            epsiR(ncomp*(kpg-1)+2) = zr(jvEpsi+ncomp*(kpg-1)+2-1)
            epsiR(ncomp*(kpg-1)+3) = zr(jvEpsi+ncomp*(kpg-1)+3-1)
            epsiR(ncomp*(kpg-1)+4) = zr(jvEpsi+ncomp*(kpg-1)+4-1)
            epsiR(ncomp*(kpg-1)+5) = zr(jvEpsi+ncomp*(kpg-1)+5-1)
            epsiR(ncomp*(kpg-1)+6) = zr(jvEpsi+ncomp*(kpg-1)+6-1)
        end do
        call dxefgi(plateCara, plateOrie, &
                    npg, epsiR, efge)
    else if (option .eq. 'CHAR_MECA_EPSI_F') then
        call jevech('PEPSINF', 'L', jvEpsi)
        epsiF(1) = zk8(jvEpsi+1-1)
        epsiF(2) = zk8(jvEpsi+2-1)
        epsiF(3) = zk8(jvEpsi+3-1)
        epsiF(4) = zk8(jvEpsi+4-1)
        epsiF(5) = zk8(jvEpsi+5-1)
        epsiF(6) = zk8(jvEpsi+6-1)
        call dxefgi_fonc(plateCara, plateOrie, &
                         nno, npg, &
                         epsiF, zr(jvGeom), zr(ivf), efge)
    else
        ASSERT(ASTER_FALSE)
    end if

! - CALCUL DES EFFORTS INTERNES D'ORIGINE THERMIQUE
    call dxbsig(plateCara, plateOrie, &
                nomte, 'FORC_NODA', &
                xyzl, pgl, efge, &
                bEfge)

! - AFFECTATION DU VECTEUR DES FORCES ELEMENTAIRES
    call jevech('PVECTUR', 'E', jvVect)
    do i = 1, nno*6
        zr(jvVect+i-1) = bEfge(i)
    end do
!
end subroutine
