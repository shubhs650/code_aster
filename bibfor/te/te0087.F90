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
subroutine te0087(option, nomte)
!
    use BehaviourStrain_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epsvmc.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: C_PLAN, D_PLAN, AXIS
!
! Options: EPME_ELGA, EPMG_ELGA, EPSG_ELGA, EPSL_ELGA, EPSI_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8) :: nharm
    integer(kind=8) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    integer(kind=8) :: ndim, nno, npg, nbsig, nbEpsi
    integer(kind=8) :: kpg, iret, iEpsi
    integer(kind=8) :: jvGeom, jvDisp, jvTime, jvEpsi, jvMaterc
    real(kind=8) :: epsi(100), time
    integer(kind=8) :: strainType
    aster_logical :: lStrainMeca, lStrainVarc
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)
!
    nbsig = nbsigm()
    nbEpsi = nbsig
    ASSERT(nbEpsi .le. 6)
    ASSERT(npg .le. 25)
    epsi = 0.d0

! - Get type of strain from option
    call getStrainType(option, strainType, lStrainMeca)

! - When compute strains from external state variables ?
    lStrainVarc = lStrainMeca .or. lteatt('C_PLAN', 'OUI')

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Current displacements (nodes)
    call jevech('PDEPLAR', 'L', jvDisp)

! - Get current time
    time = r8vide()
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jvTime)
    if (jvTime .ne. 0) then
        time = zr(jvTime)
    end if

! - Get material parametes
    if (lStrainVarc) then
        call jevech('PMATERC', 'L', jvMaterc)

! ----- Initializations of material parameters on current cell
        call initParaCell(fami, zi(jvMaterc), materPara)

! ----- Set local coordinate system from user
        call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)
    end if

! - Compute mechanical strains or total strains
    call epsvmc(nno, ndim, nbEpsi, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                zr(jvGeom), zr(jvDisp), &
                time, nharm, &
                strainType, lStrainMeca, materPara, &
                epsi)

! - Save strains
    call jevech('PDEFOPG', 'E', jvEpsi)
    do kpg = 1, npg
        do iEpsi = 1, nbEpsi
            zr(jvEpsi+nbEpsi*(kpg-1)+iEpsi-1) = epsi(nbEpsi*(kpg-1)+iEpsi)
        end do
    end do
!
end subroutine
