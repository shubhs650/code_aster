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
! aslint: disable=W1306
!
subroutine turigi(nbNode, nbFourier, nbDof)
!
    use pipeElem_type
    use pipeElem_module
    use beamElem_module
!
    implicit none
!
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/kcoude.h"
#include "asterfort/mavec.h"
#include "asterfort/moytem.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nbNode, nbFourier, nbDof
!
! --------------------------------------------------------------------------------------------------
!
! Compute rigidity matrix for pipe element
!
! In  nbNode           : number of nodes in element
! In  nbFourier        : number of Fourier modes
! In  nbDof            : number of DOF in element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: jvf, jdfde, jdfd2, jcoopg, jpoids
    integer(kind=8) :: jvMatr
    real(kind=8) :: rigi(nbDof, nbDof)
    real(kind=8) :: radiusLayer
    real(kind=8) :: poids, weightLayer(2*PIPE_MAX_LAYERS+1), weightSect(2*PIPE_MAX_SECTORS+1)
    integer(kind=8) :: jvMaterc
    real(kind=8) :: meanTemp
    real(kind=8) :: elasMatr(PIPE_TENS_SIZE, PIPE_TENS_SIZE)
    real(kind=8) :: jacobi, xpg(PIPE_MAX_NPG)
    real(kind=8) :: phi, zeta
    integer(kind=8) :: nbLayer, nbSect
    integer(kind=8) :: nspgLayer, nspgSect, npg, nspg
    integer(kind=8) :: kspgLayer, kspgSect, kpg
    integer(kind=8) :: iDof, jDof, iret, nc
    real(kind=8) :: b(PIPE_TENS_SIZE, nbDof)
    type(pipeElem_Prop) :: pipeElem
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg, &
                     jpoids=jpoids, jcoopg=jcoopg, jvf=jvf, jdfde=jdfde, jdfd2=jdfd2)
    ASSERT(npg .le. PIPE_MAX_NPG)

! - Output field
    call jevech('PMATUUR', 'E', jvMatr)

! - Get parameters about layers and sections
    call pipeGetSubPoints(nbLayer, nbSect, &
                          nspgLayer, nspgSect, nspg, &
                          weightLayer, weightSect)

! - Get properties of pipe
    call pipeGetProperties(nbNode, nbFourier, pipeElem)
    if (pipeElem%lModiMetric) then
        if (pipeElem%beamElem%sectPipe%thickness/pipeElem%beamElem%sectPipe%radiusMoy .gt. 0.25d0) &
            then
            call utmess('A', 'PIPE1_54')
        end if
    end if

! - Get coordinates of Gauss points (on segment)
    do kpg = 1, npg
        xpg(kpg) = zr(jcoopg-1+kpg)
    end do

! - Compute mean temperature (on all point and "sous-point" gauss)
    iret = 0
    call moytem('RIGI', npg, nspg, '+', meanTemp, iret)
    if (iret .ne. 0) then
        meanTemp = 0.d0
    end if

! - Get elastic properties
    call jevech('PMATERC', 'L', jvMaterc)
    call pipeGetElasProp(zi(jvMaterc), temp_=meanTemp, &
                         c_=elasMatr)

! - Loop on Gauss points (on segment)
    rigi = 0.d0
    do kpg = 1, npg
! ----- Loop on sub-points in thickness
        do kspgLayer = 1, nspgLayer
! --------- Radius of layer
            call pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)

! --------- Coordinate of sub-point in thickness
            zeta = (kspgLayer-1)*pipeElem%beamElem%sectPipe%thickness/(2.d0*nbLayer)- &
                   pipeElem%beamElem%sectPipe%thickness/2.d0

! --------- Loop on sub-points in section
            do kspgSect = 1, nspgSect
! ------------- Coordinate of sub-point in section
                phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)

! ------------- Compute B matrix (u => epsi)
                call pipeBMatr(pipeElem, nbDof, &
                               npg, xpg, &
                               phi, zeta, radiusLayer, &
                               kpg, zr(jvf), zr(jdfde), zr(jdfd2), &
                               b)

! ------------- Jacobian
                call pipeGetJacobian(pipeElem, nbSect, nbLayer, phi, radiusLayer, &
                                     jacobi)

! ------------- Weighted jacobian of current sub-point
                poids = zr(jpoids-1+kpg)* &
                        weightLayer(kspgLayer)*weightSect(kspgSect)* &
                        jacobi

! ------------- Product by elasticity matrix
                call kcoude(nbDof, poids, b, elasMatr, rigi)
            end do
        end do
    end do
!
    do iDof = 1, nbDof
        do jDof = 1, nbDof
            rigi(iDof, jDof) = rigi(iDof, jDof)
        end do
    end do

! - Change base of rigidity matrix
    call pipeBaseForMatr(pipeElem, nbDof, rigi)

! - Convert matrix to storage
    nc = nbDof*(nbDof+1)/2
    call mavec(rigi, nbDof, zr(jvMatr), nc)
!
end subroutine
