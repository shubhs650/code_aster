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
subroutine tumgamma(nbNode, nbFourier, nbDof)
!
    use pipeElem_type
    use pipeElem_module
    use beamElem_module
!
    implicit none
!
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/beamElem_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/mavec.h"
#include "asterfort/moytem.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/pmavec.h"
#include "asterfort/promat.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nbNode, nbFourier, nbDof
!
! --------------------------------------------------------------------------------------------------
!
! Compute inertial vector for pipe element
!
! In  nbNode           : number of nodes in element
! In  nbFourier        : number of Fourier modes
! In  nbDof            : number of DOF in element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "MASS"
    integer(kind=8) :: jvf, jdfde, jdfd2, jcoopg, jpoids
    integer(kind=8) :: jvAcce
    integer(kind=8) :: jvVect
    real(kind=8) :: mass(nbDof, nbDof)
    real(kind=8) :: radiusLayer
    real(kind=8) :: poids, weightLayer(2*PIPE_MAX_LAYERS+1), weightSect(2*PIPE_MAX_SECTORS+1)
    integer(kind=8) :: jvMaterc
    real(kind=8) :: meanTemp, rho
    real(kind=8) :: jacobi, xpg(PIPE_MAX_NPG)
    real(kind=8) :: phi
    integer(kind=8) :: nbLayer, nbSect
    integer(kind=8) :: nspgLayer, nspgSect, npg, nspg
    integer(kind=8) :: kspgLayer, kspgSect, kpg
    integer(kind=8) :: iDof, jDof, iret, shiftDOF
    real(kind=8) :: nvec(PIPE_NBDOF_BEAM, nbDof), tnvec(nbDof, PIPE_NBDOF_BEAM)
    type(pipeElem_Prop) :: pipeElem
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg, &
                     jpoids=jpoids, jcoopg=jcoopg, jvf=jvf, jdfde=jdfde, jdfd2=jdfd2)
    ASSERT(npg .le. PIPE_MAX_NPG)
    ASSERT(PIPE_NBDOF_BEAM .eq. BEAM_NBDOF)
    shiftDOF = (6+3+6*(nbFourier-1))

! - Output field
    call jevech('PVECTUR', 'E', jvVect)

! - Get parameters about layers and sections
    call pipeGetSubPoints(nbLayer, nbSect, &
                          nspgLayer, nspgSect, nspg, &
                          weightLayer, weightSect)

! - Get properties of pipe
    call pipeGetProperties(nbNode, nbFourier, pipeElem)

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

! - Get density
    call jevech('PMATERC', 'L', jvMaterc)
    call pipeGetDensity(jvMaterc, rho, meanTemp)

! - Loop on Gauss points (on segment)
    nvec = 0.d0
    tnvec = 0.d0
    mass = 0.d0
    do kpg = 1, npg
! ----- Loop on sub-points in thickness
        do kspgLayer = 1, nspgLayer
! --------- Radius of layer
            call pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)

! --------- Loop on sub-points in section
            do kspgSect = 1, nspgSect
! ------------- Coordinate of sub-point in section
                phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)

! ------------- Get N matrix of shape functions (beam part)
                call beamNMatr(pipeElem%beamElem, nbDof, npg, &
                               kpg, xpg, &
                               phi, radiusLayer, shiftDOF, &
                               zr(jvf), nvec, tnvec)

! ------------- Get N matrix of shape functions (shell part)
                call pipeNMatr(pipeElem, nbDof, kpg, &
                               phi, &
                               zr(jvf), nvec, tnvec)

! ------------- Jacobian
                call pipeGetJacobian(pipeElem, nbSect, nbLayer, phi, radiusLayer, &
                                     jacobi)

! ------------- Weighted jacobian of current sub-point
                poids = zr(jpoids-1+kpg)* &
                        weightLayer(kspgLayer)*weightSect(kspgSect)* &
                        jacobi*rho

! ------------- Compute mass matrix
                call promat(tnvec, nbDof, nbDof, 6, nvec, &
                            6, 6, nbDof, mass)
                do iDof = 1, nbDof
                    do jDof = 1, iDof
                        mass(iDof, jDof) = mass(iDof, jDof)+poids*mass(iDof, jDof)
                        mass(jDof, iDof) = mass(iDof, jDof)
                    end do
                end do
            end do
        end do
    end do

! - Change base of mass matrix
    call pipeBaseForMatr(pipeElem, nbDof, mass)

! - Get acceleration
    call jevech('PACCELR', 'L', jvAcce)

! - Product and storage
    call pmavec('ZERO', nbDof, mass, zr(jvAcce), zr(jvVect))
!
end subroutine
