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
subroutine tusief(nbNode, nbFourier, nbDof)
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
#include "asterfort/moytem.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/prmave.h"
#include "asterfort/promat.h"
#include "asterfort/verifg.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nbNode, nbFourier, nbDof
!
! --------------------------------------------------------------------------------------------------
!
! Compute stress for pipe element
!
! In  nbNode           : number of nodes in element
! In  nbFourier        : number of Fourier modes
! In  nbDof            : number of DOF in element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jvf, jdfde, jdfd2, jcoopg, jpoids
    integer(kind=8) :: jvDisp
    real(kind=8) :: dispLoca(nbDof)
    integer(kind=8) :: jvSigm
    real(kind=8) :: sigm(PIPE_TENS_SIZE)
    real(kind=8) :: radiusLayer
    integer(kind=8) :: jvMaterc
    real(kind=8) :: meanTemp
    real(kind=8) :: elasMatr(PIPE_TENS_SIZE, PIPE_TENS_SIZE)
    real(kind=8) :: xpg(PIPE_MAX_NPG)
    real(kind=8) :: phi, zeta
    integer(kind=8) :: nbLayer, nbSect
    integer(kind=8) :: nspgLayer, nspgSect, npg, nspg
    integer(kind=8) :: kspgLayer, kspgSect, kpg, kspg
    integer(kind=8) :: iret
    real(kind=8) :: b(PIPE_TENS_SIZE, nbDof)
    real(kind=8) :: sigmTher(2), epsiTher, prodMatrMatr(PIPE_TENS_SIZE, nbDof)
    type(pipeElem_Prop) :: pipeElem
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', npg=npg, &
                     jpoids=jpoids, jcoopg=jcoopg, jvf=jvf, jdfde=jdfde, jdfd2=jdfd2)
    ASSERT(npg .le. PIPE_MAX_NPG)

! - Output field
    call jevech('PCONTRR', 'E', jvSigm)

! - Get parameters about layers and sections
    call pipeGetSubPoints(nbLayer, nbSect, &
                          nspgLayer, nspgSect, nspg)

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

! - Get elastic properties
    call jevech('PMATERC', 'L', jvMaterc)
    call pipeGetElasProp(zi(jvMaterc), temp_=meanTemp, &
                         c_=elasMatr)

! - Get displacements  (in local base)
    call jevech('PDEPLAR', 'L', jvDisp)
    call pipeGetDisp(pipeElem, nbDof, zr(jvDisp), dispLoca)

! - Loop on Gauss points (on segment)
    kspg = 0
    sigmTher = 0.d0
    do kpg = 1, npg
! ----- Compute thermal strain and stress
        call verifg('RIGI', kpg, nspg, '+', zi(jvMaterc), epsiTher)
        sigmTher(1) = (elasMatr(1, 1)+elasMatr(1, 2))*epsiTher
        sigmTher(2) = (elasMatr(2, 1)+elasMatr(2, 2))*epsiTher

! ----- Loop on sub-points in thickness
        do kspgLayer = 1, nspgLayer
! --------- Radius of layer
            call pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)

! --------- Coordinate of sub-point in thickness
            zeta = (kspgLayer-1)*pipeElem%beamElem%sectPipe%thickness/(2.d0*nbLayer)- &
                   pipeElem%beamElem%sectPipe%thickness/2.d0

! --------- Loop on sub-points in section
            do kspgSect = 1, nspgSect
                kspg = kspg+1

! ------------- Coordinate of sub-point in section
                phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)

! ------------- Compute B matrix (u => epsi)
                call pipeBMatr(pipeElem, nbDof, &
                               npg, xpg, &
                               phi, zeta, radiusLayer, &
                               kpg, zr(jvf), zr(jdfde), zr(jdfd2), &
                               b)

! ------------- Compute [K] [B]
                call promat(elasMatr, PIPE_TENS_SIZE, PIPE_TENS_SIZE, PIPE_TENS_SIZE, b, &
                            PIPE_TENS_SIZE, PIPE_TENS_SIZE, nbDof, prodMatrMatr)

! ------------- Compute {sigm} = [K] [B] {u}
                iret = 0
                call prmave(0, prodMatrMatr, PIPE_TENS_SIZE, PIPE_TENS_SIZE, nbDof, &
                            dispLoca, nbDof, sigm, PIPE_TENS_SIZE, iret)

! ------------- Set values
                zr(jvSigm-1+PIPE_NB_SIGM*(kspg-1)+PIPE_SIGM_SIXX) = sigm(1)-sigmTher(1)
                zr(jvSigm-1+PIPE_NB_SIGM*(kspg-1)+PIPE_SIGM_SIYY) = sigm(2)-sigmTher(2)
                zr(jvSigm-1+PIPE_NB_SIGM*(kspg-1)+PIPE_SIGM_SIZZ) = 0.d0
                zr(jvSigm-1+PIPE_NB_SIGM*(kspg-1)+PIPE_SIGM_SIXY) = sigm(3)
                zr(jvSigm-1+PIPE_NB_SIGM*(kspg-1)+PIPE_SIGM_SIXZ) = sigm(4)
                zr(jvSigm-1+PIPE_NB_SIGM*(kspg-1)+PIPE_SIGM_SIYZ) = 0.d0
            end do
        end do
    end do
!
end subroutine
