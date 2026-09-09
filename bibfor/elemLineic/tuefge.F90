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
subroutine tuefge(nbNode, nbFourier)
!
    use beamElem_type
    use pipeElem_type
    use pipeElem_module
    use beamElem_module
!
    implicit none
!
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/pipeElem_type.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nbNode, nbFourier
!
! --------------------------------------------------------------------------------------------------
!
! Compute generalized forces for pipe element - EFGE_ELGA (non-linear case)
!
! In  nbNode           : number of nodes in element
! In  nbFourier        : number of Fourier modes
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: jvf, jdfde, jdfd2, jcoopg, jpoids
    real(kind=8) :: weightLayer(2*PIPE_MAX_LAYERS+1), weightSect(2*PIPE_MAX_SECTORS+1)
    real(kind=8) :: efgeElga(PIPE_MAX_NPG, PIPE_NBDOF_BEAM)
    integer(kind=8) :: jvSigm, jvEfge
    integer(kind=8) :: nbLayer, nbSect
    integer(kind=8) :: nspgLayer, nspgSect, nspg, npg
    integer(kind=8) :: iDofBeam, kpg
    type(pipeElem_Prop) :: pipeElem
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg, &
                     jpoids=jpoids, jcoopg=jcoopg, jvf=jvf, jdfde=jdfde, jdfd2=jdfd2)

! - Get parameters about layers and sections
    call pipeGetSubPoints(nbLayer, nbSect, &
                          nspgLayer, nspgSect, nspg, &
                          weightLayer, weightSect)

! - Get properties of pipe
    call pipeGetProperties(nbNode, nbFourier, pipeElem)

! - Get stresses
    call jevech('PSIEFR', 'L', jvSigm)

! - Output field
    call jevech('PEFGER', 'E', jvEfge)

! - Compute EFGE_ELGA(non-linear case)
    call pipeEfgeElgaNlin(pipeElem, &
                          jvSigm, npg, &
                          nbSect, nbLayer, &
                          nspgSect, nspgLayer, &
                          weightSect, weightLayer, &
                          efgeElga)

! - Set output values
    do iDofBeam = 1, PIPE_NBDOF_BEAM
        do kpg = 1, npg
            zr(jvEfge+6*(kpg-1)+iDofBeam-1) = efgeElga(kpg, iDofBeam)
        end do
    end do
!
end subroutine
