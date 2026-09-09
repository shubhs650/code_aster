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
subroutine tuefgeElno(lLine, nbNode, nbDof, nbFourier)
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
    aster_logical, intent(in) :: lLine
    integer(kind=8), intent(in) :: nbNode, nbDof, nbFourier
!
! --------------------------------------------------------------------------------------------------
!
! Compute generalized forces for pipe element - EFGE_ELNO
!
! In  lLine            : flag for linear case (elastic)
! In  nbNode           : number of nodes in element
! In  nbFourier        : number of Fourier modes
! In  nbDof            : number of DOF in element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: jvf, jdfde, jdfd2, jcoopg, jpoids, jgano
    integer(kind=8) :: jvDisp
    real(kind=8) :: weightLayer(2*PIPE_MAX_LAYERS+1), weightSect(2*PIPE_MAX_SECTORS+1)
    integer(kind=8) :: jvMaterCode
    integer(kind=8) :: jvEfgeElno
    real(kind=8) :: efgeElga(PIPE_MAX_NPG, PIPE_NBDOF_BEAM)
    real(kind=8) :: efgeElno(PIPE_NBDOF_BEAM*PIPE_MAX_NODE)
    integer(kind=8) :: jvSigm
    real(kind=8) :: xpg(PIPE_MAX_NPG)
    integer(kind=8) :: nbLayer, nbSect
    integer(kind=8) :: nspgLayer, nspgSect, nspg, npg
    integer(kind=8) :: kpg, iDof
    type(pipeElem_Prop) :: pipeElem
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg, jgano=jgano, &
                     jpoids=jpoids, jcoopg=jcoopg, jvf=jvf, jdfde=jdfde, jdfd2=jdfd2)

! - Get parameters about layers and sections
    call pipeGetSubPoints(nbLayer, nbSect, &
                          nspgLayer, nspgSect, nspg, &
                          weightLayer, weightSect)

! - Get properties of pipe
    call pipeGetProperties(nbNode, nbFourier, pipeElem)

! - Acces to EFGE_ELNO
    call jevech('PEFFORR', 'E', jvEfgeElno)

! - Compute EFGE_ELGA
    if (lLine) then
! ----- Get coordinates of Gauss points (on segment)
        do kpg = 1, npg
            xpg(kpg) = zr(jcoopg-1+kpg)
        end do

! ----- Get material properties
        call jevech('PMATERC', 'L', jvMaterCode)

! ----- Get displacements
        call jevech('PDEPLAR', 'L', jvDisp)

! ----- Compute EFGE_ELGA (linear case)
        call pipeEfgeElgaLine(pipeElem, &
                              jvMaterCode, jvDisp, &
                              nbDof, nspg, npg, xpg, &
                              nbSect, nbLayer, &
                              nspgSect, nspgLayer, &
                              weightSect, weightLayer, &
                              zr(jvf), zr(jdfde), zr(jdfd2), &
                              efgeElga)
    else
! ----- Get stresses
        call jevech('PCONTRR', 'L', jvSigm)

! ----- Compute EFGE_ELGA (non-linear case)
        call pipeEfgeElgaNlin(pipeElem, &
                              jvSigm, npg, &
                              nbSect, nbLayer, &
                              nspgSect, nspgLayer, &
                              weightSect, weightLayer, &
                              efgeElga)
    end if

! - Compute EFGE_ELNO from EFGE_ELGA
    call pipeEfgeElno(pipeElem, &
                      jgano, nbNode, npg, xpg, zr(jvf), &
                      efgeElga, efgeElno)

! - Save EFGE_ELNO
    call jevech('PEFFORR', 'E', jvEfgeElno)
    do iDof = 1, PIPE_NBDOF_BEAM*nbNode
        zr(jvEfgeElno-1+iDof) = efgeElno(iDof)
    end do
!
end subroutine
