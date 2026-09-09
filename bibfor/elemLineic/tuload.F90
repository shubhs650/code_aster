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
subroutine tuload(option, nbNode, nbDof, nbFourier)
!
    use pipeElem_type
    use pipeElem_module
    use beamElem_module
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    integer(kind=8), intent(in) :: nbNode, nbDof, nbFourier
!
! --------------------------------------------------------------------------------------------------
!
! Compute loads for pipe element
!
! In  option           : option to compute
! In  nbNode           : number of nodes in element
! In  nbFourier        : number of Fourier modes
! In  nbDof            : number of DOF in element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "MASS"
    integer(kind=8) :: jvf, jcoopg, jpoids
    integer(kind=8) :: jvGeom, jvTime, jvLoad
    integer(kind=8) :: jvVect
    real(kind=8) :: weightLayer(2*PIPE_MAX_LAYERS+1), weightSect(2*PIPE_MAX_SECTORS+1)
    character(len=1) :: typeScal
    integer(kind=8) :: itab(2)
    integer(kind=8) :: jvCurv
    aster_logical :: lAbsCurv
    real(kind=8) :: absCurv(PIPE_MAX_NODE)
    integer(kind=8) :: jvMaterc
    real(kind=8) :: rho
    real(kind=8) :: xpg(PIPE_MAX_NPG)
    integer(kind=8) :: nbLayer, nbSect
    integer(kind=8) :: nspgLayer, nspgSect, npg, nspg
    integer(kind=8) :: kpg
    integer(kind=8) :: iNode, iret
    aster_logical :: lGravity
    type(pipeElem_Prop) :: pipeElem
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg, &
                     jpoids=jpoids, jcoopg=jcoopg, jvf=jvf)
    ASSERT(npg .le. PIPE_MAX_NPG)

! - Output field
    call jevech('PVECTUR', 'E', jvVect)

! - Get initial coordinates of nodes
    call jevech('PGEOMER', 'L', jvGeom)

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

! - Get curvilinear parameter
    lAbsCurv = ASTER_FALSE
    if (option .eq. 'CHAR_MECA_PRES_F') then
        call tecach('ONO', 'PABSCUR', 'L', iret, iad=jvCurv)
        if (jvCurv .ne. 0) then
            ASSERT(iret .eq. 0)
            call tecach('OOO', 'PABSCUR', 'L', iret, nval=2, itab=itab)
            ASSERT(itab(1) .eq. jvCurv)
            ASSERT(itab(2) .eq. nbNode)
            do iNode = 1, nbNode
                absCurv(iNode) = zr(jvCurv-1+iNode)
            end do
            lAbsCurv = ASTER_TRUE
        end if
    end if
!
    jvTime = 0
    if (option(1:14) .eq. 'CHAR_MECA_PRES') then
        typeScal = option(16:16)
        if (typeScal .eq. 'R') then
            call jevecd('PPRESSR', jvLoad, 0.d0)
        else if (typeScal .eq. 'F') then
            call jevech('PPRESSF', 'L', jvLoad)
            call jevech('PINSTR', 'L', jvTime)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Compute load for internal pressure in pipe element
        call pipeLoadPres(typeScal, pipeElem, &
                          nbNode, nbFourier, nbSect, &
                          npg, nspgSect, &
                          jvLoad, jvTime, jvGeom, jpoids, &
                          lAbsCurv, absCurv, zr(jvf), &
                          xpg, weightSect, &
                          nbDof, jvVect)

    else if (option .eq. 'CHAR_MECA_PESA_R') then
        typeScal = "R"
        lGravity = ASTER_TRUE
        call jevech('PPESANR', 'L', jvLoad)
        call jevech('PMATERC', 'L', jvMaterc)
        call pipeGetDensity(jvMaterc, rho)
        call pipeLoadLine(typeScal, lGravity, pipeElem, &
                          nbNode, nbFourier, &
                          nbSect, nbLayer, &
                          npg, nspgSect, nspgLayer, &
                          jvLoad, jvTime, jvGeom, jpoids, &
                          rho, &
                          zr(jvf), weightSect, weightLayer, &
                          nbDof, jvVect)

    else if (option .eq. 'CHAR_MECA_FR1D1D') then
        typeScal = "R"
        lGravity = ASTER_FALSE
        call jevech('PFR1D1D', 'L', jvLoad)
        rho = -1.d0
        call pipeLoadLine(typeScal, lGravity, pipeElem, &
                          nbNode, nbFourier, &
                          nbSect, nbLayer, &
                          npg, nspgSect, nspgLayer, &
                          jvLoad, jvTime, jvGeom, jpoids, &
                          rho, &
                          zr(jvf), weightSect, weightLayer, &
                          nbDof, jvVect)

    else if (option .eq. 'CHAR_MECA_FC1D1D') then
        typeScal = "C"
        lGravity = ASTER_FALSE
        call jevech('PCR1D1D', 'L', jvLoad)
        rho = -1.d0
        call pipeLoadLine(typeScal, lGravity, pipeElem, &
                          nbNode, nbFourier, &
                          nbSect, nbLayer, &
                          npg, nspgSect, nspgLayer, &
                          jvLoad, jvTime, jvGeom, jpoids, &
                          rho, &
                          zr(jvf), weightSect, weightLayer, &
                          nbDof, jvVect)

    else if ((option .eq. 'CHAR_MECA_FF1D1D')) then
        typeScal = "F"
        lGravity = ASTER_FALSE
        call jevech('PFF1D1D', 'L', jvLoad)
        call jevech('PINSTR', 'L', jvTime)
        rho = -1.d0
        call pipeLoadLine(typeScal, lGravity, pipeElem, &
                          nbNode, nbFourier, &
                          nbSect, nbLayer, &
                          npg, nspgSect, nspgLayer, &
                          jvLoad, jvTime, jvGeom, jpoids, &
                          rho, &
                          zr(jvf), weightSect, weightLayer, &
                          nbDof, jvVect)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
