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
subroutine tuforc(option, nbNode, nbDof, nbFourier)
!
    use resi_refe_module, only: RESI_REFE
    use beamElem_type
    use pipeElem_type
    use pipeElem_module
    use beamElem_module
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    integer(kind=8), intent(in) :: nbNode, nbDof, nbFourier
!
! --------------------------------------------------------------------------------------------------
!
! Compute force vector for pipe element
!
! In  option           : option to compute
! In  nbNode           : number of nodes in element
! In  nbFourier        : number of Fourier modes
! In  nbDof            : number of DOF in element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: jvf, jdfde, jdfd2, jcoopg, jpoids, jgano
    real(kind=8) :: weightLayer(2*PIPE_MAX_LAYERS+1), weightSect(2*PIPE_MAX_SECTORS+1)
    integer(kind=8) :: jvSigm, jvVect
    real(kind=8) :: xpg(PIPE_MAX_NPG)
    integer(kind=8) :: nbLayer, nbSect
    integer(kind=8) :: nspgLayer, nspgSect, nspg, npg
    integer(kind=8) :: kpg
    integer(kind=8) :: nspgCheck, nval, itab(7)
    integer(kind=8) :: iDof, iret
    real(kind=8) :: sigmRefe
    real(kind=8) :: forcNodaMean(nbDof), forcNoda(nbDof)
    type(pipeElem_Prop) :: pipeElem
    blas_int :: b_incx, b_incy, b_n
    type(RESI_REFE):: refe
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

! - Get coordinates of Gauss points (on segment)
    do kpg = 1, npg
        xpg(kpg) = zr(jcoopg-1+kpg)
    end do

! - Compute options
    if (option .eq. 'FORC_NODA') then
! ----- Get stress
        call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, itab=itab)
        jvSigm = itab(1)
        nspgCheck = itab(7)
        if (nspgCheck .ne. nspg) then
            call utmess('F', 'ELEMENTS_4')
        end if

! ----- Output field is vector
        call jevech('PVECTUR', 'E', jvVect)

! ----- Compute vector of nodal forces
        call pipeForcNoda(pipeElem, &
                          npg, xpg, zr(jpoids), &
                          nbSect, nbLayer, &
                          nspgSect, nspgLayer, &
                          weightSect, weightLayer, &
                          zr(jvf), zr(jdfde), zr(jdfd2), &
                          nbDof, forcNoda, &
                          jvSigm_=jvSigm)

! ----- Save vector
        do iDof = 1, nbDof
            zr(jvVect-1+iDof) = forcNoda(iDof)
        end do

    else if (option .eq. 'REFE_FORC_NODA') then
! ----- Get reference stress
        call refe%Init('MECA_TUYAU')
        sigmRefe = refe%GetRef('SIGM')
        call refe%Check()

! ----- Output field is vector
        call jevech('PVECTUR', 'E', jvVect)

! ----- Compute vector of nodal forces
        call pipeForcNoda(pipeElem, &
                          npg, xpg, zr(jpoids), &
                          nbSect, nbLayer, &
                          nspgSect, nspgLayer, &
                          weightSect, weightLayer, &
                          zr(jvf), zr(jdfde), zr(jdfd2), &
                          nbDof, forcNoda, &
                          sigmRefe_=sigmRefe)

! ----- Compute mean value
        forcNodaMean = 0.d0
        nval = PIPE_TENS_SIZE*npg*(2*nbLayer+1)*(2*nbSect+1)
        b_n = to_blas_int(nbDof)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0/nval, forcNoda, b_incx, forcNodaMean, b_incy)

! ----- Save vector
        do iDof = 1, nbDof
            zr(jvVect-1+iDof) = forcNodaMean(iDof)
        end do

    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
