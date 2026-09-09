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
subroutine te0049(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epsimc.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sigimc.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
!
! Options: CHAR_MECA_EPSI_*
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8), parameter :: nharm = 0.d0, zero = 0.d0
    real(kind=8) :: sigi(162), epsi(162), bsigmEner(81)
    real(kind=8) :: time
    integer(kind=8) :: i, idfde, jvGeom, jvMaterc, ipoids, iret
    integer(kind=8) :: jvInstr, ivectu, ivf
    integer(kind=8) :: npg, nbsig, ndim, nno
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)

! - Initializations
    nbsig = nbsigm()
    epsi = zero
    sigi = zero
    bsigmEner = zero

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Get current time
    time = r8vide()
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jvInstr)
    if (jvInstr .ne. 0) then
        time = zr(jvInstr)
    end if

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Compute initial strains
    call epsimc(option, zr(jvGeom), nno, npg, ndim, &
                nbsig, zr(ivf), epsi)

! - Compute initial stresses
    call sigimc(materPara, &
                nbsig, npg, time, &
                epsi, sigi)

! - Compute CHAR_MECA_EPSI
    call bsigmc(nno, ndim, nbsig, npg, ipoids, &
                ivf, idfde, zr(jvGeom), nharm, sigi, &
                bsigmEner)

! - Set output
    call jevech('PVECTUR', 'E', ivectu)
    do i = 1, ndim*nno
        zr(ivectu+i-1) = bsigmEner(i)
    end do
!
end subroutine
