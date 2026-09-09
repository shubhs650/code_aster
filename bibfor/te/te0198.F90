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
subroutine te0198(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sigtmc.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: AXIS_FOURIER
!
! Option: CHAR_MECA_TEMP_R
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: dimmod = 3
    real(kind=8) :: nharm
    real(kind=8) :: forcTher(3*MT_NNOMAX), sigmTher(162)
    real(kind=8) :: time
    integer(kind=8) :: jvTime, jvVect, jvHarmon
    integer(kind=8) :: ndim, nno, npg, nbsig, i, iret
    integer(kind=8) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    integer(kind=8) :: jvGeom, jvMaterc
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)
    nbsig = nbsigm()
    ASSERT(nno .le. MT_NNOMAX)
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Get time
    time = r8vide()
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jvTime)
    if (jvTime .ne. 0) then
        time = zr(jvTime)
    end if

! - Get Fourier mode
    call jevech('PHARMON', 'L', jvHarmon)
    nharm = dble(zi(jvHarmon))

! - Calcul des contraintes thermiques
    call sigtmc(materPara, time, &
                nbsig, npg, ndim, &
                VARC_STRAIN_TEMP, sigmTher)

! - Compute CHAR_MECA_TEMP_R: [B]Tx{SIGTH}
    call bsigmc(nno, dimmod, nbsig, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                zr(jvGeom), nharm, sigmTher, &
                forcTher)

! - Set output vector
    call jevech('PVECTUR', 'E', jvVect)
    do i = 1, dimmod*nno
        zr(jvVect+i-1) = forcTher(i)
    end do
!
end subroutine
