! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine te0204(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/evalPressure.h"
#include "asterfort/lteatt.h"
#include "asterfort/vff2dn.h"
#include "asterfort/teattr.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D_FLUI_STRU, AXIS_FLUI_STRU
!
! Options: CHAR_MECA_PRES_*
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_func, l_time
    integer :: jv_geom, jv_pres, jv_time, jv_vect
    real(kind=8) :: nx, ny, tx, ty
    real(kind=8) :: poids, time
    real(kind=8) :: pres, cisa
    integer :: ipoids, ivf, idfde
    integer :: nno, npg, ndim, ndofbynode
    integer :: ldec, iret
    integer :: i, ipg
    aster_logical :: l_axis
    real(kind=8) :: r
    character(len=16) :: FEForm
!
! --------------------------------------------------------------------------------------------------
!
    l_func = (option .eq. 'CHAR_MECA_PRES_F')
!
! - Input fields: for pressure, no node affected -> 0
!
    call jevech('PGEOMER', 'L', jv_geom)
    if (l_func) then
        call jevech('PPRESSF', 'L', jv_pres)
    else
        call jevecd('PPRESSR', jv_pres, 0.d0)
    end if
!
! - Get time if present
!
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jv_time)
    l_time = ASTER_FALSE
    time = 0.d0
    if (jv_time .ne. 0) then
        l_time = ASTER_TRUE
        time = zr(jv_time)
    end if
!
! - Get element parameters
!
    l_axis = (lteatt('AXIS', 'OUI'))
    call teattr('S', 'FORMULATION', FEForm, iret)
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, ndim=ndim, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - Pressure are on skin elements but DOF are surfacic + FSI
!
    ASSERT(ndim .eq. 1)
    if (FEForm .eq. 'U_P_PHI') then
        ndofbynode = 3
    elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
        ndofbynode = 3
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if
!
! - Output field
!
    call jevech('PVECTUR', 'E', jv_vect)
    do i = 1, ndofbynode*nno
        zr(jv_vect+i-1) = 0.d0
    end do
!
! - Loop on Gauss points
!
    do ipg = 1, npg
        ldec = (ipg-1)*nno
! ----- Compute normal
        nx = 0.d0
        ny = 0.d0
        call vff2dn(ndim, nno, ipg, ipoids, idfde, &
                    zr(jv_geom), nx, ny, poids)
! ----- Geometry
        if (l_axis) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(jv_geom+2*(i-1))*zr(ivf+ldec+i-1)
            end do
            poids = poids*r
        end if
! ----- Get value of pressure
        call evalPressure(l_func, l_time, time, &
                          nno, ndim, ipg, &
                          ivf, jv_geom, jv_pres, &
                          pres, cisa)
! ----- Compute vector
        tx = -nx*pres-ny*cisa
        ty = -ny*pres+nx*cisa
! ----- Normal convention is fluid -> structure
        tx = -tx
        ty = -ty
        do i = 1, nno
            zr(jv_vect+ndofbynode*(i-1)-1+1) = &
                zr(jv_vect+ndofbynode*(i-1)-1+1)+tx*zr(ivf+ldec+i-1)*poids
            zr(jv_vect+ndofbynode*(i-1)-1+2) = &
                zr(jv_vect+ndofbynode*(i-1)-1+2)+ty*zr(ivf+ldec+i-1)*poids
        end do
    end do
end subroutine
