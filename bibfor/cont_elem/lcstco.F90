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
subroutine lcstco(l_upda_jaco, l_norm_smooth, i_reso_geom, &
                  lagrc_curr, gap_curr, &
                  indi_cont, &
                  gapi, nmcp, &
                  nb_poin_inte, poin_inte_sl, poin_inte_ma)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/assert.h"
#include "asterf_types.h"
#include "Contact_type.h"
!
    aster_logical, intent(out) :: l_upda_jaco
    aster_logical, intent(out) :: l_norm_smooth
    integer(kind=8), intent(out) :: i_reso_geom
    real(kind=8), intent(out) :: lagrc_curr
    real(kind=8), intent(out) :: gap_curr
    integer(kind=8), intent(out) :: indi_cont
    real(kind=8), intent(out) :: gapi
    integer(kind=8), intent(out) :: nmcp
    integer(kind=8), intent(out) :: nb_poin_inte
    real(kind=8), intent(out) :: poin_inte_sl(2, 8)
    real(kind=8), intent(out) :: poin_inte_ma(2, 8)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (LAC) - Elementary computations
!
! Get indicators
!
! --------------------------------------------------------------------------------------------------
!
! Out l_upda_jaco      : .true. to use updated jacobian
! Out l_norm_smooth    : indicator for normals smoothing
! Out i_reso_geom      : algorithm for geometry
! Out lagrc_curr       : current value of contact lagrangian
! Out gap_curr         : current value of gap
! Out indi_cont        : contact indicator
!                        -1 No pairing
!                         0 Paired - No contact
!                        +1 Paired - Contact
! Out gapi             : gap integral on patch
! Out nmcp             : number of contact elements associated to the concerned patch
! Out nb_poin_inte     : number of intersection points
! Out point_inte_sl    : intersection points for slave side
! Out point_inte_ma    : intersection points for master side
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jpcf, ino

! --------------------------------------------------------------------------------------------------
!
    call jevech('PCONFR', 'L', jpcf)
!
    poin_inte_sl = 0.d0
    poin_inte_ma = 0.d0
    indi_cont = nint(zr(jpcf-1+5))
    lagrc_curr = (zr(jpcf-1+6))
    l_upda_jaco = nint(zr(jpcf-1+2)) .eq. 1
    gap_curr = (zr(jpcf-1+15))
    gapi = (zr(jpcf-1+4))
    nmcp = nint(zr(jpcf-1+3))
    l_norm_smooth = nint(zr(jpcf-1+1)) .eq. 1
    nb_poin_inte = nint(zr(jpcf-1+8))
    i_reso_geom = nint(zr(jpcf-1+41))
    ASSERT(i_reso_geom .eq. ALGO_NEWT)
    do ino = 1, nb_poin_inte
        poin_inte_sl(1:2, ino) = zr(jpcf-1+8+2*(ino-1)+1:jpcf-1+8+2*ino)
        poin_inte_ma(1:2, ino) = zr(jpcf-1+24+2*(ino-1)+1:jpcf-1+24+2*ino)
    end do
!
end subroutine
