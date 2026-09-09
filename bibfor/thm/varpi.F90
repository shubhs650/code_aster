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
subroutine varpi(ds_thm, &
                 p1, p1m, dp1, dp2, &
                 ep, surf, shut, &
                 phi0, dpi, sbjhm, &
                 wbjhm, epm, sbjh, wbjh)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterfort/rcvala.h"
#include "asterfort/THM_type.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: p1, p1m, dp1, dp2
    real(kind=8), intent(in) :: phi0
    real(kind=8), intent(in) :: ep, surf, shut, sbjh, wbjh
    real(kind=8), intent(out) :: dpi
    real(kind=8), intent(out) :: sbjhm, wbjhm, epm
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute the variation of the hydraulic pressure
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  p1m              : capillary pressure - At beginning of step
! In  p1               : capillary pressure - At end of current step!
! In  dp1              : increment of capillary pressure
! In  dp2              : increment of gaz pressure
! In  phi0             : initial porosity (THM_INIT)
! In  ep               : thickness of the adsorbed water layer
! In  surf             : specific surface of the material
! In  shut             : shuttleworth parameter
! In  sbjh             : saturated pores volume fraction from BJH - At end of current step
! In  wbjh             : unsaturatred pores surface fraction from BJH - At end of current step
! In  sbjhm            : saturated pores volume fraction from BJH - At beginning of step
! In  wbjhm            : unsaturatred pores surface fraction from BJH - At beginning of step
! In  epm              : thickness of the adsorbed water layer  - At beginning of step
! Out dpi              : variation of the hydraulic pressure at end of current time
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropBJH = 5
    real(kind=8) :: propValeBJH(nbPropBJH)
    integer(kind=8) :: propCodeBJH(nbPropBJH)
    character(len=16), parameter :: propNameBJH(nbPropBJH) = (/'A0     ', &
                                                               'SHUTTLE', &
                                                               'EPAI   ', &
                                                               'S_BJH  ', &
                                                               'W_BJH  '/)
!
! --------------------------------------------------------------------------------------------------
!
    dpi = 0.d0

! Value of sBJH and wbjh at beginning of step
    call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                ' ', 'THM_DIFFU', &
                1, 'PCAP', [p1m], &
                nbPropBJH, propNameBJH, propValeBJH, &
                propCodeBJH, 1)
    sbjhm = propValeBJH(4)
    wbjhm = propValeBJH(5)
    epm = propValeBJH(3)

    dpi = dp2-(sbjh*p1)+(sbjhm*p1m)+((2./3.)*(0.5*(p1+p1m)*(sbjh-sbjhm))) &
          +((2./3.)*(surf/phi0)*(((wbjh*ep)+(wbjhm*epm))*0.5)*(-dp1)) &
          -((2./3.)*(surf/phi0)*shut*(wbjh-wbjhm))

end subroutine
