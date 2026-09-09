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
subroutine thmGetParaBJH(ds_thm, p1)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvala.h"
#include "asterfort/THM_type.h"
#include "asterfort/utmess.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    real(kind=8), intent(in) :: p1
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Evaluation of BJH Parameters
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  p1               : capillary pressure - At end of current step
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
    real(kind=8) :: ep, surf, sbjh, wbjh
!
! --------------------------------------------------------------------------------------------------
    propValeBJH = 0.d0
    ep = 0.d0
    surf = 0.d0

    if ((ds_thm%ds_behaviour%rela_hydr) .eq. 'HYDR_TABBAL') then
        call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                    ' ', 'THM_DIFFU', &
                    1, 'PCAP', [p1], &
                    nbPropBJH, propNameBJH, propValeBJH, &
                    propCodeBJH, 1)
        surf = propValeBJH(1)
        ds_thm%ds_material%bjh%shuttle = propValeBJH(2)

        ep = propValeBJH(3)
        sbjh = propValeBJH(4)
        wbjh = propValeBJH(5)

        if (surf .lt. 0.d0) then
            call utmess('F', 'THM1_95')
        else
            ds_thm%ds_material%bjh%A0 = surf
        end if

        if (ep .lt. 0.d0) then
            call utmess('F', 'THM1_96')
        else
            ds_thm%ds_material%bjh%epai = ep
        end if

        if (sbjh .lt. 0.d0 .or. sbjh .gt. 1.d0 .or. wbjh .lt. 0.d0 .or. wbjh .gt. 1.d0) then
            call utmess('F', 'THM1_97')
        else
            ds_thm%ds_material%bjh%SBJH = sbjh
            ds_thm%ds_material%bjh%wBJH = wbjh
        end if
    end if

end subroutine
