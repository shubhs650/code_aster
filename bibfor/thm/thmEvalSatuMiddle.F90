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
subroutine thmEvalSatuMiddle(ds_thm, &
                             p1, tempCurr, &
                             satur, dsatur, retcom)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvala.h"
#include "asterfort/THM_type.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: p1, tempCurr
    real(kind=8), intent(out) :: satur, dsatur
    integer(kind=8), intent(out) :: retcom
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Evaluation of "middle" saturation (only LIQU_VAPE)
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  p1               : capillary pressure - At end of current step
! Out satur            : saturation
! Out dsatur           : derivative of saturation (/pc)
! Out retcom           : return code for error
!                         2 - If saturation doesn't belon to ]0,1[
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 2
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'SATU_PRES  ', 'D_SATU_PRES'/)
    integer(kind=8), parameter :: nbPara = 2
    character(len=16), parameter :: paraName(nbPara) = (/'PCAP', 'TEMP'/)
!
! --------------------------------------------------------------------------------------------------
!
    retcom = 0
    propVale = 0.d0
    if (ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_UTIL' .or. &
        ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_ENDO' .or. &
        ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_TABBAL') then
        call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                    ' ', 'THM_DIFFU', &
                    nbPara, paraName, [p1, tempCurr], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        satur = propVale(1)
        dsatur = propVale(2)
        ASSERT(ds_thm%ds_behaviour%satur_type .eq. SATURATED_SPEC)
    else
        ASSERT(ASTER_FALSE)
    end if
    if (satur .gt. 1.d0 .or. satur .lt. 0.d0) then
        retcom = 2
    end if
!
end subroutine
