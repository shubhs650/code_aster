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
subroutine thmGetParaInit(ds_thm, l_check_)
!
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/rcvala.h"
#include "asterfort/THM_type.h"
#include "asterfort/utmess.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    aster_logical, optional, intent(in) :: l_check_
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get initial parameters (THM_INIT)
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  l_check          : check THM_INIT
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_temp_init, l_pre2_init
    integer(kind=8) :: nume_thmc, nume_init
    integer(kind=8), parameter :: nbProp = 5
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'TEMP     ', 'PRE1     ', &
                                                         'PRE2     ', 'PORO     ', &
                                                         'PRES_VAPE'/)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = ds_thm%ds_behaviour%BEHInteg%materPara
    propVale = r8nnem()

! - Read parameters
    call rcvala(materPara%jvMaterCode, &
                ' ', 'THM_INIT', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 0, nan='OUI')
    ds_thm%ds_parainit%temp_init = propVale(1)
    ds_thm%ds_parainit%pre1_init = propVale(2)
    ds_thm%ds_parainit%pre2_init = propVale(3)
    ds_thm%ds_parainit%poro_init = propVale(4)
    ds_thm%ds_parainit%prev_init = propVale(5)
    l_temp_init = propCode(1) .eq. 0
    l_pre2_init = propCode(3) .eq. 0

! - Check: compatibility coupling law with initial parameters
    if (present(l_check_)) then
        nume_thmc = ds_thm%ds_behaviour%nume_thmc
        call rcvala(materPara%jvMaterCode, &
                    ' ', 'THM_INIT', &
                    0, ' ', [0.d0], &
                    1, 'COMP_THM', propVale, &
                    propCode, 1)
        nume_init = nint(propVale(1))
        if (nume_init .ne. nume_thmc) then
            call utmess('F', 'THM1_34')
        end if
    end if
!
end subroutine
