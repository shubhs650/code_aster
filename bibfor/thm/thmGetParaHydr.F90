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
subroutine thmGetParaHydr(ds_thm)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "asterfort/rcvala.h"
#include "asterfort/THM_type.h"
!
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get hydraulic parameters (for Mualem-Van Genuchten)
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropVGM = 6
    real(kind=8) :: propValeVGM(nbPropVGM)
    integer(kind=8) :: propCodeVGM(nbPropVGM)
    character(len=16), parameter :: propNameVGM(nbPropVGM) = (/'VG_N    ', &
                                                               'VG_PR   ', &
                                                               'VG_SR   ', &
                                                               'VG_SMAX ', &
                                                               'VG_SATUR', &
                                                               'VG_PENTR'/)
    integer(kind=8), parameter :: nbPropEmmag = 1
    real(kind=8) :: paraValeEmmag(nbPropEmmag)
    integer(kind=8) :: paraCodeEmmag(nbPropEmmag)
    character(len=16), parameter :: paraNameEmmag(nbPropEmmag) = (/'EMMAG'/)
    character(len=16) :: hydr
!
! --------------------------------------------------------------------------------------------------
!
    propValeVGM = 0.d0
    hydr = ds_thm%ds_behaviour%rela_hydr
    if ((hydr .eq. 'HYDR_VGM') .or. (hydr .eq. 'HYDR_VGC')) then
        call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                    ' ', 'THM_DIFFU', &
                    0, ' ', [0.d0], &
                    nbPropVGM, propNameVGM, propValeVGM, &
                    propCodeVGM, 1)
        ds_thm%ds_material%hydr%n = propValeVGM(1)
        ds_thm%ds_material%hydr%pr = propValeVGM(2)
        ds_thm%ds_material%hydr%sr = propValeVGM(3)
        ds_thm%ds_material%hydr%smax = propValeVGM(4)
        ds_thm%ds_material%hydr%satuma = propValeVGM(5)
        ds_thm%ds_material%hydr%pentree = propValeVGM(6)
        if (propCodeVGM(1) .eq. 1) then
            call utmess('F', 'THM1_94')
        end if
    end if

! - For storing coefficient
    paraValeEmmag = 0.d0
    call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                ' ', 'THM_DIFFU', &
                0, ' ', [0.d0], &
                nbPropEmmag, paraNameEmmag, paraValeEmmag, &
                paraCodeEmmag, 0, nan='NON')
    if (paraCodeEmmag(1) .eq. 0) then
        ds_thm%ds_material%hydr%l_emmag = ASTER_TRUE
        if (ds_thm%ds_elem%l_dof_meca) then
            call utmess('F', 'THM1_5')
        end if
    end if
    ds_thm%ds_material%hydr%emmag = paraValeEmmag(1)
!
end subroutine
