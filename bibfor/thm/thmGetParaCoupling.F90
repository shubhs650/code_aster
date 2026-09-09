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
subroutine thmGetParaCoupling(ds_thm, temp)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvala.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    real(kind=8), intent(in) :: temp
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get coupling parameters
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  temp             : current temperature
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropL1 = 4
    real(kind=8) :: propValeL1(nbPropL1)
    integer(kind=8) :: propCodeL1(nbPropL1)
    character(len=16), parameter :: propNameL1(nbPropL1) = (/'RHO        ', 'UN_SUR_K   ', &
                                                             'VISC       ', 'D_VISC_TEMP'/)
    integer(kind=8), parameter :: nbPropL2 = 6
    real(kind=8) :: propValeL2(nbPropL2)
    integer(kind=8) :: propCodeL2(nbPropL2)
    character(len=16), parameter :: propNameL2(nbPropL2) = (/'RHO        ', 'UN_SUR_K   ', &
                                                             'ALPHA      ', 'CP         ', &
                                                             'VISC       ', 'D_VISC_TEMP'/)
    integer(kind=8), parameter :: nbPropG1 = 3
    real(kind=8) :: propValeG1(nbPropG1)
    integer(kind=8) :: propCodeG1(nbPropG1)
    character(len=16), parameter :: propNameG1(nbPropG1) = (/'MASS_MOL   ', &
                                                             'VISC       ', 'D_VISC_TEMP'/)
    integer(kind=8), parameter :: nbPropG2 = 4
    real(kind=8) :: propValeG2(nbPropG2)
    integer(kind=8) :: propCodeG2(nbPropG2)
    character(len=16), parameter :: propNameG2(nbPropG2) = (/'MASS_MOL   ', 'CP         ', &
                                                             'VISC       ', 'D_VISC_TEMP'/)
    integer(kind=8), parameter :: nbPropS = 4
    real(kind=8) :: propValeS(nbPropS)
    integer(kind=8) :: propCodeS(nbPropS)
    character(len=16), parameter :: propNameS(nbPropS) = (/'MASS_MOL   ', 'CP         ', &
                                                           'VISC       ', 'D_VISC_TEMP'/)
    integer(kind=8), parameter :: nbPropAd = 2
    real(kind=8) :: propValeAd(nbPropAd)
    integer(kind=8) :: propCodeAd(nbPropAd)
    character(len=16), parameter :: propNameAd(nbPropAd) = (/'COEF_HENRY ', 'CP         '/)
    integer(kind=8), parameter :: nbPropS1 = 1
    real(kind=8) :: propValeS1(nbPropS1)
    integer(kind=8) :: propCodeS1(nbPropS1)
    character(len=16), parameter :: propNameS1(nbPropS1) = (/'RHO        '/)
    integer(kind=8), parameter :: nbPropS2 = 2
    real(kind=8) :: propValeS2(nbPropS2)
    integer(kind=8) :: propCodeS2(nbPropS2)
    character(len=16), parameter :: propNameS2(nbPropS2) = (/'RHO        ', 'R_GAZ      '/)
    integer(kind=8), parameter :: nbPropS3 = 1
    real(kind=8) :: propValeS3(nbPropS3)
    integer(kind=8) :: propCodeS3(nbPropS3)
    character(len=16), parameter :: propNameS3(nbPropS3) = (/'CP         '/)
    integer(kind=8) :: jvMaterCode
!
! --------------------------------------------------------------------------------------------------
!
    jvMaterCode = ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode
    propValeL1 = 0.d0
    propValeL2 = 0.d0
    propValeG1 = 0.d0
    propValeG2 = 0.d0
    propValeS = 0.d0
    propValeAd = 0.d0
    propValeS1 = 0.d0
    propValeS2 = 0.d0
!
    if (ds_thm%ds_material%l_liquid) then
        if (ds_thm%ds_behaviour%l_temp) then
            call rcvala(jvMaterCode, &
                        ' ', 'THM_LIQU', &
                        1, 'TEMP', [temp], &
                        nbPropL2, propNameL2, propValeL2, &
                        propCodeL2, 1, nan='NON')
            ds_thm%ds_material%liquid%rho = propValeL2(1)
            ds_thm%ds_material%liquid%unsurk = propValeL2(2)
            ds_thm%ds_material%liquid%alpha = propValeL2(3)
            ds_thm%ds_material%liquid%cp = propValeL2(4)
            ds_thm%ds_material%liquid%visc = propValeL2(5)
            ds_thm%ds_material%liquid%dvisc_dtemp = propValeL2(6)
        else
            call rcvala(jvMaterCode, &
                        ' ', 'THM_LIQU', &
                        0, ' ', [0.d0], &
                        nbPropL1, propNameL1, propValeL1, &
                        propCodeL1, 1, nan='NON')
            ds_thm%ds_material%liquid%rho = propValeL1(1)
            ds_thm%ds_material%liquid%unsurk = propValeL1(2)
            ds_thm%ds_material%liquid%visc = propValeL1(3)
            ds_thm%ds_material%liquid%dvisc_dtemp = propValeL1(4)
        end if
    end if
    if (ds_thm%ds_material%l_gaz) then
        if (ds_thm%ds_behaviour%l_temp) then
            call rcvala(jvMaterCode, &
                        ' ', 'THM_GAZ', &
                        1, 'TEMP', [temp], &
                        nbPropG2, propNameG2, propValeG2, &
                        propCodeG2, 1, nan='NON')
            ds_thm%ds_material%gaz%mass_mol = propValeG2(1)
            ds_thm%ds_material%gaz%cp = propValeG2(2)
            ds_thm%ds_material%gaz%visc = propValeG2(3)
            ds_thm%ds_material%gaz%dvisc_dtemp = propValeG2(4)
        else
            call rcvala(jvMaterCode, &
                        ' ', 'THM_GAZ', &
                        0, ' ', [0.d0], &
                        nbPropG1, propNameG1, propValeG1, &
                        propCodeG1, 1, nan='NON')
            ds_thm%ds_material%gaz%mass_mol = propValeG1(1)
            ds_thm%ds_material%gaz%visc = propValeG1(2)
            ds_thm%ds_material%gaz%dvisc_dtemp = propValeG1(3)
        end if
    end if
    if (ds_thm%ds_material%l_steam) then
        call rcvala(jvMaterCode, &
                    ' ', 'THM_VAPE_GAZ', &
                    0, ' ', [0.d0], &
                    nbPropS, propNameS, propValeS, &
                    propCodeS, 1, nan='NON')
        ds_thm%ds_material%steam%mass_mol = propValeS(1)
        ds_thm%ds_material%steam%cp = propValeS(2)
        ds_thm%ds_material%steam%visc = propValeS(3)
        ds_thm%ds_material%steam%dvisc_dtemp = propValeS(4)
    end if
    if (ds_thm%ds_material%l_ad) then
        call rcvala(jvMaterCode, &
                    ' ', 'THM_AIR_DISS', &
                    1, 'TEMP', [temp], &
                    nbPropAd, propNameAd, propValeAd, &
                    propCodeAd, 1, nan='NON')
        ds_thm%ds_material%ad%coef_henry = propValeAd(1)
        ds_thm%ds_material%ad%cp = propValeAd(2)
    end if
    if (ds_thm%ds_material%l_r_gaz) then
        call rcvala(jvMaterCode, &
                    ' ', 'THM_DIFFU', &
                    1, 'TEMP', [temp], &
                    nbPropS2, propNameS2, propValeS2, &
                    propCodeS2, 1, nan='NON')
        ds_thm%ds_material%solid%rho = propValeS2(1)
        ds_thm%ds_material%solid%r_gaz = propValeS2(2)
    else
        call rcvala(jvMaterCode, &
                    ' ', 'THM_DIFFU', &
                    1, 'TEMP', [temp], &
                    nbPropS1, propNameS1, propValeS1, &
                    propCodeS1, 1, nan='NON')
        ds_thm%ds_material%solid%rho = propValeS1(1)
    end if
    if (ds_thm%ds_behaviour%l_temp) then
        call rcvala(jvMaterCode, &
                    ' ', 'THM_DIFFU', &
                    1, 'TEMP', [temp], &
                    nbPropS3, propNameS3, propValeS3, &
                    propCodeS3, 1, nan='NON')
        ds_thm%ds_material%solid%cp = propValeS3(1)
    end if
!
end subroutine
