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
subroutine thmGetParaTher(temp, ds_thm)
!
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/get_elasth_para.h"
#include "asterfort/rcvala.h"
#include "asterfort/THM_type.h"
!
    real(kind=8), intent(in) :: temp
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get thermic parameters
!
! --------------------------------------------------------------------------------------------------
!
! In  temp             : current temperature
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: alpha(2)
    integer(kind=8) :: biot_type
    integer(kind=8), parameter :: nbProp1 = 4
    character(len=16), parameter :: propName1(nbProp1) = (/'LAMB_T ', 'LAMB_TL', &
                                                           'LAMB_TN', 'LAMB_TT'/)
    real(kind=8) :: propVale1(nbProp1)
    integer(kind=8) :: propCode1(nbProp1)
    integer(kind=8), parameter :: nbProp2 = 4
    character(len=16), parameter :: propName2(nbProp2) = (/'D_LB_T ', 'D_LB_TL', &
                                                           'D_LB_TN', 'D_LB_TT'/)
    real(kind=8) :: propVale2(nbProp2)
    integer(kind=8) :: propCode2(nbProp2)
    integer(kind=8), parameter :: nbProp3 = 4
    character(len=16), parameter :: propName3(nbProp3) = (/'LAMB_CT ', 'LAMB_C_L', &
                                                           'LAMB_C_N', 'LAMB_C_T'/)
    real(kind=8) :: propVale3(nbProp3)
    integer(kind=8) :: propCode3(nbProp3)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = ds_thm%ds_behaviour%BEHInteg%materPara

! - Read parameters (mechanic dilatation)
    if (ds_thm%ds_elem%l_dof_ther .and. ds_thm%ds_elem%l_dof_meca) then
        call get_elasth_para(materPara%schemePara%fami, &
                             materPara%jvMaterCode, '+', &
                             materPara%schemePara%kpg, &
                             materPara%schemePara%ksp, &
                             materPara%elasID, materPara%elasKeyword, &
                             temp_vale_=temp, &
                             alpha=alpha, &
                             alpha_l=ds_thm%ds_material%ther%alpha_l, &
                             alpha_t=ds_thm%ds_material%ther%alpha_t, &
                             alpha_n=ds_thm%ds_material%ther%alpha_n)
        ds_thm%ds_material%ther%alpha = alpha(1)
    else
        ds_thm%ds_material%ther%alpha = 0.d0
        ds_thm%ds_material%ther%alpha_l = 0.d0
        ds_thm%ds_material%ther%alpha_t = 0.d0
        ds_thm%ds_material%ther%alpha_n = 0.d0
    end if

! - Read parameters for conductivity
    biot_type = ds_thm%ds_material%biot%type
    if (ds_thm%ds_elem%l_dof_ther) then
        propVale1(:) = 0.d0
        propVale2(:) = 0.d0
        propVale3(:) = 0.d0
        call rcvala(materPara%jvMaterCode, &
                    ' ', 'THM_DIFFU', &
                    1, 'TEMP', [temp], &
                    nbProp1, propName1, propVale1, &
                    propCode1, 0, nan='NON')
        ds_thm%ds_material%ther%lambda = propVale1(1)
        ds_thm%ds_material%ther%lambda_tl = propVale1(2)
        ds_thm%ds_material%ther%lambda_tn = propVale1(3)
        ds_thm%ds_material%ther%lambda_tt = propVale1(4)
        call rcvala(materPara%jvMaterCode, ' ', 'THM_DIFFU', &
                    1, 'TEMP', [temp], &
                    nbProp2, propName2, propVale2, &
                    propCode2, 0, nan='NON')
        ds_thm%ds_material%ther%dlambda = propVale2(1)
        ds_thm%ds_material%ther%dlambda_tl = propVale2(2)
        ds_thm%ds_material%ther%dlambda_tn = propVale2(3)
        ds_thm%ds_material%ther%dlambda_tt = propVale2(4)
        call rcvala(materPara%jvMaterCode, ' ', 'THM_DIFFU', &
                    1, 'TEMP', [temp], &
                    nbProp3, propName3, propVale3, &
                    propCode3, 0, nan='NON')
        ds_thm%ds_material%ther%lambda_ct = propVale3(1)
        ds_thm%ds_material%ther%lambda_ct_l = propVale3(2)
        ds_thm%ds_material%ther%lambda_ct_n = propVale3(3)
        ds_thm%ds_material%ther%lambda_ct_t = propVale3(4)
        if (propCode1(1) .eq. 0) then
            ds_thm%ds_material%ther%cond_type = THER_COND_ISOT
            ASSERT(propCode1(2) .eq. 1)
            ASSERT(propCode1(3) .eq. 1)
            ASSERT(propCode1(4) .eq. 1)
        else
            if (propCode1(4) .eq. 0) then
                ds_thm%ds_material%ther%cond_type = THER_COND_ORTH
            else
                ds_thm%ds_material%ther%cond_type = THER_COND_ISTR
            end if
        end if
    else
        ds_thm%ds_material%ther%cond_type = THER_COND_ISOT
        ds_thm%ds_material%ther%lambda = 0
        ds_thm%ds_material%ther%lambda_tl = 0
        ds_thm%ds_material%ther%lambda_tn = 0
        ds_thm%ds_material%ther%lambda_tt = 0
        ds_thm%ds_material%ther%dlambda = 0
        ds_thm%ds_material%ther%dlambda_tl = 0
        ds_thm%ds_material%ther%dlambda_tn = 0
        ds_thm%ds_material%ther%dlambda_tt = 0
        ds_thm%ds_material%ther%lambda_ct = 0
        ds_thm%ds_material%ther%lambda_ct_l = 0
        ds_thm%ds_material%ther%lambda_ct_n = 0
        ds_thm%ds_material%ther%lambda_ct_t = 0
    end if
!
end subroutine
