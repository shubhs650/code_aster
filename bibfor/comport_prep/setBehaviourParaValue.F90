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
subroutine setBehaviourParaValue(prepCrit, parm_theta_thm, parm_alpha_thm, &
                                 iFactorKeyword_, carcriList_, carcriMap_)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/setMFrontPara.h"
!
    type(BehaviourPrep_Crit), pointer :: prepCrit(:)
    real(kind=8), intent(in) :: parm_theta_thm, parm_alpha_thm
    integer(kind=8), optional, intent(in) :: iFactorKeyword_
    real(kind=8), intent(out), optional :: carcriList_(:)
    real(kind=8), pointer, optional :: carcriMap_(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Set values in the map or in list
!
! --------------------------------------------------------------------------------------------------
!
! Ptr prepCrit         : pointer to behaviour criteria
! In  iFactorKeyword   : index of factor keyword (for map)
! In  carcriList       : list for parameters for integration of constitutive law
! In  carcriMap        : map for parameters for integration of constitutive law
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
    iFactorKeyword = 1
    if (present(iFactorKeyword_)) then
        iFactorKeyword = iFactorKeyword_
    end if
!
    if (present(carcriMap_)) then
        if (associated(prepCrit(iFactorKeyword)%iter_inte_maxi)) &
            carcriMap_(ITER_INTE_MAXI) = dble(prepCrit(iFactorKeyword)%iter_inte_maxi)
        carcriMap_(TYPE_MATR_T) = prepCrit(iFactorKeyword)%type_matr_t
        if (associated(prepCrit(iFactorKeyword)%resi_inte)) &
            carcriMap_(RESI_INTE) = prepCrit(iFactorKeyword)%resi_inte
        carcriMap_(PARM_THETA) = prepCrit(iFactorKeyword)%parm_theta
        carcriMap_(ITER_INTE_PAS) = prepCrit(iFactorKeyword)%iter_inte_pas
        carcriMap_(ALGO_INTE_R) = prepCrit(iFactorKeyword)%algo_inte_r
        carcriMap_(VALE_PERT_RELA) = prepCrit(iFactorKeyword)%vale_pert_rela
        carcriMap_(RESI_DEBORST_MAX) = prepCrit(iFactorKeyword)%resi_deborst_max
        carcriMap_(ITER_DEBORST_MAX) = prepCrit(iFactorKeyword)%iter_deborst_max
        carcriMap_(RESI_RADI_RELA) = prepCrit(iFactorKeyword)%resi_radi_rela
        carcriMap_(IVARIEXT1) = prepCrit(iFactorKeyword)%jvariext1
        carcriMap_(IVARIEXT2) = prepCrit(iFactorKeyword)%jvariext2
        carcriMap_(PARM_THETA_THM) = parm_theta_thm
        carcriMap_(PARM_ALPHA_THM) = parm_alpha_thm
        carcriMap_(IPOSTITER) = prepCrit(iFactorKeyword)%ipostiter
        if (prepCrit(iFactorKeyword)%l_matr_unsymm) then
            carcriMap_(CARCRI_MATRSYME) = 1
        else
            carcriMap_(CARCRI_MATRSYME) = 0
        end if
        carcriMap_(EXTE_TYPE) = prepCrit(iFactorKeyword)%prepExte%solvBehavType
        carcriMap_(EXTE_PTR) = prepCrit(iFactorKeyword)%prepExte%adrsUMAT
        carcriMap_(EXTE_STRAIN) = prepCrit(iFactorKeyword)%prepExte%strainMGIS
    end if
    if (present(carcriList_)) then
        if (associated(prepCrit(iFactorKeyword)%iter_inte_maxi)) &
            carcriList_(ITER_INTE_MAXI) = dble(prepCrit(iFactorKeyword)%iter_inte_maxi)
        carcriList_(TYPE_MATR_T) = prepCrit(iFactorKeyword)%type_matr_t
        if (associated(prepCrit(iFactorKeyword)%resi_inte)) &
            carcriList_(RESI_INTE) = prepCrit(iFactorKeyword)%resi_inte
        carcriList_(PARM_THETA) = prepCrit(iFactorKeyword)%parm_theta
        carcriList_(ITER_INTE_PAS) = prepCrit(iFactorKeyword)%iter_inte_pas
        carcriList_(ALGO_INTE_R) = prepCrit(iFactorKeyword)%algo_inte_r
        carcriList_(VALE_PERT_RELA) = prepCrit(iFactorKeyword)%vale_pert_rela
        carcriList_(RESI_DEBORST_MAX) = prepCrit(iFactorKeyword)%resi_deborst_max
        carcriList_(ITER_DEBORST_MAX) = prepCrit(iFactorKeyword)%iter_deborst_max
        carcriList_(RESI_RADI_RELA) = prepCrit(iFactorKeyword)%resi_radi_rela
        carcriList_(IVARIEXT1) = prepCrit(iFactorKeyword)%jvariext1
        carcriList_(IVARIEXT2) = prepCrit(iFactorKeyword)%jvariext2
        carcriList_(PARM_THETA_THM) = parm_theta_thm
        carcriList_(PARM_ALPHA_THM) = parm_alpha_thm
        carcriList_(IPOSTITER) = prepCrit(iFactorKeyword)%ipostiter
        if (prepCrit(iFactorKeyword)%l_matr_unsymm) then
            carcriList_(CARCRI_MATRSYME) = 1
        else
            carcriList_(CARCRI_MATRSYME) = 0
        end if
        carcriList_(EXTE_TYPE) = prepCrit(iFactorKeyword)%prepExte%solvBehavType
        carcriList_(EXTE_PTR) = prepCrit(iFactorKeyword)%prepExte%adrsUMAT
        carcriList_(EXTE_STRAIN) = prepCrit(iFactorKeyword)%prepExte%strainMGIS
    end if

! - Set values for MFRONT
    call setMFrontPara(prepCrit, iFactorKeyword)
!
end subroutine
