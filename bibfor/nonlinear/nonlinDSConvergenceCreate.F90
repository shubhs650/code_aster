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
subroutine nonlinDSConvergenceCreate(ds_conv)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8nnem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/SetResi.h"
!
    type(NL_DS_Conv), intent(out) :: ds_conv
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Create convergence management datastructure
!
! --------------------------------------------------------------------------------------------------
!
! Out ds_conv          : datastructure for convergence management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_resi
!
! --------------------------------------------------------------------------------------------------
!

!
! - Checks
!
    nb_resi = 7
    ds_conv%nb_resi = nb_resi
    ASSERT(nb_resi .le. ds_conv%nb_resi_maxi)
!
! - Set type of residuals
!
    ds_conv%list_resi(1)%type = 'RESI_GLOB_RELA'
    ds_conv%list_resi(2)%type = 'RESI_GLOB_MAXI'
    ds_conv%list_resi(3)%type = 'RESI_REFE_RELA'
    ds_conv%list_resi(4)%type = 'RESI_COMP_RELA'
    ds_conv%list_resi(5)%type = 'RESI_FROT'
    ds_conv%list_resi(6)%type = 'RESI_GEOM'
    ds_conv%list_resi(7)%type = 'RESI_PENE'
!
! - Set name of columns in convergence table (for values)
!
    ds_conv%list_resi(1)%col_name = 'RESI_RELA'
    ds_conv%list_resi(2)%col_name = 'RESI_MAXI'
    ds_conv%list_resi(3)%col_name = 'RESI_REFE'
    ds_conv%list_resi(4)%col_name = 'RESI_COMP'
    ds_conv%list_resi(5)%col_name = 'FROT_NEWT'
    ds_conv%list_resi(6)%col_name = 'GEOM_NEWT'
    ds_conv%list_resi(7)%col_name = 'PENE_MAXI'
!
! - Set name of columns in convergence table (for locus)
!
    ds_conv%list_resi(1)%col_name_locus = 'RELA_NOEU'
    ds_conv%list_resi(2)%col_name_locus = 'MAXI_NOEU'
    ds_conv%list_resi(3)%col_name_locus = 'REFE_NOEU'
    ds_conv%list_resi(4)%col_name_locus = 'COMP_NOEU'
    ds_conv%list_resi(5)%col_name_locus = 'FROT_NOEU'
    ds_conv%list_resi(6)%col_name_locus = 'GEOM_NOEU'
    ds_conv%list_resi(7)%col_name_locus = '         '
!
! - Set event for divergence
!
    ds_conv%list_resi(1)%eventType = 'DIVE_RELA'
    ds_conv%list_resi(2)%eventType = 'DIVE_MAXI'
    ds_conv%list_resi(3)%eventType = 'DIVE_REFE'
    ds_conv%list_resi(4)%eventType = 'DIVE_COMP'
    ds_conv%list_resi(5)%eventType = 'DIVE_FROT'
    ds_conv%list_resi(6)%eventType = 'DIVE_GEOM'
    ds_conv%list_resi(7)%eventType = 'DIVE_PENE'
!
! - Initializations for all residuals
!
    call SetResi(ds_conv, &
                 vale_calc_=r8vide(), locus_calc_=' ', user_para_=r8vide(), &
                 l_conv_=ASTER_FALSE, l_resi_test_=ASTER_FALSE)
!
! - Initializations for reference residual
!
    ds_conv%cresiref = ' '
    ds_conv%cresicmp = ' '
!
! - Other convergence parameters
!
    ds_conv%iter_glob_maxi = 0
    ds_conv%iter_glob_elas = 0
    ds_conv%l_stop = ASTER_TRUE
    ds_conv%l_stop_pene = ASTER_TRUE
!
! - Parameters for automatic swap of convergence criterias
!
    ds_conv%swap_trig = 0.d0

!
! - Parameters for line search
!
    ds_conv%line_sear_coef = r8vide()
    ds_conv%line_sear_iter = 1
!
end subroutine
