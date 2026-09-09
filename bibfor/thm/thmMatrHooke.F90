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
subroutine thmMatrHooke(ds_thm)
!
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/matrHooke3d.h"
#include "asterfort/separ_RI_elas_3D.h"
!
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute Hooke elastic matrix
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: h(6), hi(6), g, e, nu
    real(kind=8) :: e1i, e2i, e3i, gi
    real(kind=8) :: nu12i, nu13i, nu23i, nui
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = ds_thm%ds_behaviour%BEHInteg%materPara

! - Prepare Hooke matrix coefficient
    if (materPara%elasID .eq. ELAS_ISOT) then
        e = ds_thm%ds_material%elas%e
        nu = ds_thm%ds_material%elas%nu
        g = e/(2.d0*(1.d0+nu))
    else
        g = ds_thm%ds_material%elas%g
    end if

    call separ_RI_elas_3D(materPara%elasID, &
                          nu=ds_thm%ds_material%elas%nu, &
                          g=g, nui=nui, gi=gi, &
                          e1=ds_thm%ds_material%elas%e_l, &
                          e2=ds_thm%ds_material%elas%e_t, &
                          e3=ds_thm%ds_material%elas%e_n, &
                          nu12=ds_thm%ds_material%elas%nu_lt, &
                          nu13=ds_thm%ds_material%elas%nu_ln, &
                          nu23=ds_thm%ds_material%elas%nu_tn, &
                          e1i=e1i, e2i=e2i, e3i=e3i, &
                          nu12i=nu12i, nu13i=nu13i, nu23i=nu23i, &
                          hr=h, hi=hi)
! - Compute matrix
    call matrHooke3d(materPara, &
                     h=h, g=g, &
                     g1=ds_thm%ds_material%elas%g_lt, &
                     g2=ds_thm%ds_material%elas%g_ln, &
                     g3=ds_thm%ds_material%elas%g_tn, &
                     matr_elas=ds_thm%ds_material%elas%d)
!
end subroutine
