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
subroutine dmat3d(materPara, poum, time, &
                  dr_, di_)
!
    use MaterialPara_type
    implicit none
!
#include "asterfort/get_elas_para.h"
#include "asterfort/matrHooke3d.h"
#include "asterfort/separ_RI_elas_3D.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=*), intent(in) :: poum
    real(kind=8), intent(in) :: time
    real(kind=8), optional, intent(out) :: dr_(6, 6), di_(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
! Hooke matrix for iso-parametric elements
!
! 3D and Fourier
!
! --------------------------------------------------------------------------------------------------
!
! In  materPara        : parameters of material
! In  poum             : '-' or '+' for parameters evaluation (previous or current temperature)
! In  time             : current time
! Out dr               : Hooke matrix (real part)
! Out di               : Hooke matrix (imaginary part, for viscoelasticity)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: nur, nui, nu12r, nu13r, nu23r, nu12i, nu13i, nu23i
    real(kind=8) :: e1r, e2r, e3r, e1i, e2i, e3i, er, ei
    real(kind=8) :: g1r, g2r, g3r, g1i, g2i, g3i, gr, gi
    real(kind=8) :: di(6, 6), dr(6, 6), hr(6), hi(6)
!
! --------------------------------------------------------------------------------------------------
!

! - Get elastic parameters
    call get_elas_para(materPara%schemePara%fami, materPara%jvMaterCode, poum, &
                       materPara%schemePara%kpg, materPara%schemePara%ksp, &
                       materPara%elasID, materPara%elasKeyword, &
                       time=time, &
                       e_=er, nu_=nur, g_=gr, &
                       e1_=e1r, e2_=e2r, e3_=e3r, &
                       nu12_=nu12r, nu13_=nu13r, nu23_=nu23r, &
                       g1_=g1r, g2_=g2r, g3_=g3r, &
                       ei_=ei, nui_=nui, gi_=gi, &
                       e1i_=e1i, e2i_=e2i, e3i_=e3i, &
                       nu12i_=nu12i, nu13i_=nu13i, nu23i_=nu23i, &
                       g1i_=g1i, g2i_=g2i, g3i_=g3i)

! - Prepare Hooke matrix coefficient
    call separ_RI_elas_3D(materPara%elasID, &
                          nur, gr, nui, gi, &
                          e1r, e2r, e3r, &
                          nu12r, nu13r, nu23r, &
                          e1i, e2i, e3i, &
                          nu12i, nu13i, nu23i, &
                          hr, hi)

! - Compute Hooke matrix
    if (present(di_)) then
        call matrHooke3d(materPara, &
                         hi, gi, g1i, g2i, g3i, &
                         di)
        di_ = di
    end if
    if (present(dr_)) then
        call matrHooke3d(materPara, &
                         hr, gr, g1r, g2r, g3r, &
                         dr)
        dr_ = dr
    end if
!
end subroutine
