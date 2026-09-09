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

subroutine get_elasth_para(fami, jvMaterCode, poum, kpg, ksp, &
                           elasID, elasKeyword, materi_, temp_vale_, &
                           alpha, alpha_l, alpha_t, alpha_n, &
                           z_h_r_, deps_ch_tref_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/ElasticityMaterial_type.h"
!
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: elasID
    character(len=16), intent(in) :: elasKeyword
    character(len=8), optional, intent(in) :: materi_
    real(kind=8), optional, intent(in) :: temp_vale_
    real(kind=8), optional, intent(out) :: alpha(2)
    real(kind=8), optional, intent(out) :: alpha_l, alpha_t, alpha_n
    real(kind=8), optional, intent(out) :: z_h_r_
    real(kind=8), optional, intent(out) :: deps_ch_tref_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility
!
! Get elastic parameters for thermic dilatation
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  jvMaterCode      : coded material address
! In  poum         : '-' or '+' for parameters evaluation (previous or current temperature)
! In  kpg          : current point gauss
! In  ksp         : current "sous-point" gauss
! In  elasID    : Type of elasticity
!                       1 - Isotropic
!                       2 - Orthotropic
!                       3 - Transverse isotropic
! In  elasKeyword : keyword factor linked to type of elasticity parameters
! In  materi       : name of material if multi-material Gauss point (PMF)
! In  temp_vale    : specifi temperature (example: mean temperature for structural elements)
! Out alpha        : thermic dilatation ratio (isotropic)
!                     if   META -> alpha(1) for hot phasis and alpha(2) for cold phasis
!                     else alpha(1) only
! Out alpha_l      : thermic dilatation ratio - Direction L (Orthotropic/Transverse isotropic)
! Out alpha_t      : thermic dilatation ratio - Direction T (Orthotropic)
! Out alpha_n      : thermic dilatation ratio - Direction N (Orthotropic/Transverse isotropic)
! Out z_h_r        : characterizes the reference metallurgical phase :
!                    z_h_r = 1 --> reference phasis = hot phase
!                    z_h_r = 0 --> reference phasis = cold phase
! Out deps_ch_tref : Compactness difference between the hot phase and the cold phase
!                    at the reference temperature
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropMaxi = 4
    integer(kind=8) :: propCode(nbPropMaxi)
    character(len=16) :: propName(nbPropMaxi)
    real(kind=8) :: propVale(nbPropMaxi)
!
    character(len=8) :: paraName, materi
    character(len=24) :: valk(3)
    real(kind=8) :: paraVale
    integer(kind=8) :: nbProp, nbPara, i
    real(kind=8) :: alpha_c, alpha_f, alpha_a
    integer(kind=8) :: iadzi, iazk24
    real(kind=8) :: z_h_r
    real(kind=8) :: deps_ch_tref
!
! --------------------------------------------------------------------------------------------------
!
    nbPara = 0
    paraName = ' '
    paraVale = 0.d0
    materi = ' '
    if (present(materi_)) then
        materi = materi_
    end if
    if (present(temp_vale_)) then
        nbPara = 1
        paraVale = temp_vale_
        paraName = 'TEMP'
    end if
!
! - Get parameters
!
    if (elasID .eq. ELAS_ISOT .or. elasID .eq. ELAS_SHELL .or. &
        elasID .eq. ELAS_MEMBRANE) then
        if (elasKeyword .eq. 'ELAS_HYPER') then
            call utmess('F', 'COMPOR5_6')
        elseif (elasKeyword .eq. 'ELAS_META') then
            nbProp = 4
            propName(1) = 'C_ALPHA'
            propName(2) = 'F_ALPHA'
            propName(3) = 'PHASE_REFE'
            propName(4) = 'EPSF_EPSC_TREF'
            call rcvalb(fami, kpg, ksp, poum, &
                        jvMaterCode, materi, elasKeyword, &
                        nbPara, paraName, [paraVale], &
                        nbProp, propName, propVale, &
                        propCode, 1)
            alpha_c = propVale(1)
            alpha_f = propVale(2)
            z_h_r = propVale(3)
            deps_ch_tref = propVale(4)
            if (present(alpha)) then
                alpha(1) = alpha_c
                alpha(2) = alpha_f
            end if
            if (present(z_h_r_)) then
                z_h_r_ = z_h_r
            end if
            if (present(deps_ch_tref_)) then
                deps_ch_tref_ = deps_ch_tref
            end if
        else
            nbProp = 1
            propName(1) = 'ALPHA'
            call rcvalb(fami, kpg, ksp, poum, &
                        jvMaterCode, materi, elasKeyword, &
                        nbPara, paraName, [paraVale], &
                        nbProp, propName, propVale, &
                        propCode, 1)
            alpha_a = propVale(1)
            alpha(1) = alpha_a
            alpha(2) = 0.d0
        end if
    elseif (elasID .eq. ELAS_ORTH) then
        nbProp = 3
        propName(1) = 'ALPHA_L'
        propName(2) = 'ALPHA_T'
        propName(3) = 'ALPHA_N'
        call rcvalb(fami, kpg, ksp, poum, &
                    jvMaterCode, materi, elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        alpha_l = propVale(1)
        alpha_t = propVale(2)
        alpha_n = propVale(3)
    elseif (elasID .eq. ELAS_ISTR) then
        nbProp = 2
        propName(1) = 'ALPHA_L'
        propName(2) = 'ALPHA_N'
        call rcvalb(fami, kpg, ksp, poum, &
                    jvMaterCode, materi, elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        alpha_l = propVale(1)
        alpha_n = propVale(2)
    else
        WRITE (6, *) "ELAS: ", elasID, elasKeyword
        ASSERT(.false.)
    end if
!
! - Test
!
    do i = 1, nbProp
        if (propCode(i) .ne. 0) then
            call tecael(iadzi, iazk24)
            valk(1) = zk24(iazk24-1+3)
            valk(2) = 'TEMP'
            valk(3) = propName(i)
            call utmess('F', 'COMPOR5_32', nk=3, valk=valk)
        end if
    end do
!
end subroutine
