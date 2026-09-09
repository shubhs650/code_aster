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
subroutine thmTherElas(ds_thm, mdal, dalal)
!
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/matrot.h"
#include "asterfort/utbtab.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(out) :: mdal(6)
    real(kind=8), intent(out) :: dalal
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute thermic quantities
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! Out mdal             : product [Elas] {alpha}
! Out dalal            : product <alpha> [Elas] {alpha}
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: al(6), tal(3, 3), talg(3, 3), work(6, 6), pass(3, 3), anglNaut(3)
    integer(kind=8) :: i, j
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    mdal = 0.d0
    dalal = 0.d0
    al = 0.d0
    tal = 0.d0
    talg = 0.d0
    work = 0.d0
    pass = 0.d0

    materPara = ds_thm%ds_behaviour%BEHInteg%materPara
    anglNaut = materPara%lcsPara%lcsAngle

! - Get dilatation coefficient
    if (materPara%elasID .eq. ELAS_ISOT) then
        al(1) = ds_thm%ds_material%ther%alpha
        al(2) = ds_thm%ds_material%ther%alpha
        al(3) = ds_thm%ds_material%ther%alpha

    else if (materPara%elasID .eq. ELAS_ISTR) then
        call matrot(anglNaut, pass)
        tal(1, 1) = ds_thm%ds_material%ther%alpha_l
        tal(2, 2) = ds_thm%ds_material%ther%alpha_l
        tal(3, 3) = ds_thm%ds_material%ther%alpha_n
        al(1) = talg(1, 1)
        al(2) = talg(2, 2)
        al(3) = talg(3, 3)
        al(4) = talg(1, 2)
        al(5) = talg(1, 3)
        al(6) = talg(2, 3)

    else if (materPara%elasID .eq. ELAS_ORTH) then
        call matrot(anglNaut, pass)
        tal(1, 1) = ds_thm%ds_material%ther%alpha_l
        tal(2, 2) = ds_thm%ds_material%ther%alpha_t
        tal(3, 3) = ds_thm%ds_material%ther%alpha_n
        call utbtab('ZERO', 3, 3, tal, pass, work, talg)
        al(1) = talg(1, 1)
        al(2) = talg(2, 2)
        al(3) = talg(3, 3)
        al(4) = talg(1, 2)
        al(5) = talg(1, 3)
        al(6) = talg(2, 3)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Compute
    do i = 1, 6
        do j = 1, 6
            mdal(i) = mdal(i)+ds_thm%ds_material%elas%d(i, j)*al(j)
        end do
    end do
    do i = 1, 6
        dalal = dalal+mdal(i)*al(i)
    end do
!
end subroutine
