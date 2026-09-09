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
subroutine fonoda(ds_thm, &
                  ndim, fnoevo, &
                  mecani, press1, press2, tempe, second, &
                  dimdef, dimcon, dt, congem, congep, &
                  r)
!
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/thmEvalGravity.h"
!
    type(THM_DS), intent(in) :: ds_thm
    integer(kind=8), intent(in) :: ndim
    aster_logical, intent(in) :: fnoevo
    integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
    integer(kind=8), intent(in) :: dimdef, dimcon
    real(kind=8), intent(in) :: dt
    real(kind=8), intent(inout) :: congem(dimcon)
    real(kind=8), intent(inout) :: congep(dimcon)
    real(kind=8), intent(out) :: r(dimdef+1)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute stress vector {R}
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  ndim             : dimension of element (2 ou 3)
! In  fnoevo           : .true. if compute in non-linear operator (transient terms)
! In  mecani           : parameters for mechanic
! In  press1           : parameters for hydraulic (first pressure)
! In  press1           : parameters for hydraulic (second pressure)
! In  tempe            : parameters for thermic
! In  second           : parameters for second gradient
! In  dimdef           : number of generalized strains
! In  dimcon           : number of generalized stresses
! In  dt               : time increment
! IO  congem           : generalized stresses at the beginning of time step
!                    => output sqrt(2) on SIG_XY, SIG_XZ, SIG_YZ
! IO  congep           : generalized stresses at the end of time step
!                    => output sqrt(2) on SIG_XY, SIG_XZ, SIG_YZ
! Out r                : stress vector
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbpha1, nbpha2
    integer(kind=8) :: addeme, addete, addep1, addep2, adde2nd
    integer(kind=8) :: adcome, adcote, adcp11, adcp12, adcp21, adcp22, adco2nd
    integer(kind=8) :: i_dim
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: gravity(3)
!
! --------------------------------------------------------------------------------------------------
!
    r(1:dimdef+1) = 0.d0
!
! - Compute gravity
!
    call thmEvalGravity(ds_thm, 0.d0, gravity)
!
! - Get active physics
!
    nbpha1 = press1(2)
    nbpha2 = press2(2)
!
! - Get adresses in generalized vectors
!
    addeme = mecani(2)
    addep1 = press1(3)
    adcome = mecani(3)
    addete = tempe(2)
    adcote = tempe(3)
    addep2 = press2(3)
    adcp11 = press1(4)
    adcp12 = press1(5)
    addep2 = press2(3)
    adcp21 = press2(4)
    adcp22 = press2(5)
    adde2nd = second(2)
    adco2nd = second(3)
!
! - Transforms stress with sqrt(2)
!
    if (ds_thm%ds_elem%l_dof_meca) then
        do i_dim = 4, 6
            congep(adcome+6+i_dim-1) = congep(adcome+6+i_dim-1)*rac2
            congep(adcome+i_dim-1) = congep(adcome+i_dim-1)*rac2
        end do
    end if
!
! - Compute residual {R}
!
    if (ds_thm%ds_elem%l_dof_meca) then
! ----- {R} from mechanic
        do i_dim = 1, 6
            r(addeme+ndim+i_dim-1) = r(addeme+ndim+i_dim-1)+congep(adcome-1+i_dim)
        end do
        do i_dim = 1, 6
            r(addeme+ndim-1+i_dim) = r(addeme+ndim-1+i_dim)+congep(adcome+6+i_dim-1)
        end do
! ----- {R} from hydraulic (first)
        if (ds_thm%ds_elem%l_dof_pre1) then
            do i_dim = 1, ndim
                r(addeme+i_dim-1) = r(addeme+i_dim-1)-gravity(i_dim)*congep(adcp11)
            end do
            if (nbpha1 .gt. 1) then
                do i_dim = 1, ndim
                    r(addeme+i_dim-1) = r(addeme+i_dim-1)-gravity(i_dim)*congep(adcp12)
                end do
            end if
        end if
! ----- {R} from hydraulic (second)
        if (ds_thm%ds_elem%l_dof_pre2) then
            do i_dim = 1, ndim
                r(addeme+i_dim-1) = r(addeme+i_dim-1)-gravity(i_dim)*congep(adcp21)
            end do
            if (nbpha2 .gt. 1) then
                do i_dim = 1, ndim
                    r(addeme+i_dim-1) = r(addeme+i_dim-1)-gravity(i_dim)*congep(adcp22)
                end do
            end if
        end if
    end if
! - For transient terms
    if (fnoevo) then
! ----- {R(t)} from hydraulic (first)
        if (ds_thm%ds_elem%l_dof_pre1) then
            r(addep1) = r(addep1)-congep(adcp11)+congem(adcp11)
            if (nbpha1 .gt. 1) then
                r(addep1) = r(addep1)-congep(adcp12)+congem(adcp12)
            end if
            do i_dim = 1, ndim
                r(addep1+i_dim) = r(addep1+i_dim)+dt*congep(adcp11+i_dim)
            end do
            if (nbpha1 .gt. 1) then
                do i_dim = 1, ndim
                    r(addep1+i_dim) = r(addep1+i_dim)+dt*congep(adcp12+i_dim)
                end do
            end if
            if (ds_thm%ds_elem%l_dof_ther) then
                do i_dim = 1, ndim
                    r(addete) = r(addete)+dt*congep(adcp11+i_dim)*gravity(i_dim)
                end do
                if (nbpha1 .gt. 1) then
                    do i_dim = 1, ndim
                        r(addete) = r(addete)+dt*congep(adcp12+i_dim)*gravity(i_dim)
                    end do
                end if
                do i_dim = 1, ndim
                    r(addete+i_dim) = r(addete+i_dim)+ &
                                      dt*congep(adcp11+ndim+1)*congep(adcp11+i_dim)
                end do
                if (nbpha1 .gt. 1) then
                    do i_dim = 1, ndim
                        r(addete+i_dim) = r(addete+i_dim)+ &
                                          dt*congep(adcp12+ndim+1)*congep(adcp12+i_dim)
                    end do
                end if
!
            end if
        end if
! ----- {R(t)} from hydraulic (second)
        if (ds_thm%ds_elem%l_dof_pre2) then
            r(addep2) = r(addep2)-congep(adcp21)+congem(adcp21)
            if (nbpha2 .gt. 1) then
                r(addep2) = r(addep2)-congep(adcp22)+congem(adcp22)
            end if
            do i_dim = 1, ndim
                r(addep2+i_dim) = r(addep2+i_dim)+dt*congep(adcp21+i_dim)
            end do
            if (nbpha2 .gt. 1) then
                do i_dim = 1, ndim
                    r(addep2+i_dim) = r(addep2+i_dim)+dt*congep(adcp22+i_dim)
                end do
            end if
            if (ds_thm%ds_elem%l_dof_ther) then
                do i_dim = 1, ndim
                    r(addete) = r(addete)+dt*congep(adcp21+i_dim)*gravity(i_dim)
                end do
                if (nbpha2 .gt. 1) then
                    do i_dim = 1, ndim
                        r(addete) = r(addete)+dt*congep(adcp22+i_dim)*gravity(i_dim)
                    end do
                end if
                do i_dim = 1, ndim
                    r(addete+i_dim) = r(addete+i_dim)+ &
                                      dt*congep(adcp21+ndim+1)*congep(adcp21+i_dim)
                end do
                if (nbpha2 .gt. 1) then
                    do i_dim = 1, ndim
                        r(addete+i_dim) = r(addete+i_dim)+ &
                                          dt*congep(adcp22+ndim+1)*congep(adcp22+i_dim)
                    end do
                end if
            end if
        end if
! ----- {R(t)} from thermic
        if (ds_thm%ds_elem%l_dof_ther) then
            r(dimdef+1) = r(dimdef+1)-(congep(adcote)-congem(adcote))
            do i_dim = 1, ndim
                r(addete+i_dim) = r(addete+i_dim)+dt*congep(adcote+i_dim)
            end do
        end if
    end if
! - Second gradient terms
    if (ds_thm%ds_elem%l_dof_2nd) then
        do i_dim = 1, ndim+3
            r(adde2nd-1+i_dim) = r(adde2nd-1+i_dim)+congep(adco2nd-1+i_dim)
        end do
    end if
!
end subroutine
