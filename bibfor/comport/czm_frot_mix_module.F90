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

module czm_frot_mix_module

    implicit none
    private
    public:: CONSTITUTIVE_LAW, Init, Integrate

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"

! --------------------------------------------------------------------------------------------------

    ! Material characteristics
    type MATERIAL
        real(kind=8) :: kn, kt, cohe, frot, k, tau
        aster_logical :: ad, regu_visc
    end type MATERIAL

    ! CZM_FROT_MIX class
    type CONSTITUTIVE_LAW
        integer(kind=8) :: exception = 0
        aster_logical :: elas, rigi, pred
        integer(kind=8) :: ndim
        real(kind=8) :: r
        real(kind=8), dimension(:), allocatable:: phi, deltap, deltav
        real(kind=8) :: dt
        real(kind=8) :: statep
        real(kind=8)  :: cvuser
        type(MATERIAL):: mat
    end type CONSTITUTIVE_LAW

contains

! =====================================================================
!  OBJECT CREATION AND INITIALISATION
! =====================================================================

    function Init(ndim, option, fami, kpg, ksp, imate, t, su, vim, dt, precvg) result(self)

        implicit none

        integer(kind=8), intent(in) :: kpg, ksp, imate, ndim
        real(kind=8), intent(in) :: t(:), su(:), vim(:), dt, precvg
        character(len=16), intent(in) :: option
        character(len=*), intent(in) :: fami
        type(CONSTITUTIVE_LAW) :: self
! --------------------------------------------------------------------------------------------------
! ndim      displacement jump dimension
! option    computation option
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! t         cohesive forces
! su        displacement jump
! vim       internal variables at the beginning of the time step
! dt        time increment
! precvg    precision on stress
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter :: nbel = 2, nbpl = 3, nbvi = 2, nblg = 1
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: iokel(nbel), iokpl(nbpl), iokvi(nbvi), ioklg(nblg)
        real(kind=8) :: valel(nbel), valpl(nbpl), valvi(nbvi), vallg(nblg)
        character(len=16) :: nomel(nbel), nompl(nbpl), nomvi(nbvi), nomlg(nblg)
! --------------------------------------------------------------------------------------------------
        data nomel/'RIGI_NOR', 'RIGI_TAN'/
        data nompl/'ADHE', 'COHESION', 'COEF_FROT'/
        data nomvi/'RIGI_REGU_VISC', 'TAU_REGU_VISC'/
        data nomlg/'PENA_LAGR_ABSO'/
! --------------------------------------------------------------------------------------------------

        ! Dimension controls
        ASSERT(size(t) .eq. ndim)
        ASSERT(size(su) .eq. ndim)

        ! General parameters
        self%ndim = ndim
        self%cvuser = precvg

        ! Options
        self%elas = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'FULL_MECA_ELAS'
        self%rigi = option .eq. 'RIGI_MECA_TANG' .or. option .eq. 'RIGI_MECA_ELAS' &
                    .or. option .eq. 'FULL_MECA' .or. option .eq. 'FULL_MECA_ELAS'
        self%pred = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'RIGI_MECA_TANG'

        ! Elastic parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_FROT_MIX', 0, ' ', [0.d0], nbel, nomel, &
                    valel, iokel, 0)

        ! Plastic parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_FROT_MIX', 0, ' ', [0.d0], nbpl, nompl, &
                    valpl, iokpl, 2)

        ! Viscous regularisation parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_FROT_MIX', 0, ' ', [0.d0], nbvi, nomvi, &
                    valvi, iokvi, 0)

        ! Augmentation coefficient
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_FROT_MIX', 0, ' ', [0.d0], nblg, nomlg, &
                    vallg, ioklg, 2)

        self%mat%kn = valel(1)
        self%mat%kt = valel(2)

        select case (nint(valpl(1)))
        case (0)
            self%mat%ad = ASTER_FALSE
        case (1)
            self%mat%ad = ASTER_TRUE
        end select
        self%mat%cohe = valpl(2)
        self%mat%frot = valpl(3)

        if (iokvi(1) .eq. 0 .and. iokvi(2) .eq. 0) then
            self%mat%regu_visc = ASTER_TRUE
            self%mat%k = valvi(1)
            self%mat%tau = valvi(2)
        else
            self%mat%regu_visc = ASTER_FALSE
        end if

        self%r = vallg(1)

        ! Check that RIGI_NOR and RIGI_TAN are found if ADHE = 'ELAS'
        if ((.not. self%mat%ad) .and. &
            (iokel(1) .eq. 1 .or. iokel(2) .eq. 1)) call utmess('F', 'MECANONLINE3_4')

        ! Constitutive input phi = t + r*su
        allocate (self%phi(ndim))
        self%phi = t+self%r*su

        ! Previous plastic jump
        allocate (self%deltap(ndim-1))
        self%deltap = vim(1:ndim-1)

        ! Previous viscous jump
        allocate (self%deltav(ndim))
        if (self%mat%regu_visc) then
            self%deltav = vim(8:7+ndim)
        else
            self%deltav = 0.d0
        end if

        ! Previous stick/slip state
        self%statep = vim(6)

        ! Time increment
        self%dt = dt

    end function Init

! =====================================================================
!  INTEGRATION OF THE CONSTITUTIVE LAW
! =====================================================================

    subroutine Integrate(self, delta, dphi_delta, vi)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout) :: self
        real(kind=8), intent(out) :: delta(:), dphi_delta(:, :), vi(:)
! --------------------------------------------------------------------------------------------------
! delta         Gauss point jump
! dphi_delta    derivative d(delta)/d(phi)
! vi            internal variables
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: i, j
        real(kind=8) :: tel(self%ndim), nel(self%ndim), id(self%ndim, self%ndim)
        real(kind=8) :: fslip, fopen, teln, telq
        real(kind=8) :: dka, delta_nl
        real(kind=8) :: alpha_n, alpha_t, alpha_v
        integer(kind=8) :: state
!--------------------------------------------------------------------------------------------------

        ! Initialisation
        delta = 0.d0
        dphi_delta = 0.d0
        dka = 0.d0
        tel = 0.d0
        id = 0.d0
        if (self%mat%ad) then
            alpha_n = 1.d0
            alpha_t = 1.d0
        else
            alpha_n = 1.d0/(1.d0+self%r/self%mat%kn)
            alpha_t = 1.d0/(1.d0+self%r/self%mat%kt)
        end if
        if (self%mat%regu_visc) then
            alpha_v = self%mat%k/(1.d0+self%dt/self%mat%tau)
        else
            alpha_v = 0.d0
        end if
        forall (i=1:self%ndim) id(i, i) = 1.d0

! ======================================================================
!  COMPUTATION OF THE NORMAL NONLINEAR JUMP
! ======================================================================

        if (self%mat%frot .gt. 0.d0) then
            fopen = self%phi(1)*alpha_n-self%mat%cohe/self%mat%frot+alpha_v*self%deltav(1)
            delta_nl = max(fopen, 0.d0)/(self%r*alpha_n+alpha_v)
        else
            delta_nl = 0.d0
        end if

! ======================================================================
!  COMPUTATION OF THE TANGENTIAL PLASTIC JUMP
! ======================================================================

        ! Elastic prediction
        tel(1) = (self%phi(1)-self%r*delta_nl)*alpha_n
        do i = 2, self%ndim
            tel(i) = (self%phi(i)-self%r*self%deltap(i-1))*alpha_t
        end do
        teln = tel(1)
        telq = sqrt(dot_product(tel(2:self%ndim), tel(2:self%ndim)))
        fslip = telq+self%mat%frot*(teln-alpha_v*(delta_nl-self%deltav(1)))-self%mat%cohe

        ! Stick or slip
        if (fslip .le. self%cvuser*self%mat%cohe) then
            state = 0
            vi(1:self%ndim-1) = self%deltap(1:self%ndim-1)
        else
            state = 1
            dka = fslip/(self%r*alpha_t+alpha_v)
            do i = 2, self%ndim
                nel(i) = tel(i)/telq
                vi(i-1) = self%deltap(i-1)+dka*nel(i)
            end do
        end if

! ======================================================================
!  COMPUTATION OF THE VISCOUS JUMP
! ======================================================================

        if (self%mat%regu_visc) then
            vi(8) = self%deltav(1)+(delta_nl-self%deltav(1))/(1.d0+self%mat%tau/self%dt)
            do i = 2, self%ndim
                vi(7+i) = self%deltav(i)+(vi(i-1)-self%deltav(i))/(1.d0+self%mat%tau/self%dt)
            end do
        end if

! ======================================================================
!  COMPUTATION OF DELTA
! ======================================================================

        if (.not. self%mat%ad) then
            delta(1) = (self%phi(1)/self%mat%kn+delta_nl)*alpha_n
            do i = 2, self%ndim
                delta(i) = (self%phi(i)/self%mat%kt+vi(i-1))*alpha_t
            end do
        else
            delta(1) = delta_nl
            do i = 2, self%ndim
                delta(i) = vi(i-1)
            end do
        end if

! ======================================================================
!  COMPUTATION OF THE POSTPROCESSING INTERNAL VARIABLES
! ======================================================================

        ! Jump
        vi(3:self%ndim+2) = delta(1:self%ndim)

        ! Stick or slip
        vi(6) = state

        ! Contact or gap
        if (delta(1) .le. 0.d0) then
            vi(7) = 0
        else
            vi(7) = 1
        end if

        ! Viscous stress
        vi(11) = self%mat%k*self%mat%tau &
                 *sqrt(dot_product( &
                       vi(8:7+self%ndim)-self%deltav(1:self%ndim), &
                       vi(8:7+self%ndim)-self%deltav(1:self%ndim)) &
                       )/self%dt

! ======================================================================
!  COMPUTATION OF THE TANGENT MATRIX
! ======================================================================

        if (.not. self%rigi) goto 999

        ! d(delta(1))/d(phi(1))
        if (.not. self%mat%ad) then
            dphi_delta(1, 1) = 1.d0/(self%mat%kn+self%r)
        end if
        if (delta_nl .gt. 0.d0) then
            dphi_delta(1, 1) = dphi_delta(1, 1)+alpha_n**2/(self%r*alpha_n+alpha_v)
        end if

        ! Tangent matrix selection for the prediction
        if (self%pred) then
            state = self%statep
        end if

        if (state .eq. 0 .or. self%elas) then
            ! d(delta(i))/d(phi(i)), for i=2 to ndim
            if (.not. self%mat%ad) then
                do i = 2, self%ndim
                    dphi_delta(i, i) = 1.d0/(self%mat%kt+self%r)
                end do
            end if
        else
            ! d(delta(i))/d(phi(1)), for i=2 to ndim
            if (delta_nl .le. 0.d0) then
                do i = 2, self%ndim
                    dphi_delta(i, 1) = self%mat%frot*nel(i) &
                                       *alpha_t*alpha_n/(self%r*alpha_t+alpha_v)
                end do
            end if
            ! d(delta(i))/d(phi(j)), for i,j=2 to ndim
            do i = 2, self%ndim
                if (.not. self%mat%ad) then
                    dphi_delta(i, i) = 1.d0/(self%mat%kt+self%r)
                end if
                do j = 2, self%ndim
                    dphi_delta(i, j) = dphi_delta(i, j) &
                                       +(nel(i)*nel(j)+fslip/telq*(id(i, j)-nel(i)*nel(j))) &
                                       *alpha_t**2/(self%r*alpha_t+alpha_v)
                end do
            end do
        end if

999     continue

    end subroutine Integrate

end module czm_frot_mix_module
