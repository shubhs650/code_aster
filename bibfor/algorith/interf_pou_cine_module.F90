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

module interf_pou_cine_module

    use Behaviour_type

    use scalar_newton_module, only: &
        newton_state, &
        utnewt

    implicit none
    private
    ! public:: CONSTITUTIVE_LAW, Init, Integrate, PathFollowing
    public:: CONSTITUTIVE_LAW, Init, Integrate

#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcprt2.h"
#include "asterfort/matinv.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"

! --------------------------------------------------------------------------------------------------

    ! Material parameters
    type MATERIAL
        real(kind=8) :: kn
        real(kind=8) :: kt
        real(kind=8) :: fny
        real(kind=8) :: fty
        real(kind=8) :: gn
        real(kind=8) :: gt
    end type MATERIAL

    ! Shared attibutes through the global variable self
    type CONSTITUTIVE_LAW
        integer(kind=8) :: exception = 0
        aster_logical :: rigi, pred
        integer(kind=8) :: ndimsi, itemax
        real(kind=8) :: cvuser
        real(kind=8)  :: deltat
        type(MATERIAL):: mat
        real(kind=8), pointer :: rot(:) => null()
    end type CONSTITUTIVE_LAW

contains

! ==================================================================================================
!  OBJECT CREATION AND INITIALISATION
! ==================================================================================================

    function Init(ndimsi, option, fami, kpg, ksp, imate, nomat, deltat, &
                  itemax, precvg, BEHinteg_) result(self)

        implicit none

        integer(kind=8), intent(in) :: kpg, ksp, imate, itemax, ndimsi
        real(kind=8), intent(in)    :: precvg, deltat
        character(len=16), intent(in) :: option
        character(len=*), intent(in) :: fami
        character(len=8), intent(in) :: nomat
        type(Behaviour_Integ), optional:: BEHinteg_
        type(CONSTITUTIVE_LAW) :: self
! --------------------------------------------------------------------------------------------------
! ndimsi    symmetric tensor dimension (2*ndim)
! option    computation option
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! deltat    time increment (instap - instam)
! itemax    max number of iterations for the solver
! precvg    required accuracy (with respect to stress level))
! BEHinteg_ behaviour
! --------------------------------------------------------------------------------------------------
        self%rigi = option .eq. 'RIGI_MECA_TANG' .or. option .eq. 'RIGI_MECA_ELAS' &
                    .or. option .eq. 'FULL_MECA' .or. option .eq. 'FULL_MECA_ELAS'
        self%pred = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'RIGI_MECA_TANG'
        self%ndimsi = ndimsi
        self%itemax = itemax
        self%cvuser = precvg
        self%deltat = deltat
        self%mat = GetMaterial(self, fami, kpg, ksp, imate, nomat, BEHinteg_=BEHinteg_)

        AS_ALLOCATE(vr=self%rot, size=ndimsi*ndimsi)
        self%rot = GetMatrot(self, BEHinteg_)

    end function Init

! =====================================================================
!  INTEGRATION OF THE CONSTITUTIVE LAW (MAIN ROUTINE)
! =====================================================================

    subroutine Integrate(self, eps, vim, sig, vip, dsde)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)    :: eps(:), vim(:)
        real(kind=8), intent(out)    :: sig(:), vip(:), dsde(:, :)
! --------------------------------------------------------------------------------------------------
! eps   current strain
! vim   internal variables at the beginning of the time step
! sig   stress at the end of the time step
! vip   internal variables at the end of the time step
! dsde  tangent matrix (ndimsi,ndimsi)
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: state
        real(kind=8) :: epm(self%ndimsi), ep(self%ndimsi)
        real(kind=8) :: invrot(self%ndimsi, self%ndimsi), r8bid
! --------------------------------------------------------------------------------------------------

!   Unpack internal variables
        call pmavec('ZERO', self%ndimsi, self%rot, vim(1:self%ndimsi), epm)

!   Behaviour integration
        call ComputePlasticity(self, eps, epm, state, ep, sig, dsde)
        if (self%exception .ne. 0) goto 999

!   Pack internal variables
        call matinv('S', self%ndimsi, self%rot, invrot, r8bid)
        call pmavec('ZERO', self%ndimsi, invrot, ep, vip(1:self%ndimsi))
        AS_DEALLOCATE(vr=self%rot)

999     continue
    end subroutine Integrate

! =====================================================================
!  PLASTICITY COMPUTATION AND TANGENT OPERATORS
! =====================================================================

    subroutine ComputePlasticity(self, eps, epm, state, ep, sig, dsde)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in) :: eps(:)
        real(kind=8), intent(in) :: epm(:)
        integer(kind=8), intent(out) :: state
        real(kind=8), intent(out) :: ep(:)
        real(kind=8), intent(out) :: sig(:)
        real(kind=8), intent(out) :: dsde(:, :)
! --------------------------------------------------------------------------------------------------
! eps   current strain
! epm   plastic strain at the beginning of the time step
! state 0=elastic, 1=plastic
! ep    plastic strain at the end of the time step
! sig   stress at the end of the time step
! dsde  tangent matrix (ndimsi,ndimsi)
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: ite
        real(kind=8) :: dka
        real(kind=8) :: fel, tel(size(eps)), ksiel(size(eps))
        real(kind=8) :: dep(size(eps))
        real(kind=8) :: cel(size(eps)), gl(size(eps))
        real(kind=8) :: sigYinv(size(eps)), Al(size(eps)), Bl(size(eps))
        real(kind=8) :: Dl(size(eps))
        real(kind=8) :: cpl(size(eps), size(eps))
        real(kind=8) :: dq_dsig(size(eps)), dq_dr(size(eps))
        real(kind=8) :: res, dres
        type(newton_state):: mem
! --------------------------------------------------------------------------------------------------

!   Initialization
        cel = (/self%mat%kt, self%mat%kn, self%mat%kn/)
        gl = (/self%mat%gt, self%mat%gn, self%mat%gn/)
        sigYinv = (/1.d0/self%mat%fty, 1.d0/self%mat%fny, 1.d0/self%mat%fny/)
        Bl = 2.d0*sigYinv**2*cel
        Dl = 2.d0*sigYinv**2*gl
        Al = Bl+Dl

!   Elastic prediction
        tel = cel*(eps-epm)
        ksiel = sigYinv*(tel-gl*epm)
        fel = q_dka(self, tel, ksiel, epm, sigYinv, gl, Al, 0.d0)

! ==============================================
!               STRESS COMPUTATION
! ==============================================

!   Elastic regime
        if (fel .le. self%cvuser) then

            state = 0
            dka = 0.d0
            dep = 0.d0

        end if

!   Plastic regime
        if (fel .gt. self%cvuser) then

            ! Initialialization
            dka = 0.d0

            do ite = 1, self%itemax
                ! Residual
                res = -q_dka(self, tel, ksiel, epm, sigYinv, gl, Al, dka)
                ! Convergence check
                if (abs(res) .le. self%cvuser) exit
                ! Iteration
                dres = -dq_dka(self, tel, ksiel, epm, sigYinv, gl, Al, dka)
                dka = utnewt(dka, res, dres, ite, mem, xmin=0.d0)
            end do

            if (ite .gt. self%itemax) then
                self%exception = 1
                goto 999
            end if

            state = 1
            dep = dka*2*sigYinv*ksiel/(1.d0+Al*dka)

        end if

        ep = epm+dep
        sig = tel-cel*dep

! ==========================================
!               TANGENT MATRIX
! ==========================================

        if (.not. self%rigi) goto 999

!   Special treatment for the prediction
        if (self%pred) then

            if (fel .lt. -self%cvuser) then
                state = 0
            else
                state = 1
            end if

        end if

!   Elastic regime
        if (state .eq. 0) then

            dsde = diag(cel)

        end if

!   Plastic regime
        if (state .eq. 1) then

            ! Derivative of q(sig, r) with respect to sig
            dq_dsig = 2*sigYinv**2*(sig-gl*ep)

            ! Derivative of q(sig, r) with respect to r
            dq_dr = -dq_dsig

            ! Tangent operator
            call lcprt2(self%ndimsi, cel*dq_dsig, dq_dsig*cel, cpl)
            dsde = diag(cel)-(cpl/(dot_product(dq_dsig, cel*dq_dsig) &
                                   -dot_product(dq_dr, gl*dq_dsig)))

        end if

999     continue

    end subroutine ComputePlasticity

! ==================================================================================================
!  MATERIAL CHARACTERISTICS
! ==================================================================================================

    function GetMaterial(self, fami, kpg, ksp, imate, nomat, BEHinteg_) result(mat)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout) :: self
        integer(kind=8), intent(in) :: kpg, ksp, imate
        character(len=*), intent(in) :: fami
        character(len=8), intent(in) :: nomat
        type(Behaviour_Integ), optional:: BEHinteg_
        type(MATERIAL) :: mat
! --------------------------------------------------------------------------------------------------
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! nomat     material name
! BEHinteg_ behaviour
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter:: nbel = 2, nbpl = 4, nbco = 3
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: iokel(nbel), iokpl(nbpl), nb_para
        real(kind=8) :: valel(nbel), valpl(nbpl), para_vale(nbco)
        character(len=16) :: nomel(nbel), nompl(nbpl)
        character(len=8) :: para_name(nbco)
        character(len=1) :: poum
! --------------------------------------------------------------------------------------------------
        data nomel/'K_N', 'K_T'/
        data nompl/'F_NY', 'F_TY', 'G_N', 'G_T'/
! --------------------------------------------------------------------------------------------------

        if (self%pred) then
            poum = '-'
        else
            poum = '+'
        end if

        nb_para = 0
        para_name = ' '
        para_vale = 0.d0

!   Parameter evolution
        if (present(BEHinteg_)) then
            nb_para = nb_para+1
            para_name(nb_para) = 'X'
            para_vale(nb_para) = BEHinteg_%behavESVA%behavESVAGeom%coorElga(kpg, 1)
            nb_para = nb_para+1
            para_name(nb_para) = 'Y'
            para_vale(nb_para) = BEHinteg_%behavESVA%behavESVAGeom%coorElga(kpg, 2)
            nb_para = nb_para+1
            para_name(nb_para) = 'Z'
            para_vale(nb_para) = BEHinteg_%behavESVA%behavESVAGeom%coorElga(kpg, 3)
        end if

!   Elasticity
        call rcvalb(fami, kpg, ksp, '+', imate, nomat, &
                    'INTERF_POU_ELAS', nb_para, para_name, [para_vale], &
                    nbel, nomel, valel, iokel, 2)

!   Plasticity
        call rcvalb(fami, kpg, ksp, '+', imate, nomat, &
                    'INTERF_POU_CINE', nb_para, para_name, [para_vale], &
                    nbpl, nompl, valpl, iokpl, 2)

        mat%kn = valel(1)
        mat%kt = valel(2)
        mat%fny = valpl(1)
        mat%fty = valpl(2)
        mat%gn = valpl(3)
        mat%gt = valpl(4)

    end function GetMaterial

! ==================================================================================================
!  ORIENTATION
! ==================================================================================================

    function GetMatrot(self, BEHinteg) result(rot)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout) :: self
        type(Behaviour_Integ), intent(in) :: BEHinteg
        real(kind=8) :: rot(self%ndimsi*self%ndimsi)
! --------------------------------------------------------------------------------------------------
! BEHinteg_ behaviour
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: i
! --------------------------------------------------------------------------------------------------

!   Rotation matrix : global to local
        rot = 0.d0
        do i = 1, self%ndimsi*self%ndimsi
            rot(i) = BEHinteg%behavESVA%behavESVAOther%rotpg(i)
        end do

    end function GetMatrot

! ----------------------------------------------------------------------------------------
!  Useful intermediate functions
! ----------------------------------------------------------------------------------------

    function q_dka(self, tel, ksiel, ep, sigYinv, gl, Al, dka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in) :: self
        real(kind=8) :: res
        real(kind=8), intent(in) :: dka
        real(kind=8), intent(in) :: tel(:), ksiel(:), ep(:)
        real(kind=8), intent(in) :: sigYinv(:), gl(:), Al(:)
        res = sum((ksiel/(1.d0+Al*dka))**2)-1.d0
    end function q_dka

    function dq_dka(self, tel, ksiel, ep, sigYinv, gl, Al, dka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in) :: self
        real(kind=8) :: res
        real(kind=8), intent(in) :: dka
        real(kind=8), intent(in) :: tel(:), ksiel(:), ep(:)
        real(kind=8), intent(in) :: sigYinv(:), gl(:), Al(:)
        res = sum(-2.d0*Al*ksiel**2/(1.d0+Al*dka)**3)
    end function dq_dka

    function diag(v) result(D)
        real(kind=8), intent(in) :: v(:)
        real(kind=8) :: D(size(v), size(v))
        integer(kind=8) :: i
        D(:, :) = 0.d0
        do i = 1, size(v)
            D(i, i) = v(i)
        end do
    end function diag

end module interf_pou_cine_module
