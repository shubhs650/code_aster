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
module HHO_compor_module
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
    private
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - mechanics - Behaviour
!
! --------------------------------------------------------------------------------------------------
!
    public :: HHO_Compor_State
    public :: isLargeStrain
    private :: initialize_compor, typmod1
!
    type HHO_Compor_State
!
        aster_logical :: l_debug = ASTER_FALSE
        aster_logical :: lNonLinearOption = ASTER_FALSE
        character(len=16) :: option = " "
        character(len=8) :: typmod(2) = [" ", " "]
        character(len=16) :: multComp = " "
        aster_logical :: l_largestrain = ASTER_FALSE
        aster_logical :: c_plan = ASTER_FALSE
        aster_logical :: axis = ASTER_FALSE
        aster_logical :: matsym = ASTER_FALSE
        integer(kind=8) :: nbsigm = 0
        integer(kind=8) :: lgpg = 0
        integer(kind=8) :: codret = 0
        character(len=16), pointer :: compor(:) => null()
        real(kind=8), pointer :: carcri(:) => null()
        real(kind=8), pointer :: vari_prev(:) => null()
        real(kind=8), pointer :: sig_prev(:) => null()
        real(kind=8), pointer :: vari_curr(:) => null()
        real(kind=8), pointer :: sig_curr(:) => null()

! ----- Behaviour parameters
        type(Behaviour_Integ) :: BEHInteg

! ----- member function
    contains
        procedure, pass :: initialize => initialize_compor
!
    end type HHO_Compor_State
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_compor(this, fami, option, ndim, bary, &
                                 typmod2)
!
        implicit none
!
        class(HHO_Compor_State), intent(inout)  :: this
        character(len=4), intent(in) :: fami
        character(len=16), intent(in) :: option
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: bary(3)
        character(len=8), intent(in) :: typmod2
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a HHO_Compor_State type
!   In this     : a HHo Compor
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: timeCurr, timePrev
        integer(kind=8) :: iret, jvMaterc, jtab(7)
        character(len=16), pointer :: mulcom(:) => null()
        character(len=8) :: typmod(2)
        aster_logical :: lNonLinearOption
        type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
        typmod(1) = typmod1(ndim)
        typmod(2) = typmod2
        this%option = option
        this%typmod = typmod
        this%axis = this%typmod(1) .eq. 'AXIS'
        this%c_plan = this%typmod(1) .eq. 'C_PLAN'
        this%nbsigm = nbsigm()

! ----- Non-linear option ?
        lNonLinearOption = this%option .ne. "RIGI_MECA" .and. &
                           this%option .ne. "FORC_NODA" .and. &
                           this%option .ne. "REFE_FORC_NODA"

! ----- Get material parameters
        call jevech('PMATERC', 'L', jvMaterc)
        call initParaCell(fami, zi(jvMaterc), materPara)
        call getUserLCSWithBaryCenter(ndim, bary, materPara%lcsPara)

! ----- Get time
        timePrev = r8vide()
        timeCurr = r8vide()
        call tecach('NNO', 'PINSTMR', 'L', iret, nval=7, itab=jtab)
        if (iret .eq. 0) then
            timePrev = zr(jtab(1)-1+1)
        end if
        call tecach('NNO', 'PINSTPR', 'L', iret, nval=7, itab=jtab)
        if (iret .eq. 0) then
            timeCurr = zr(jtab(1)-1+1)
        end if

        if (lNonLinearOption) then
! --------- Prepare behaviour
            call behaviourInit(this%BEHInteg)
            call jevech('PCOMPOR', 'L', vk16=this%compor)
            call jevech('PCARCRI', 'L', vr=this%carcri)
            if (this%option(1:4) .ne. "PILO") then
                call jevech('PMULCOM', 'L', vk16=mulcom)
                this%multComp = mulcom(1)
            end if
            this%l_largestrain = isLargeStrain(this%compor(DEFO))
            this%matsym = nint(this%carcri(CARCRI_MATRSYME)) .le. 0

! --------- Set main parameters for behaviour (on cell)
            call behaviourSetParaCell(typmod, option, &
                                      this%compor, this%carcri, &
                                      timePrev, timeCurr, &
                                      materPara, this%BEHInteg)

! --------- Other fields
            call jevech('PCONTMR', 'L', vr=this%sig_prev)
            call jevech('PVARIMR', 'L', vr=this%vari_prev)
            call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
            ASSERT(iret .eq. 0)
            this%lgpg = max(jtab(6), 1)*jtab(7)

        else
            this%matsym = ASTER_TRUE
            this%l_largestrain = ASTER_FALSE
            if (this%option == "FORC_NODA") then
                call jevech('PSIEFR', 'L', vr=this%sig_prev)
                call jevech('PCOMPOR', 'L', vk16=this%compor)
                this%l_largestrain = isLargeStrain(this%compor(DEFO))
            end if
            call behaviourInit(this%BEHInteg)
            this%BEHInteg%materPara = materPara

        end if

        if (L_SIGM(option)) then
            call jevech('PCONTPR', 'E', vr=this%sig_curr)
        end if
        if (L_VARI(option)) then
            call jevech('PVARIPR', 'E', vr=this%vari_curr)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function isLargeStrain(defo_comp) result(bool)
!
        implicit none
!
        character(len=16), intent(in) :: defo_comp
        aster_logical :: bool
!
! --------------------------------------------------------------------------------------------------
!
!   To know if the model is in large deformations
!   In defo_comp : type of deformation
! --------------------------------------------------------------------------------------------------
!
        bool = ASTER_FALSE
!
        if (defo_comp(1:5) .eq. 'PETIT') then
            bool = ASTER_FALSE
        elseif (defo_comp .eq. 'GDEF_LOG') then
            bool = ASTER_TRUE
        elseif (defo_comp .eq. 'GREEN_LAGRANGE') then
            bool = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function typmod1(ndim) result(typmod)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        character(len=8) :: typmod
!
! --------------------------------------------------------------------------------------------------
!
!   Select typmod 1
! --------------------------------------------------------------------------------------------------
!
        select case (ndim)
        case (3)
            typmod = '3D'
        case (2)
            if (lteatt('AXIS', 'OUI')) then
                typmod = 'AXIS'
            else if (lteatt('C_PLAN', 'OUI')) then
                ASSERT(ASTER_FALSE)
                typmod = 'C_PLAN'
            else if (lteatt('D_PLAN', 'OUI')) then
                typmod = 'D_PLAN'
            else
                ASSERT(ASTER_FALSE)
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
    end function
!
end module
