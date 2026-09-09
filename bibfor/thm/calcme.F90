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
subroutine calcme(BEHInteg, &
                  ndim, option, typmod, &
                  compor, carcri, &
                  instam, instap, &
                  addeme, adcome, dimdef, dimcon, &
                  defgem, deps, &
                  congem, vintm, &
                  congep, vintp, &
                  dsdeme, retcom)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmcomp.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam, instap
    integer(kind=8), intent(in) :: ndim, dimdef, dimcon, addeme, adcome
    real(kind=8), intent(in) :: vintm(*)
    real(kind=8), intent(in) :: defgem(dimdef), deps(6), congem(dimcon)
    real(kind=8), intent(inout) :: congep(dimcon)
    real(kind=8), intent(inout) :: vintp(*)
    real(kind=8), intent(out) :: dsdeme(6, 6)
    integer(kind=8), intent(out) :: retcom
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Standard mechanical behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : option to compute
! In  ndim             : dimension of space (2 or 3)
! In  typmod           : type of modelization (TYPMOD2)
! In  compor           : behaviour
! In  carcri           : parameters for comportment
! In  instam           : time at beginning of time step
! In  instap           : time at end of time step
! In  addeme           : adress of mechanic dof in vector and matrix (generalized quantities)
! In  adcome           : adress of mechanic stress in generalized stresses vector
! In  dimdef           : dimension of generalized strains vector
! In  dimcon           : dimension of generalized stresses vector
! In  defgem           : generalized strains - At begin of current step
! In  deps             : increment of mechanic strains
! In  congem           : generalized stresses - At begin of current step
! In  vintm            : internal state variables - At begin of current step
! IO  congep           : generalized stresses - At end of current step
! IO  vintp            : internal state variables - At end of current step
! Out dsdeme           : derivative of sigma/stress for mechanics
! Out retcom           : return code
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: nsig = 6, neps = 6, ndsdeme = 36
!
! --------------------------------------------------------------------------------------------------
!
    dsdeme = 0.d0
    retcom = 0

! - Integration of mechanical behaviour
    call nmcomp(BEHInteg, &
                ndim, option, typmod, &
                instam, instap, &
                compor, carcri, multComp, &
                neps, defgem(addeme+ndim), deps, nsig, congem(adcome), &
                vintm, &
                congep(adcome), vintp, ndsdeme, dsdeme, retcom)
!
end subroutine
