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
! aslint: disable=W1504,W0104
!
subroutine lc0031(BEHInteg, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, instam, instap, neps, &
                  epsm, deps, sigm, nvi, vim, option, &
                  sigp, vip, typmod, &
                  dsidep, codret)
!
    use Behaviour_type
    use MaterialPara_type
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmveei.h"
#include "asterfort/nmvprk.h"
#include "asterfort/utlcal.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    integer(kind=8) :: jvMaterCode, ndim, kpg, ksp, codret, nvi, neps
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(6), deps(6), sigm(6), sigp(6), vim(nvi), vip(nvi)
    real(kind=8) :: dsidep(6, 6)
    character(len=16) :: option
    character(len=8) :: typmod(*)
    character(len=*) :: fami
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! VENDOCHAB / VISC_ENDO_LEMA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: algoInte
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = BEHInteg%materPara
    call utlcal('VALE_NOM', algoInte, carcri(6))
    if (algoInte .eq. 'RUNGE_KUTTA') then
        call nmvprk(BEHInteg, &
                    option, typmod, ndim, &
                    compor, carcri, &
                    instam, instap, &
                    neps, epsm, deps, sigm, nvi, vim, &
                    sigp, vip, dsidep, &
                    codret)
    else
        call nmveei(materPara, &
                    carcri, compor, ndim, typmod, &
                    instam, instap, &
                    epsm, deps, sigm, nvi, vim, option, &
                    sigp, vip, dsidep, codret)
    end if
!
end subroutine
