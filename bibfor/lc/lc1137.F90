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
! aslint: disable=W1504
!
subroutine lc1137(BEHInteg, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, multComp, carcri, instam, instap, &
                  neps, epsm, deps, sigm, nvi, vim, option, &
                  sigp, vip, &
                  typmod, &
                  dsidep, codret)
!
    use Behaviour_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmvprk.h"
#include "asterfort/plasti.h"
#include "asterfort/utlcal.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: jvMaterCode, nvi
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    character(len=16), intent(in) :: multComp
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam, instap
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
    real(kind=8), intent(in) :: sigm(6)
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(6)
    real(kind=8), intent(out) :: vip(nvi)
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour - Special SIMO_MIEHE
!
! polycristal, monocristal
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHInteg       : parameters for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: algoInte, relaComp
    character(len=11) :: meting
    common/meti/meting
!
! --------------------------------------------------------------------------------------------------
!
    relaComp = COMPOR(RELA_NAME)
    if (relaComp .eq. 'POLYCRISTAL') then
        call nmvprk(BEHInteg, &
                    option, typmod, ndim, &
                    compor, carcri, &
                    instam, instap, &
                    neps, epsm, deps, sigm, nvi, vim, &
                    sigp, vip, dsidep, &
                    codret, multComp)

    elseif (relaComp .eq. 'MONOCRISTAL') then
        call utlcal('VALE_NOM', algoInte, carcri(6))
        if (algoInte(1:6) .eq. 'NEWTON') then
            meting = algoInte(1:11)
            call plasti(BEHInteg, &
                        option, typmod, &
                        fami, kpg, ksp, jvMaterCode, &
                        compor, carcri, instam, instap, &
                        epsm, deps, &
                        sigm, &
                        nvi, vim, &
                        sigp, vip, &
                        dsidep, codret, &
                        multComp)

        else if (algoInte .eq. 'RUNGE_KUTTA') then
            meting = 'RUNGE_KUTTA'
            call nmvprk(BEHInteg, &
                        option, typmod, ndim, &
                        compor, carcri, &
                        instam, instap, &
                        neps, epsm, deps, sigm, nvi, vim, &
                        sigp, vip, dsidep, &
                        codret, multComp)
        else
            write (6, *) 'ALGOInte:', algoInte
            ASSERT(ASTER_FALSE)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
