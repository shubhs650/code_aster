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
! aslint: disable=W1504,W0104,C1505,W1306

subroutine lc5001(BEHInteg, &
                  fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, typmod, ndsde, &
                  dsidep, codret)

    use Behaviour_type
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/ElasticityMaterial_type.h"
! --------------------------------------------------------------------------------------------------
    type(Behaviour_Integ)        :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
    character(len=16), intent(in) :: compor(COMPOR_SIZE), option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    real(kind=8)                 :: sigp(nsig)
    real(kind=8)                 :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out):: codret
! --------------------------------------------------------------------------------------------------
!   RELATIONS ELAS
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: ndimsi, iok(1)
    real(kind=8) :: vi(nvi), carac(1), young_m, young_p, depsth, epsth_p
    real(kind=8) :: sig(BEHInteg%behavPara%ndimsi), eps(BEHInteg%behavPara%ndimsi)
    real(kind=8) :: dsde(BEHInteg%behavPara%ndimsi, BEHInteg%behavPara%ndimsi)
! --------------------------------------------------------------------------------------------------
    ndimsi = BEHInteg%behavPara%ndimsi
    ASSERT(ndimsi .eq. 1)
    ASSERT(BEHInteg%materPara%elasID .eq. ELAS_ISOT)

    sig = 0
    vi = 0
    dsde = 0
    eps = epsm(1:ndimsi)+deps(1:ndimsi)

    if (BEHInteg%behavPara%lVari) vip = 0

    if (compor(INCRELAS) .eq. 'COMP_INCR') then

        call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ELAS', 0, ' ', [0.d0], 1, ["E"], carac, &
                    iok, 1)
        young_m = carac(1)

        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.d0], 1, ["E"], carac, &
                    iok, 1)
        young_p = carac(1)

        call verift(fami, kpg, ksp, 'T', imate, epsth_=depsth)

        sig = young_p*(sigm(1:ndimsi)/young_m+deps-depsth)
        dsde = merge(young_m, young_p, BEHInteg%behavPara%lPred)

    else if (compor(INCRELAS) .eq. 'COMP_ELAS') then

        ! call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.d0], 1, ["E"], carac, &
        !         iok, 2)
        ! young_p = carac(1)

        ! call verift(fami, kpg, ksp, '+', imate, epsth_=epsth_p)

        ! sig  = young_p*(eps - epsth_p)
        ! dsde = merge(young_m, young_p, BEHInteg%behavPara%lPred)

        ! On passe quand même en incrémental pour rester compatible avec un éventuel
        ! chargement de précontrainte (ça disparaîtra en même temps que cette notion
        ! de charge de précontrainte). On émet une alarme
        call utmess('A', 'COMPOR7_3')

        call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ELAS', 0, ' ', [0.d0], 1, ["E"], carac, &
                    iok, 1)
        young_m = carac(1)

        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.d0], 1, ["E"], carac, &
                    iok, 1)
        young_p = carac(1)

        call verift(fami, kpg, ksp, 'T', imate, epsth_=depsth)

        sig = young_p*(sigm(1:ndimsi)/young_m+deps-depsth)
        dsde = merge(young_m, young_p, BEHInteg%behavPara%lPred)

    else
        ASSERT(ASTER_FALSE)
    end if

    codret = 0
    if (BEHInteg%behavPara%lSigm) sigp(1:ndimsi) = sig
    if (BEHInteg%behavPara%lVari) vip(1:nvi) = vi
    if (BEHInteg%behavPara%lMatr) dsidep(1:ndimsi, 1:ndimsi) = dsde

end subroutine
