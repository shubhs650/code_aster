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
! aslint: disable=W1306,W1504,C1505
!
subroutine redece(BEHInteg, &
                  ndim, option, typmod, &
                  instam, instap, &
                  compor, carcri, multComp, &
                  neps, epsm, deps, &
                  nsig, sigm, &
                  nvi, vim, &
                  sigp, vip, &
                  ndsde, dsidep, codret, &
                  l_epsi_varc, numlc)
!
    use calcul_module, only: ca_iredec_, ca_td1_, ca_tf1_, ca_timed1_, ca_timef1_
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lc0000.h"
#include "asterfort/utmess.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8), intent(in) :: ndim
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(in) :: instam, instap
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: multComp
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps), deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    real(kind=8), intent(inout) :: sigp(nsig), vip(nvi)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(inout) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                          merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(inout) :: codret
    aster_logical, intent(in) :: l_epsi_varc
    integer(kind=8), intent(in)::  numlc
!
! --------------------------------------------------------------------------------------------------
!
! Mechanical non-linear behaviours
!
! Management of local time step division
!
! --------------------------------------------------------------------------------------------------
!
! IO  BEHInteg         : parameters for integration of behaviour
! In  option           : option to compute
! In  typmod           : type of modeling (3D, 2D, etc.)
! In  jvMaterCode      : adress for material parameters
! In  instam           : time at beginning of current time step
! In  instap           : time at end of current time step
! In  compor           : description of behaviour
! In  carcri           : parameters for integration of behaviour
! In  multComp         : name of JEVEUX object for multi-behaviour (DEFI_COMPOR)
! In  neps             : size of strain tensor
! In  epsm_inp         : strain tensor at beginning of current time step
! In  deps_inp         : increment of strain tensor from beginning of current time step
! In  nsig             : size of stress tensor
! In  sigm             : stress tensor at beginning of current time step
! In  vim              : internal state variables at beginning of current time step
! IO  sigp             : stress tensor at end of current time step
! IO  vip              : internal state variables at end of current time step
! In  ndsde            : size of jacobian matrix (dSig/dEps)
! IO  dsidep           : jacobian matrix (dSig/dEps)
! IO  codret           : return code from integration of behaviour
!     LDC_ERROR_NONE => No problem
!     LDC_ERROR_NCVG => convergence default
!     LDC_ERROR_QUAL => quality problem
!     LDC_ERROR_CPLA => stress plane algorithm not converged
!     LDC_ERROR_DVAL => out of bound for validity
! In  l_epsi_varc      : flag to compute non-mechanical strains (from external state variables)
! In  num_lc           : index of behaviour to integrate
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter:: NBR_DECOUP_MAX = 5
    aster_logical:: lMatrPred, lMatr, lSigm, lVari
    integer(kind=8) :: cutStrategy, cutLevel, niv_ini, npas, iterIntePas, codret_sub, pas
    real(kind=8) :: epsm_sub(neps), deps_sub(neps), sigm_sub(nsig), vim_sub(nvi)
    real(kind=8) :: deltat, tm, tp
    real(kind=8) :: dsidep_sub(merge(nsig, 6, nsig*neps .eq. ndsde), &
                               merge(neps, 6, nsig*neps .eq. ndsde))
    character(len=16) :: defoComp
!
! --------------------------------------------------------------------------------------------------
!
    iterIntePas = nint(carcri(ITER_INTE_PAS))
    defoComp = compor(DEFO)

! - Option (operators) to compute
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)
    lVari = L_VARI(option)
    lMatrPred = L_MATR_PRED(option)

! - To manage external state variable
    ca_iredec_ = 1
    ca_timed1_ = instam
    ca_timef1_ = instap
    ca_td1_ = instam
    ca_tf1_ = instap

! - Select strategy of local time step division
    if (abs(iterIntePas) .le. 1 .or. lMatrPred) then
        cutStrategy = LDC_TIMEDIV_NONE
    else if (iterIntePas .le. -2) then
        cutStrategy = LDC_TIMEDIV_AUTO
    else if (iterIntePas .ge. 2) then
        cutStrategy = LDC_TIMEDIV_FORC
    end if

! --------------------------------------------------------------------------------------------------
!  Integration du comportement sans redecoupage
! --------------------------------------------------------------------------------------------------
    if (cutStrategy .eq. LDC_TIMEDIV_NONE) then
        codret = LDC_ERROR_NONE
        BEHInteg%behavPara%cutLevel = LDC_TIMEDIV_AUTO
        call lc0000(BEHInteg, &
                    ndim, option, typmod, &
                    instam, instap, &
                    compor, carcri, multComp, &
                    neps, epsm, deps, &
                    nsig, sigm, &
                    nvi, vim, &
                    sigp, vip, &
                    ndsde, dsidep, codret, &
                    l_epsi_varc, numlc)
        goto 999
    end if

! --------------------------------------------------------------------------------------------------
!  Integration du comportement avec redecoupage
! --------------------------------------------------------------------------------------------------
! - Some checks
    if (typmod(2) .eq. 'GRADVARI') then
        call utmess('F', 'COMPOR2_10', sk=typmod(2))
    end if
    if (numlc .ge. 8000 .and. numlc .lt. 9000) then
        call utmess('F', 'COMPOR2_10', sk='KIT_DDI')
    end if
    if (defoComp .eq. 'SIMO_MIEHE') then
        call utmess('F', 'COMPOR2_10', sk='SIMO_MIEHE')
    end if
    ASSERT(lSigm)
    ASSERT(lVari)

! - Initial level of local time step division
    niv_ini = merge(0, 1, cutStrategy .eq. LDC_TIMEDIV_AUTO)

! - Boucle sur les niveaux de decoupage
    do cutLevel = niv_ini, NBR_DECOUP_MAX
        codret = LDC_ERROR_NONE

! ----- Number of divisions
        npas = max(1, merge(1, abs(iterIntePas)*(2**(cutLevel-1)), cutLevel .eq. 0))
        deltat = (instap-instam)/npas

! ----- Prepare input parameters
        deps_sub = deps/npas
        vim_sub = vim
        sigm_sub = sigm
        if (lMatr) then
            dsidep = 0
        end if

! ----- Boucle sur les sous-pas
        do pas = 1, npas
! --------- Prepare input parameters
            epsm_sub = epsm+(pas-1)*deps_sub
            tm = instam+(pas-1)*deltat
            tp = instam+pas*deltat
            ca_td1_ = tm
            ca_tf1_ = tp
            codret_sub = LDC_ERROR_NONE
            vip = vim_sub
            BEHInteg%behavPara%cutLevel = cutLevel

! --------- Main switch to the integration of behaviour laws
            call lc0000(BEHInteg, &
                        ndim, option, typmod, &
                        tm, tp, &
                        compor, carcri, multComp, &
                        neps, epsm_sub, deps_sub, &
                        nsig, sigm_sub, &
                        nvi, vim_sub, &
                        sigp, vip, &
                        ndsde, dsidep_sub, codret_sub, &
                        l_epsi_varc, numlc)

            select case (codret_sub)
            case (LDC_ERROR_NONE)
                continue
            case (LDC_ERROR_NCVG)
                codret = LDC_ERROR_NCVG
                exit
            case (LDC_ERROR_QUAL)
                codret = LDC_ERROR_QUAL
                if (pas .ne. npas) exit
            case default
                ASSERT(ASTER_FALSE)
            end select

! --------- Jacobian matrix is an average of Jacobian matrices at each sub-step
            if (lMatr) then
                dsidep = dsidep+dsidep_sub/npas
            end if

! --------- Update input parameters
            sigm_sub = sigp
            vim_sub = vip

        end do
        if (codret .eq. LDC_ERROR_NONE) then
            goto 999
        end if

    end do

999 continue
    ASSERT(codret .ge. LDC_ERROR_NONE .and. codret .le. LDC_ERROR_QUAL)
end subroutine
