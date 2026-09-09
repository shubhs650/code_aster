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
subroutine pmf_vmis(for_pmf, nf, nbvalc, &
                    pmfCompor, materPara, &
                    varim, contm, ddefp, modf, &
                    sigf, varip, codret)
!
    use pmfcom_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nm1dci.h"
#include "asterfort/nm1dis.h"
#include "asterfort/paeldt.h"
#include "asterfort/rcexistvarc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/vmci1d.h"
#include "MultiFiber_type.h"
!
    type(pmfcom_user), intent(in) :: for_pmf
    integer(kind=8) :: nf, nbvalc
    character(len=24) :: pmfCompor(*)
    type(Material_Para), intent(inout) :: materPara
    real(kind=8) :: varim(nbvalc*nf), contm(nf), ddefp(nf), modf(nf)
    real(kind=8) :: sigf(nf), varip(nbvalc*nf)
    integer(kind=8) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!               COMPORTEMENT DES ÉLÉMENTS DE POUTRE MULTI-FIBRES
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       pmfCompor  : information sur le comportement du groupe de fibres
!       carcri    : critères de convergence locaux
!       nf      : nombre de fibres du groupe
!       nbvalc  : nombre de variable internes
!       defam   : déformations anélastiques a l'instant précédent
!       defap   : déformations anélastiques a l'instant du calcul
!       varim   : variables internes moins
!       varimp  : variables internes itération précédente (pour DE BORST)
!       contm   : contraintes moins par fibre
!       defm    : déformation à l'instant du calcul précédent
!       ddefp   : incrément de déformation
!
!   OUT
!       modf    : module tangent des fibres
!       sigf    : contrainte a l'instant actuel des fibres
!       varip   : variables internes a l'instant actuel
!       codret  : code retour (0 c'est ok)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: nbProp = 2
    integer(kind=8)  :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: ksp, fib, ivari, nbvari_grfibre
    real(kind=8) :: ep, em, depsth
    real(kind=8) :: depsm, nu
    character(len=8)  :: materPoin
    character(len=16) :: relaComp, algoInte
    aster_logical :: istemp
    integer(kind=8) :: kpg, debsp, jvMaterCode
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm
    character(len=16) :: option
!
! --------------------------------------------------------------------------------------------------
!
    kpg = for_pmf%kpg
    jvMaterCode = for_pmf%icdmat
    option = for_pmf%option
    debsp = for_pmf%debsp
    instam = for_pmf%instam
    instap = for_pmf%instap
    epsm = for_pmf%epsm
    codret = 0
    materPoin = pmfCompor(MULTI_FIBER_MATER) (1:8)
    relaComp = pmfCompor(MULTI_FIBER_RELA) (1:16)
    algoInte = pmfCompor(MULTI_FIBER_ALGO) (1:16)
    read (pmfCompor(MULTI_FIBER_NBVARI), '(I24)') nbvari_grfibre
!   Vérification du nombre de fibre
    ASSERT(nbvari_grfibre .le. nbvalc)

    if (relaComp .eq. 'VMIS_CINE_GC') then
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(fami, 1, 1, '+', jvMaterCode, &
                        materPoin, 'ELAS', &
                        0, '', [0.d0], &
                        1, ['E'], propVale, &
                        propCode, 1)
            ep = propVale(1)
            em = ep
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                depsm = ddefp(fib)
                call initParaPoin(kpg, ksp, materPara)
                call vmci1d(materPara, &
                            option, materPoin, &
                            em, ep, &
                            contm(fib), depsm, varim(ivari), &
                            sigf(fib), varip(ivari), modf(fib))
            end do
        else
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', jvMaterCode, materPoin, em, ep, nu, depsth)
                depsm = ddefp(fib)-depsth
                call initParaPoin(kpg, ksp, materPara)
                call vmci1d(materPara, &
                            option, materPoin, &
                            em, ep, &
                            contm(fib), depsm, varim(ivari), &
                            sigf(fib), varip(ivari), modf(fib))
            end do
        end if

    else if (relaComp .eq. 'VMIS_CINE_LINE') then
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(fami, 1, 1, '+', jvMaterCode, materPoin, 'ELAS', &
                        0, '', [0.d0], 1, ['E'], propVale, propCode, 1)
            ep = propVale(1)
            em = ep
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                depsm = ddefp(fib)
                call initParaPoin(kpg, ksp, materPara)
                call nm1dci(materPara, &
                            option, materPoin, &
                            em, ep, &
                            contm(fib), depsm, varim(ivari), &
                            sigf(fib), varip(ivari), modf(fib))
            end do
        else
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', jvMaterCode, materPoin, em, ep, nu, depsth)
                depsm = ddefp(fib)-depsth
                call initParaPoin(kpg, ksp, materPara)
                call nm1dci(materPara, &
                            option, materPoin, &
                            em, ep, &
                            contm(fib), depsm, varim(ivari), &
                            sigf(fib), varip(ivari), modf(fib))
            end do
        end if

    else if ((relaComp .eq. 'VMIS_ISOT_LINE') .or. &
             (relaComp .eq. 'VMIS_ISOT_TRAC')) then
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(fami, 1, 1, '+', jvMaterCode, materPoin, 'ELAS', &
                        0, '', [0.d0], 1, ['E'], propVale, propCode, 1)
            ep = propVale(1)
            em = ep
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                depsm = ddefp(fib)
                call initParaPoin(kpg, ksp, materPara)
                call nm1dis(materPara, &
                            option, relaComp, materPoin, &
                            em, ep, &
                            contm(fib), depsm, varim(ivari), &
                            sigf(fib), varip(ivari), modf(fib))
            end do
        else
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', jvMaterCode, materPoin, em, ep, nu, depsth)
                depsm = ddefp(fib)-depsth
                call initParaPoin(kpg, ksp, materPara)
                call nm1dis(materPara, &
                            option, relaComp, materPoin, &
                            em, ep, &
                            contm(fib), depsm, varim(ivari), &
                            sigf(fib), varip(ivari), modf(fib))
            end do
        end if
    else
        call utmess('F', 'ELEMENTS2_39', sk=relaComp)
    end if
!
end subroutine
