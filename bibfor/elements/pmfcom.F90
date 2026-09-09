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
subroutine pmfcom(materPara, &
                  option, carcri, &
                  kpg, debsp, pmfCompor, &
                  nf, instam, instap, nbvalc, &
                  defam, defap, varim, varimp, contm, &
                  defm, ddefp, epsm, modf, sigf, &
                  varip, codret)
!
    use pmfcom_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nm1dco.h"
#include "asterfort/nm1vil.h"
#include "asterfort/paeldt.h"
#include "asterfort/pmf_mazars_unilater.h"
#include "asterfort/pmf_vmis.h"
#include "asterfort/rcexistvarc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "MultiFiber_type.h"
!
    type(Material_Para), intent(inout) :: materPara
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    integer(kind=8) :: nf, nbvalc, kpg, debsp, codret
    real(kind=8) :: contm(nf), defm(nf), ddefp(nf), modf(nf), sigf(nf)
    real(kind=8) :: varimp(nbvalc*nf), varip(nbvalc*nf), varim(nbvalc*nf)
    real(kind=8) :: instam, instap, epsm
    real(kind=8) :: defap(*), defam(*)
    character(len=24) :: pmfCompor(*)
!
! --------------------------------------------------------------------------------------------------
!
! COMPORTEMENT DES ÉLÉMENTS DE POUTRE MULTI-FIBRES
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       kpg     : numéro de point de gauss
!       debsp   : numéro de sous-point de la première fibre du groupe
!       option  : option de calcul
!       pmfCompor  : information sur le comportement du groupe de fibres
!       carcri    : critères de convergence locaux
!       nf      : nombre de fibres du groupe
!       instam  : instant du calcul précédent
!       instap  : instant du calcul
!       icdmat  : code matériau
!       nbvalc  : nombre de variable internes
!       defam   : déformations anélastiques a l'instant précédent
!       defap   : déformations anélastiques a l'instant du calcul
!       varim   : variables internes moins
!       varimp  : variables internes itération précédente (pour DE BORST)
!       contm   : contraintes moins par fibre
!       defm    : déformation  a l'instant du calcul precedent
!       ddefp   : incrément de déformation
!       epsm    : déformation a l'instant préck, édent sur l'élément de structure
!
!   OUT
!       modf    : module tangent des fibres
!       sigf    : contrainte a l'instant actuel des fibres
!       varip   : variables internes a l'instant actuel
!       codret :
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: kpgTemp = 1, kspTemp = 1
    character(len=8), parameter :: famiTemp = 'RIGI'
    type(Material_Para) :: materParaTemp
    integer(kind=8), parameter :: nbProp = 1
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'E'/)
    integer(kind=8) ::  codrep, ksp, fib, ivari, nbvari_grfibre
    real(kind=8) :: ep, em, depsth, tempm, tempp
    real(kind=8) :: depsm, nu
    character(len=8) :: materPoin
    character(len=16) :: relaComp, algoInte
    character(len=30) :: valkm(3)
    aster_logical :: istemp
    type(pmfcom_user) :: for_pmf
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    codrep = 0

! - Behaviour on current cell
    materPoin = pmfCompor(MULTI_FIBER_MATER) (1:8)
    relaComp = pmfCompor(MULTI_FIBER_RELA) (1:16)
    algoInte = pmfCompor(MULTI_FIBER_ALGO) (1:16)
    read (pmfCompor(MULTI_FIBER_NBVARI), '(I24)') nbvari_grfibre

! - Set local coordinate system
    if ((relaComp .eq. 'GRAN_IRRA_LOG') .or. &
        (relaComp .eq. 'VISC_IRRA_LOG')) then
        call initLCSZero(materPara)
    else
        call initLCSNone(materPara)
    end if

! - Copy material parameters with other scheme parameters
    call copyMaterPara(materPara, &
                       famiTemp, kpgTemp, kspTemp, &
                       materParaTemp)

!   Vérification du nombre de fibre
    ASSERT(nbvari_grfibre .le. nbvalc)
!   Attention :
!       nbvari_grfibre : nombre de Vint du comportement
!       nbvalc         : nombre de Vint sur les pts de Gauss de la PMF
!                        nbvalc = max( Vint des comportements sur la PMF )
!       ==> Le décalage c'est donc avec "nbvalc" et pas le nb de Vint du comportement
!
!   Initialisation
    sigf(1:nf) = 0.d0

    if (relaComp .eq. 'ELAS') then
        istemp = rcexistvarc('TEMP')
! ----- Initializations of material parameters on current integration point

        if (.not. istemp) then
            call rcvalb(materParaTemp%schemePara%fami, &
                        materParaTemp%schemePara%kpg, &
                        materParaTemp%schemePara%ksp, &
                        '+', &
                        materParaTemp%jvMaterCode, &
                        materPoin, 'ELAS', &
                        0, '', [0.d0], &
                        nbProp, propName, propVale, &
                        propCode, 1)
            ep = propVale(1)
            em = ep
            do fib = 1, nf
                ksp = debsp-1+fib
                modf(fib) = ep
                sigf(fib) = ep*(contm(fib)/em+ddefp(fib))
            end do
        else
            do fib = 1, nf
                ksp = debsp-1+fib
                call initParaPoin(kpg, ksp, materPara)
! ------------- Get thermal strain
                call paeldt(materPara%schemePara%kpg, &
                            materPara%schemePara%ksp, &
                            materPara%schemePara%fami, &
                            'T', materPara%jvMaterCode, &
                            materPoin, &
                            em, ep, nu, depsth)
                modf(fib) = ep
                sigf(fib) = ep*(contm(fib)/em+ddefp(fib)-depsth)
            end do
        end if
    else if (relaComp .eq. 'MAZARS_UNIL') then
        for_pmf%kpg = kpg
        for_pmf%icdmat = materPara%jvMaterCode
        for_pmf%option = option
        for_pmf%debsp = debsp
        for_pmf%instam = instam
        for_pmf%instap = instap
        for_pmf%epsm = epsm
        call pmf_mazars_unilater(for_pmf, nf, nbvalc, &
                                 pmfCompor, carcri, defam, defap, varim, &
                                 varimp, contm, defm, ddefp, modf, &
                                 sigf, varip, codret)

    else if ((relaComp .eq. 'VMIS_CINE_GC') .or. &
             (relaComp .eq. 'VMIS_CINE_LINE') .or. &
             (relaComp .eq. 'VMIS_ISOT_LINE') .or. &
             (relaComp .eq. 'VMIS_ISOT_TRAC')) then
        for_pmf%kpg = kpg
        for_pmf%icdmat = materPara%jvMaterCode
        for_pmf%option = option
        for_pmf%debsp = debsp
        for_pmf%instam = instam
        for_pmf%instap = instap
        for_pmf%epsm = epsm
        call pmf_vmis(for_pmf, nf, nbvalc, &
                      pmfCompor, materPara, &
                      varim, contm, &
                      ddefp, modf, &
                      sigf, varip, codret)

    else if (relaComp .eq. 'CORR_ACIER') then
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(materParaTemp%schemePara%fami, &
                        materParaTemp%schemePara%kpg, &
                        materParaTemp%schemePara%ksp, &
                        '+', &
                        materParaTemp%jvMaterCode, &
                        materPoin, 'ELAS', &
                        0, '', [0.d0], &
                        nbProp, propName, propVale, &
                        propCode, 1)
            ep = propVale(1)
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                depsm = ddefp(fib)
! ------------- Initializations of material parameters on current integration point
                call initParaPoin(kpg, fib, materPara)

! ------------- Integration
                call nm1dco(materPara, option, carcri, &
                            materPoin, &
                            ep, contm(fib), defm(fib), depsm, &
                            varim(ivari), sigf(fib), varip(ivari), modf(fib), &
                            codret)
                if (codret .ne. 0) goto 999
            end do
        else
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
! ------------- Initializations of material parameters on current integration point
                call initParaPoin(kpg, ksp, materPara)

! ------------- Get thermal strain
                call paeldt(materPara%schemePara%kpg, &
                            materPara%schemePara%ksp, &
                            materPara%schemePara%fami, &
                            '+', materPara%jvMaterCode, &
                            materPoin, &
                            em, ep, nu, depsth)
                depsm = ddefp(fib)-depsth

! ------------- Integration
                call nm1dco(materPara, option, carcri, &
                            materPoin, &
                            ep, contm(fib), defm(fib), depsm, &
                            varim(ivari), sigf(fib), varip(ivari), modf(fib), &
                            codret)
                if (codret .ne. 0) goto 999
            end do
        end if

    else if ((relaComp .eq. 'GRAN_IRRA_LOG') .or. &
             (relaComp .eq. 'VISC_IRRA_LOG')) then

        if (algoInte(1:10) .ne. 'ANALYTIQUE') then
            valkm(1) = relaComp
            valkm(2) = 'DEFI_COMPOR/MULTIFIBRE'
            valkm(3) = algoInte(1:10)
            call utmess('F', 'COMPOR5_81', nk=3, valk=valkm)
        end if
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call utmess('F', 'COMPOR5_40', sk=relaComp)
        end if

        do fib = 1, nf
            ivari = nbvalc*(fib-1)+1
            ksp = debsp-1+fib

! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, materPara)

! --------- Get thermal strain
            call paeldt(materPara%schemePara%kpg, &
                        materPara%schemePara%ksp, &
                        materPara%schemePara%fami, &
                        'T', materPara%jvMaterCode, &
                        materPoin, &
                        em, ep, nu, depsth, &
                        tmoins=tempm, tplus=tempp)
            depsm = ddefp(fib)-depsth

! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, fib, materPara)

! --------- Integration
            call nm1vil(materPara, &
                        relaComp, carcri, &
                        materPoin, &
                        instam, instap, tempm, tempp, &
                        depsm, contm(fib), varim(ivari), &
                        defam(1), defap(1), &
                        sigf(fib), varip(ivari), &
                        modf(fib), codret, nbvalc)

            if (codret .ne. 0) goto 999
        end do

    else
        call utmess('F', 'ELEMENTS2_39', sk=relaComp)
    end if
!
999 continue
end subroutine
