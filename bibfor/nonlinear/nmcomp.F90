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
! aslint: disable=W1504,C1505,W1306,W0413
!
subroutine nmcomp(BEHInteg, &
                  ndim, option, typmod, &
                  instam, instap, &
                  compor, carcri, multComp, &
                  neps, epsm_inp, deps_inp, &
                  nsig, sigm, &
                  vim, &
                  sigp, vip, &
                  ndsde, dsidep, &
                  codret, &
                  l_epsi_varc_)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcvali.h"
#include "asterfort/redece.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8), intent(in) :: ndim
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(in) :: instam, instap
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: multComp
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm_inp(neps), deps_inp(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    real(kind=8), intent(in) :: vim(*)
    character(len=16), intent(in) :: option
    real(kind=8), intent(inout) :: sigp(nsig), vip(*)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(inout) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                          merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(inout) :: codret
    aster_logical, optional, intent(in) :: l_epsi_varc_
!
! --------------------------------------------------------------------------------------------------
!
! Mechanical non-linear behaviours
!
! Main factory for integration
!
! --------------------------------------------------------------------------------------------------
!
! IO  BEHInteg         : parameters for integration of behaviour
! In  option           : option to compute
! In  typmod           : type of modeling (3D, 2D, etc.)
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
!
! --------------------------------------------------------------------------------------------------
    aster_logical :: conv_cp, l_epsi_varc, lMatr, lVari, lSigm, lMatrPred, lPred, invert
    aster_logical :: lStrainMeca, l_deborst
    integer(kind=8) :: icp, numlc, nvi_all, nvi, k, l, ndimsi
    integer(kind=8) :: codret_vali, codret_ldc, codret_cp
    real(kind=8):: prec
    real(kind=8):: epsm_meca(neps), deps_meca(neps), epsm(neps), deps(neps)
    real(kind=8) :: dsidep_cp(merge(nsig, 6, nsig*neps .eq. ndsde), &
                              merge(neps, 6, nsig*neps .eq. ndsde))
    real(kind=8), allocatable:: vip_cp(:), ka3_min, k3a_min, c_min
    character(len=8)  :: typmod_cp(2), typ_crit
    character(len=16) :: option_cp, defoComp
    type(Behaviour_Integ) :: BEHintegCP
!
! --------------------------------------------------------------------------------------------------
!
    l_epsi_varc = ASTER_TRUE
    if (present(l_epsi_varc_)) then
        l_epsi_varc = l_epsi_varc_
    end if

! - Initialisations
    codret_ldc = LDC_ERROR_NONE
    codret_cp = 0
    codret_vali = 0

! - Variables protegees (in)
    epsm = epsm_inp
    deps = deps_inp

    ! Size controls
    ndimsi = BEHInteg%behavPara%ndimsi
    if (BEHInteg%behavPara%lNdimsiEnabled) then
        ASSERT(ndimsi .le. 6)
        ASSERT(ndimsi .le. neps)
        ASSERT(ndimsi .le. nsig)
    end if
    ASSERT(neps*nsig .eq. ndsde .or. (ndsde .eq. 36 .and. neps .le. 9 .and. nsig .le. 6))

! - Parameters of behaviour of the current integration point
    numlc = BEHInteg%behavPara%numlc
    l_deborst = compor(PLANESTRESS) (1:7) .eq. 'DEBORST'
    ASSERT(.not. l_deborst .or. BEHInteg%behavPara%lDeborstEnabled)
    if (l_deborst) then
        read (compor(NVAR), '(I16)') nvi_all
    else
        nvi_all = BEHInteg%behavPara%nvi
    end if
    lStrainMeca = BEHInteg%behavPara%lStrainMeca
    defoComp = compor(DEFO)

! - Option (operators) to compute
    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)
    lMatrPred = L_MATR_PRED(option)
    lPred = L_PRED(option)

! --------------------------------------------------------------------------------------------------
!   Modification des parametres en entree
! --------------------------------------------------------------------------------------------------

! En contraintes planes, EPZZ est stocke dans les variables internes
! a noter que le mecanisme vip_k contient vip_(k-1) est utilise en cours d'iterations
    if (l_deborst) then
        epsm(3) = vim(nvi_all)
        if (.not. lVari) then
            deps(3) = 0.d0
        else
            deps(3) = vip(nvi_all)-vim(nvi_all)
        end if
    end if

! En phase de prediction / defo_meca, deps est tel que deps_meca = 0 (structure additive defos)
    if (lStrainMeca .and. lPred) then
! ----- Detect external state variables
        call detectVarc(BEHInteg)

! ----- Prepare external state variables at Gauss point
        call behaviourPrepESVAPoin(BEHInteg)

! ----- Prepare input strains for the behaviour law
        epsm_meca = epsm
        deps_meca = 0
        call behaviourPrepStrain(BEHInteg, epsm_meca, deps_meca)
        deps = -deps_meca
    end if

! --------------------------------------------------------------------------------------------------
!   Integration standard du comportement
! --------------------------------------------------------------------------------------------------

    if (.not. l_deborst) then
        call redece(BEHInteg, &
                    ndim, option, typmod, &
                    instam, instap, &
                    compor, carcri, multComp, &
                    neps, epsm, deps, &
                    nsig, sigm, &
                    nvi_all, vim, &
                    sigp, vip, &
                    ndsde, dsidep, codret_ldc, &
                    l_epsi_varc, numlc)

        if (codret_ldc .eq. LDC_ERROR_NCVG) goto 900

! --------------------------------------------------------------------------------------------------
!  Resolution des contraintes planes sizz=0 par une methode de Newton pour les lois non equipees
! --------------------------------------------------------------------------------------------------
    else
        ASSERT(ndim .eq. 2)
        ASSERT(compor(DEFO) .eq. 'PETIT')

!------ Modification des parametres
        BEHintegCP = BEHInteg
        typmod_cp(1) = 'AXIS'
        typmod_cp(2) = typmod(2)
        call behaviourPrepModel(typmod_cp, BEHintegCP)
        nvi = nvi_all-1

! ----- Definition du critere de convergence
        prec = carcri(RESI_DEBORST_MAX)
        ASSERT(prec .ne. r8vide())
        if (prec .ge. 0.d0) then
            typ_crit = 'ABSOLU'
        else
            typ_crit = 'RELATIF'
            prec = -prec
        end if

! ----- S'il faut calculer les contraintes, determination de epzz par methode de Newton
        if (lSigm) then
! --------- Creation de l'espace des variables internes si necessaire
            allocate (vip_cp(nvi))
            if (lVari) then
                vip_cp(1:nvi) = vip(1:nvi)
            else
                vip_cp(1:nvi) = vim(1:nvi)
            end if
            BEHintegCP%behavPara%nvi = nvi

            do icp = 1, nint(carcri(ITER_DEBORST_MAX))

! ------------- Choix de l'option pour accéder à la matrice tangente pour methode de Newton
                if (icp .eq. 1 .and. lMatrPred) then
                    option_cp = 'RIGI_MECA_TANG'
                else
                    option_cp = 'FULL_MECA'
                end if

! ------------- Integration du comportement
                call redece(BEHintegCP, &
                            ndim, option_cp, typmod_cp, &
                            instam, instap, &
                            compor, carcri, multComp, &
                            neps, epsm, deps, &
                            nsig, sigm, &
                            nvi, vim, &
                            sigp, vip_cp, &
                            ndsde, dsidep_cp, codret_ldc, &
                            l_epsi_varc, numlc)

                if (codret_ldc .eq. LDC_ERROR_NCVG) then
                    deallocate (vip_cp)
                    goto 900
                end if

! ------------- Test de convergence
                if (typ_crit .eq. 'ABSOLU') then
                    conv_cp = abs(sigp(3)) .le. prec
                else
                    conv_cp = abs(sigp(3)) .le. prec*maxval(abs(sigp(1:ndimsi)))
                end if
                if (conv_cp) exit

! ------------- Reactualisation de la deformation EPZZ en verifiant l'inversibilite
                if (abs(dsidep_cp(3, 3)) .eq. 0 .and. abs(sigp(3)) .eq. 0) then
                    invert = ASTER_FALSE
                else if (abs(dsidep_cp(3, 3)) .gt. abs(sigp(3))) then
                    invert = ASTER_TRUE
                else
                    invert = abs(dsidep_cp(3, 3))/abs(sigp(3)) .gt. r8prem()
                end if

                if (invert) then
                    deps(3) = deps(3)-sigp(3)/dsidep_cp(3, 3)
                else
                    ! Pivot nul
                    exit
                end if
            end do
            deallocate (vip_cp)
        end if

! ----- Integration du comportement avec le bon epzz et l'option reelle
        BEHInteg%behavPara%nvi = nvi
        call redece(BEHInteg, &
                    ndim, option, typmod_cp, &
                    instam, instap, &
                    compor, carcri, multComp, &
                    neps, epsm, deps, &
                    nsig, sigm, &
                    nvi, vim, &
                    sigp, vip, &
                    ndsde, dsidep, codret_ldc, &
                    l_epsi_varc, numlc)

        if (codret_ldc .eq. LDC_ERROR_NCVG) goto 900

! ----- Test de convergence des contraintes planes pour le code retour (0=OK, 1=NON CVG)
        if (lSigm) then
            if (typ_crit .eq. 'ABSOLU') then
                codret_cp = merge(0, 1, abs(sigp(3)) .le. prec)
            else
                codret_cp = merge(0, 1, abs(sigp(3)) .le. prec*maxval(abs(sigp(1:ndimsi))))
            end if
        end if

! ----- Correction de la matrice tangente pour tenir compte des contraintes planes
        if (lMatr) then
            ! pivot nul -> on ne corrige pas la matrice
            ka3_min = min(minval(abs(dsidep(1:2, 3))), abs(dsidep(4, 3)))
            k3a_min = min(minval(abs(dsidep(3, 1:2))), abs(dsidep(3, 4)))
            c_min = ka3_min*k3a_min
            if (abs(dsidep(3, 3)) .eq. 0 .and. c_min .eq. 0) then
                invert = ASTER_FALSE
            else if (abs(dsidep(3, 3)) .gt. c_min) then
                invert = ASTER_TRUE
            else
                invert = abs(dsidep(3, 3))/c_min .gt. r8prem()
            end if

            if (invert) then
                do k = 1, 4
                    if (k .eq. 3) cycle
                    do l = 1, 4
                        if (l .eq. 3) cycle
                        dsidep(k, l) = dsidep(k, l)-dsidep(k, 3)*dsidep(3, l)/dsidep(3, 3)
                    end do
                end do
                dsidep(:, 3) = 0
                dsidep(3, :) = 0
            end if
        end if

! -----  Actualisation de la deformation epzz dans les variables internes
        if (lVari) then
            vip(nvi_all) = epsm(3)+deps(3)
        end if

    end if

! - Prediction: contribution of the thermal stress to the Taylor expansion if needed
    if (lStrainMeca .and. lPred) then
        ASSERT(lSigm .and. lMatr)
        call behaviourPredictionStress(BEHInteg%behavESVA, dsidep(1:ndimsi, 1:ndimsi), &
                                       sigp(1:ndimsi))
    end if

! - Examen du domaine de validité
    if (BEHInteg%behavPara%lCtrlEpsEnabled) then
        call lcvali(BEHInteg%materPara, &
                    defoComp, epsm(1:ndimsi), deps(1:ndimsi), &
                    instam, instap, codret_vali)
    end if

900 continue

! Traitement du code retour par ordre de gravite

! La loi de comportement a echoue (les resultats n'ont pas de signification)
    if (codret_ldc .eq. LDC_ERROR_NCVG) then
        codret = LDC_ERROR_NCVG

! Les contraintes planes n'ont pas converge : le resultat n'est pas acceptable
    else if (codret_cp .eq. 1) then
        codret = LDC_ERROR_CPLA

! Certaines qualites (criteres) ne sont pas respectees
    else if (codret_ldc .eq. LDC_ERROR_QUAL) then
        codret = LDC_ERROR_QUAL

! Le domaine de validite du comportement n'est pas respecte
    else if (codret_vali .eq. LDC_ERROR_DVAL) then
        codret = LDC_ERROR_DVAL

! Tout est satisfaisant
    else if (codret_ldc .eq. LDC_ERROR_NONE .and. codret_cp .eq. 0 .and. codret_vali .eq. 0) then
        codret = LDC_ERROR_NONE
    else
        ASSERT(ASTER_FALSE)
    end if

end subroutine
