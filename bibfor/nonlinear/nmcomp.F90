! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1504,C1505,W1306
!
subroutine nmcomp(BEHinteg, &
                  fami, kpg, ksp, ndim, typmod, &
                  imate, compor, carcri, instam, instap, &
                  neps, epsm_inp, deps_inp, nsig, sigm, &
                  vim, option, angmas, sigp, vip, &
                  ndsde, dsidep, codret, mult_comp_, l_epsi_varc_, &
                  materi_)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterc/r8vide.h"
#include "asterc/r8gaem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcvali.h"
#include "asterfort/redece.h"
!
    type(Behaviour_Integ) :: BEHinteg
    integer :: kpg, ksp, ndim, imate, codret, neps, nsig, ndsde
    character(len=*)    :: fami
    character(len=8)    :: typmod(*)
    character(len=16)   :: compor(*), option
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm_inp(neps), deps_inp(neps)
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    real(kind=8) :: carcri(*), sigm(nsig), vim(*), sigp(nsig), vip(*), angmas(*)
    character(len=8), optional, intent(in) :: materi_
    character(len=16), optional, intent(in) :: mult_comp_
    aster_logical, optional, intent(in) :: l_epsi_varc_
! --------------------------------------------------------------------------------------------------
!     INTEGRATION DES LOIS DE COMPORTEMENT NON LINEAIRE
! --------------------------------------------------------------------------------------------------
!
! in  behinteg       : parameters for integration of behaviour
! in  fami,kpg,ksp  : famille et numero du (sous)point de gauss
!     ndim    : dimension de l'espace
!               3 : 3d , 2 : d_plan ,axis ou  c_plan
!     typmod(2): modelisation ex: 1:3d, 2:inco
!     imate   : adresse du materiau code
!     compor  : comportement :  (1) = type de relation comportement
!                               (2) = nb variables internes / pg
!                               (3) = hypothese sur les deformations
!                               (4) etc... (voir grandeur compor)
!     crit    : criteres de convergence locaux (voir grandeur carcri)
!     instam  : instant du calcul precedent
!     instap  : instant du calcul
!     neps    : nombre de cmp de epsm et deps (suivant modelisation)
!     epsm    : deformations a l'instant du calcul precedent
!     deps    : increment de deformation totale :
!                deps(t) = deps(mecanique(t)) + deps(dilatation(t))
!     nsig    : nombre de cmp de sigm et sigp (suivant modelisation)
!     sigm    : contraintes a l'instant du calcul precedent
!     vim     : variables internes a l'instant du calcul precedent
!     option  : option demandee : rigi_meca_tang , full_meca , raph_meca
!     angmas  : les trois angles du mot_clef massif (affe_cara_elem),
!               + un reel qui vaut 0 si nautiquies ou 2 si euler
!               + les 3 angles d'euler
!
! out sigp    : contraintes a l'instant actuel
! var vip     : variables internes
!                in  : estimation (iteration precedente ou lag. augm.)
!                out : en t+
!     ndsde   : dimension de dsidep
!     dsidep  : operateur tangent dsig/deps ou dsig/df
!     codret  : code retour loi de comporment, par gravite decroissante
!                   1 : echec fatal dans l'integration de la loi (resultats non utilisables)
!                   3 : contraintes planes deborst non convergees (interdit la convergence)
!                   2 : criteres de qualite de la loi non respectes (decoupage si convergence)
!                   4 : domaine de validite de la loi non respecte (emission d'une alarme)
!                   0 : tout va bien
!
! precisions :
! -----------
!  les tenseurs et matrices sont ranges dans l'ordre :
!         xx yy zz sqrt(2)*xy sqrt(2)*xz sqrt(2)*yz
!
! -si deformation = simo_miehe
!   epsm(3,3)    gradient de la transformation en t-
!   deps(3,3)    gradient de la transformation de t- a t+
!
!  output si resi (raph_meca, full_meca_*)
!   vip      variables internes en t+
!   sigp(6)  contrainte de kirchhoff en t+ ranges dans l'ordre
!         xx yy zz sqrt(2)*xy sqrt(2)*xz sqrt(2)*yz
!
!  output si rigi (rigi_meca_*, full_meca_*)
!   dsidep(6,3,3) matrice tangente d(tau)/d(fd) * (fd)t
!                 (avec les racines de 2)
!
! -sinon (deformation = petit ou petit_reac ou gdef_...)
!   epsm(6), deps(6)  sont les deformations (linearisees ou green ou ..)
!
! --------------------------------------------------------------------------------------------------
    aster_logical :: conv_cp, l_epsi_varc, lMatr, lVari, lSigm, lMatrPred, lPred, invert
    aster_logical :: l_defo_meca, l_czm, l_large, l_deborst, l_grad_vari
    integer :: icp, numlc, nvi_all, nvi, k, l, ndimsi
    integer :: codret_vali, codret_ldc, codret_cp
    real(kind=8):: prec
    real(kind=8):: epsm_meca(neps), deps_meca(neps), epsm(neps), deps(neps)
    real(kind=8) :: dsidep_cp(merge(nsig, 6, nsig*neps .eq. ndsde), &
                              merge(neps, 6, nsig*neps .eq. ndsde))
    real(kind=8), allocatable:: vip_cp(:), ka3_min, k3a_min, c_min
    character(len=8)  :: materi
    character(len=8)  :: typmod_cp(2), typ_crit
    character(len=16) :: option_cp, mult_comp, defo_ldc, defo_comp
! --------------------------------------------------------------------------------------------------

    ! Controles
    ASSERT(neps*nsig .eq. ndsde .or. (ndsde .eq. 36 .and. neps .le. 9 .and. nsig .le. 6))

!   Les paramètres optionnels
    mult_comp = ' '
    materi = ' '
    l_epsi_varc = ASTER_TRUE
    if (present(mult_comp_)) mult_comp = mult_comp_
    if (present(materi_)) materi = materi_
    if (present(l_epsi_varc_)) l_epsi_varc = l_epsi_varc_

    ! Variables protegees (in)
    epsm = epsm_inp
    deps = deps_inp

    ! Initialisation
    codret_ldc = 0
    codret_cp = 0
    codret_vali = 0
    read (compor(NUME), '(I16)') numlc
    read (compor(NVAR), '(I16)') nvi_all
    read (compor(DEFO_LDC), '(A16)') defo_ldc
    read (compor(DEFO), '(A16)') defo_comp
    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)
    lMatrPred = L_MATR_PRED(option)
    lPred = L_PRED(option)
    l_defo_meca = defo_ldc .eq. 'MECANIQUE'
    l_czm = typmod(2) .eq. 'ELEMJOIN'
    l_grad_vari = typmod(2) .eq. 'GRADVARI'
    l_large = defo_comp .eq. 'SIMO_MIEHE' .or. defo_comp .eq. 'GROT_GDEP'
    l_deborst = compor(PLANESTRESS) (1:7) .eq. 'DEBORST'

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

    ! En phase de prediction / defo_meca, deps est tel que deps_meca = 0 (strucure additive defos)
    if (l_defo_meca .and. lPred) then

        call behaviourPrepESVAGauss(carcri, defo_ldc, imate, fami, kpg, ksp, neps, instap, BEHinteg)

        epsm_meca = epsm
        deps_meca = 0
        call behaviourPrepStrain(l_czm, l_large, l_defo_meca, l_grad_vari, imate, fami, &
                                 kpg, ksp, neps, BEHinteg%esva, epsm_meca, deps_meca)
        deps = -deps_meca
    end if

! --------------------------------------------------------------------------------------------------
!   Integration standard du comportement
! --------------------------------------------------------------------------------------------------

    if (.not. l_deborst) then

        call redece(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    l_epsi_varc, imate, materi, compor, mult_comp, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi_all, vim, option, &
                    angmas, numlc, sigp, vip, &
                    ndsde, dsidep, codret_ldc)

        if (codret_ldc .eq. 1) goto 900

! --------------------------------------------------------------------------------------------------
!  Resolution des contraintes planes sizz=0 par une methode de Newton pour les lois non equipees
! --------------------------------------------------------------------------------------------------
    else

        ! Controles
        ASSERT(ndim .eq. 2)
        ASSERT(nsig .ge. 2*ndim)
        ASSERT(neps .ge. 2*ndim)
        ASSERT(compor(DEFO) .eq. 'PETIT')

        ! Modification des parametres
        nvi = nvi_all-1
        typmod_cp(1) = 'AXIS'
        typmod_cp(2) = typmod(2)

        ! Definition du critere de convergence
        prec = carcri(RESI_DEBORST_MAX)
        ASSERT(prec .ne. r8vide())
        if (prec .ge. 0.d0) then
            typ_crit = 'ABSOLU'
        else
            typ_crit = 'RELATIF'
            prec = -prec
        end if

        ! S'il faut calculer les contraintes, determination de epzz par methode de Newton
        if (lSigm) then

            ! Creation de l'espace des variables internes si necessaire
            allocate (vip_cp(nvi))
            if (lVari) then
                vip_cp(1:nvi) = vip(1:nvi)
            else
                vip_cp(1:nvi) = vim(1:nvi)
            end if

            do icp = 1, nint(carcri(ITER_DEBORST_MAX))

                ! Choix de l'option pour accéder à la matrice tangente pour methode de Newton
                if (icp .eq. 1 .and. lMatrPred) then
                    option_cp = 'RIGI_MECA_TANG'
                else
                    option_cp = 'FULL_MECA'
                end if

                ! Integration du comportement
                call redece(BEHinteg, &
                            fami, kpg, ksp, ndim, typmod_cp, &
                            l_epsi_varc, imate, materi, compor, mult_comp, &
                            carcri, instam, instap, neps, epsm, &
                            deps, nsig, sigm, nvi, vim, option_cp, &
                            angmas, numlc, sigp, vip_cp, &
                            ndsde, dsidep_cp, codret_ldc)

                if (codret_ldc .eq. 1) then
                    deallocate (vip_cp)
                    goto 900
                end if

                ! Test de convergence
                if (typ_crit .eq. 'ABSOLU') then
                    conv_cp = abs(sigp(3)) .le. prec
                else
                    conv_cp = abs(sigp(3)) .le. prec*maxval(abs(sigp(1:2*ndim)))
                end if
                if (conv_cp) exit

                ! Reactualisation de la deformation EPZZ en verifiant l'inversibilite
                if (abs(dsidep_cp(3, 3)) .eq. 0 .and. abs(sigp(3)) .eq. 0) then
                    invert = ASTER_FALSE
                else if (abs(dsidep_cp(3, 3)) .gt. abs(sigp(3))) then
                    invert = ASTER_TRUE
                else
                    invert = abs(dsidep_cp(3, 3))/abs(sigp(3)) .gt. 1.d0/r8gaem()
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

        ! Integration du comportement avec le bon epzz et l'option reelle
        call redece(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod_cp, &
                    l_epsi_varc, imate, materi, compor, mult_comp, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    angmas, numlc, sigp, vip, &
                    ndsde, dsidep, codret_ldc)

        if (codret_ldc .eq. 1) goto 900

        ! Test de convergence des contraintes planes pour le code retour (0=OK, 1=NON CVG)
        if (lSigm) then
            if (typ_crit .eq. 'ABSOLU') then
                codret_cp = merge(0, 1, abs(sigp(3)) .le. prec)
            else
                codret_cp = merge(0, 1, abs(sigp(3)) .le. prec*maxval(abs(sigp(1:2*ndim))))
            end if
        end if

        ! Correction de la matrice tangente pour tenir compte des contraintes planes
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
                invert = abs(dsidep(3, 3))/c_min .gt. 1.d0/r8gaem()
            end if

            if (invert) then
                do k = 1, 4
                    if (k .eq. 3) goto 136
                    do l = 1, 4
                        if (l .eq. 3) goto 137
                        dsidep(k, l) = dsidep(k, l)-dsidep(k, 3)*dsidep(3, l)/dsidep(3, 3)
137                     continue
                    end do
136                 continue
                end do
                dsidep(:, 3) = 0
                dsidep(3, :) = 0
            end if
        end if

        ! Actualisation de la deformation epzz dans les variables internes
        if (lVari) then
            vip(nvi_all) = epsm(3)+deps(3)
        end if

    end if

! - Prediction: contribution of the thermal stress to the Taylor expansion if needed
    if (l_defo_meca .and. lPred) then
        if (.not. l_czm) then
            ndimsi = 2*ndim
            ! A remettre suite à la fiche issue32329
            !ASSERT(.not. l_large)
            ASSERT(typmod(2) .eq. ' ' .or. typmod(2) .eq. 'GRADVARI' .or. typmod(2) .eq. 'HHO')
            ASSERT(nsig .ge. ndimsi)
            ASSERT(size(dsidep, 1) .ge. ndimsi)
            ASSERT(size(dsidep, 2) .ge. ndimsi)
            ASSERT(lSigm .and. lMatr)

            call behaviourPredictionStress(BEHinteg%esva, dsidep, sigp(1:ndimsi))
        end if
    end if

!   Examen du domaine de validité
    call lcvali(fami, kpg, ksp, imate, materi, &
                compor, ndim, epsm, deps, instam, &
                instap, codret_vali)

900 continue

! Traitement du code retour par ordre de gravite

    ! La loi de comportement a echoue (les resultats n'ont pas de signification)
    if (codret_ldc .eq. 1) then
        codret = 1

        ! Les contraintes planes n'ont pas converge : le resultat n'est pas acceptable
    else if (codret_cp .eq. 1) then
        codret = 3

        ! Certaines qualites (criteres) ne sont pas respectees
    else if (codret_ldc .eq. 2) then
        codret = 2

        ! Le domaine de validite du comportement n'est pas respecte
    else if (codret_vali .eq. 4) then
        codret = 4

        ! Tout est satisfaisant
    else if (codret_ldc .eq. 0 .and. codret_cp .eq. 0 .and. codret_vali .eq. 0) then
        codret = 0

        ! Situation non prevue
    else
        ASSERT(ASTER_FALSE)
    end if

end subroutine
