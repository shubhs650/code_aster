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
! aslint: disable=W1306, W1504
!
subroutine nmsfin(BEHInteg, &
                  option, typmod, &
                  compor, carcri, &
                  ndim, nno, npg, nddl, &
                  iw, vff, idff, &
                  geomi, ddlm, ddld, &
                  instm, instp, &
                  lgpg, siefm, vim, &
                  siefp, vip, &
                  fint, matr, &
                  lMatr, lVect, lSigm, &
                  codret)
!
    use Behaviour_type
    use Behaviour_module
    use bloc_fe_module, only: prod_bd, prod_sb, prod_bkb, add_fint, add_matr
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmbeps.h"
#include "asterfort/nmcomp.h"
#include "asterfort/rcvalb.h"
#include "asterfort/teattr.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: option
    integer(kind=8), intent(in) :: ndim, nno, npg, nddl, lgpg
    integer(kind=8), intent(in) :: iw, idff
    real(kind=8), intent(in) :: geomi(ndim, nno), instm, instp
    real(kind=8), intent(in) :: vff(nno, npg)
    real(kind=8), intent(in) :: ddlm(nddl), ddld(nddl), siefm(4*ndim, npg)
    real(kind=8), intent(in) :: vim(lgpg, npg)
    real(kind=8), intent(out) :: fint(nddl), matr(nddl, nddl)
    real(kind=8), intent(out) :: siefp(4*ndim, npg), vip(lgpg, npg)
    aster_logical, intent(in) :: lMatr, lVect, lSigm
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     RAPH_MECA, RIGI_MECA_* ET FULL_MECA_* , ELEMENTS MIX_STA
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION  : OPTION DE CALCUL
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO     : NOMBRE DE NOEUDS STANDARDS D'UN ELEMENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  NDDL    : DEGRES DE LIBERTE D'UN ELEMENT ENRICHI
! IN  IW      : PTR. POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME DE DEPLACEMENT
! IN  IDFF    : PTR. DERIVEE DES FONCTIONS DE FORME DE DEPLACEMENT ELEMENT DE REF.
! IN  GEOMI   : COORDONNEES DES NOEUDS (CONFIGURATION INITIALE)
! IN  COMPOR  : COMPORTEMENT
! IN  LGPG    : DIMENSION DU VECTEUR DES VAR. INTERNES POUR 1 PT GAUSS
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTM   : VALEUR DE L'INSTANT T-
! IN  INSTP   : VALEUR DE L'INSTANT T+
! IN  DDLM    : DDL AU PAS T-
! IN  DDLD    : INCREMENT DE DDL ENTRE T- ET T+
! IN  SIGMG    : CONTRAINTES GENERALISEES EN T-
!                SIGMG(1:2*NDIM) CAUCHY
!                SIGMG(2*NDIM,NPES) : SIG_A, SIG_LAM
! IN  VIM     : VARIABLES INTERNES EN T-
! OUT SIGPG    : CONTRAINTES GENERALIEES (RAPH_MECA ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA_*)
! OUT FINT    : FORCES INTERIEURES (RAPH_MECA ET FULL_MECA_*)
! OUT MATR   : MATR. DE RIGIDITE NON SYM. (RIGI_MECA_* ET FULL_MECA_*)
! OUT CODRET  : CODE RETOUR DE L'INTEGRATION DE LA LDC
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: ksp = 1
    real(kind=8), parameter:: rac2 = sqrt(2.d0)
    real(kind=8), dimension(6), parameter  :: vrac2 = (/1.d0, 1.d0, 1.d0, &
                                                        sqrt(2.d0), sqrt(2.d0), sqrt(2.d0)/)
    real(kind=8), dimension(6, 6), parameter:: id = reshape( &
                                               (/1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                                                 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                                                 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, &
                                                 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
                                                 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
                                                 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/6, 6/))
    aster_logical :: axi
    character(len=16) :: formulation
    integer(kind=8) :: kpg, n, i, j
    integer(kind=8) :: xu(ndim, nno), xe(2*ndim, nno)
    integer(kind=8) :: cod(npg)
    integer(kind=8) :: ndu, nde, neu, nee
    real(kind=8)  :: r, dff(nno, ndim), poids
    real(kind=8)  :: dum(ndim, nno), dup(ndim, nno)
    real(kind=8)  :: dem(2*ndim, nno), dep(2*ndim, nno)
    real(kind=8)  :: bu(2*ndim, ndim, nno), be(2*ndim, 2*ndim, nno)
    real(kind=8)  :: b_dev(2*ndim, ndim, nno), b_vol(2*ndim, ndim, nno)
    real(kind=8)  :: b(2*ndim, ndim, nno), bbar_vol(2*ndim, ndim, nno)
    real(kind=8)  :: epum(2*ndim), epup(2*ndim)
    real(kind=8)  :: epem(2*ndim), epep(2*ndim)
    real(kind=8)  :: siefup(2*ndim), siefep(2*ndim)
    real(kind=8)  :: epsim(6), epsip(6), tau(1), coeff_trace
    real(kind=8)  :: siefcm(6), siefcp(6), dsdepsi(6, 6)
    real(kind=8)  :: kefee(2*ndim, 2*ndim), kefeu(2*ndim, 2*ndim)
    real(kind=8)  :: kefue(2*ndim, 2*ndim), kefuu(2*ndim, 2*ndim)

! --------------------------------------------------------------------------------------------------
!
    axi = typmod(1) .eq. 'AXIS'
    call teattr('S', 'FORMULATION', formulation)
    coeff_trace = 1.d0/ndim

    ! Nombre de ddls
    ndu = ndim
    nde = 2*ndim
    ! Nombre de deformations generalisees
    neu = 2*ndim
    nee = 2*ndim

    if (lVect) fint = 0
    if (lMatr) matr = 0
    dsdepsi = 0
    cod = 0

    ! tableaux de reference bloc (depl,epsi) -> numero du ddl
    forall (i=1:ndu, n=1:nno) xu(i, n) = (n-1)*(ndu+nde)+i
    forall (i=1:nde, n=1:nno) xe(i, n) = (n-1)*(ndu+nde)+ndu+i

    ! Decompactage des ddls en t- et t+
    forall (i=1:ndu, n=1:nno) dum(i, n) = ddlm(xu(i, n))
    forall (i=1:nde, n=1:nno) dem(i, n) = ddlm(xe(i, n))
    forall (i=1:ndu, n=1:nno) dup(i, n) = dum(i, n)+ddld(xu(i, n))
    forall (i=1:nde, n=1:nno) dep(i, n) = dem(i, n)+ddld(xe(i, n))

    ! Formulation STA_INCO : calcul de Bbar_vol
    if (formulation .eq. "STA_INCO") then
        bbar_vol = 0.d0
        do kpg = 1, npg
            b_vol = 0.0d0
            b = 0.d0
            ! Calcul des derivees des fonctions de forme P1, du rayon r et des poids
            call dfdmip(ndim, nno, axi, geomi, kpg, iw, vff(1, kpg), idff, r, poids, dff)
            ! Calcul de la partie volumique de B
            forall (i=1:ndim, j=1:ndim) b_vol(i, j, :) = coeff_trace*dff(:, j)
            ! Calcul de Bbar_vol
            bbar_vol = bbar_vol+1.d0/npg*b_vol
        end do
    end if

    gauss: do kpg = 1, npg
! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- Lecture du paramètre TAU_EPSI
        tau = 0.d0
        call rcvalb(BEHInteg%materPara%schemePara%fami, &
                    BEHInteg%materPara%schemePara%kpg, &
                    BEHInteg%materPara%schemePara%ksp, &
                    '+', &
                    BEHInteg%materPara%jvMaterCode, &
                    ' ', 'NON_LOCAL', &
                    0, ' ', [0.d0], &
                    1, 'TAU_EPSI', tau(1), cod(kpg), 1)
        ! -----------------------!
        !  ELEMENTS CINEMATIQUES !
        ! -----------------------!
        ! Calcul des derivees des fonctions de forme P1, du rayon r et des poids
        call dfdmip(ndim, nno, axi, geomi, kpg, iw, vff(1, kpg), idff, r, poids, dff)
        bu = 0.d0
        ! Formulation STA : calcul de la matrice Bu = B
        if (formulation .eq. "STA") then
            call nmbeps(axi, r, vff(:, kpg), dff, bu)

            ! Formulation STA_INCO : calcul de la matrice Bu = B_dev + Bbar_vol
        else if (formulation .eq. "STA_INCO") then
            b = 0.d0
            b_vol = 0.d0
            b_dev = 0.d0
            ! Calcul de la matrice B
            call nmbeps(axi, r, vff(:, kpg), dff, b)
            ! Calcul de la partie volumique de B
            forall (i=1:ndim, j=1:ndim) b_vol(i, j, :) = coeff_trace*dff(:, j)
            ! Calcul de la partie déviatorique de B
            b_dev = b-b_vol
            ! Calcul de Bu = B_dev + Bbar_vol
            bu = b_dev+bbar_vol
        else
            ASSERT(ASTER_FALSE)
        end if

        ! Calcul des déformations B.U en t- et t+
        epum = prod_bd(bu, dum)
        epup = prod_bd(bu, dup)

        ! Calcul de la matrice N et de la déformation N.E en t- et t+
        be = 0.d0
        if (ndim .eq. 2) then
            be(1, 1, :) = vff(:, kpg)
            be(2, 2, :) = vff(:, kpg)
            be(3, 3, :) = vff(:, kpg)
            be(4, 4, :) = rac2*vff(:, kpg)
        else if (ndim .eq. 3) then
            be(1, 1, :) = vff(:, kpg)
            be(2, 2, :) = vff(:, kpg)
            be(3, 3, :) = vff(:, kpg)
            be(4, 4, :) = rac2*vff(:, kpg)
            be(5, 5, :) = rac2*vff(:, kpg)
            be(6, 6, :) = rac2*vff(:, kpg)
        end if
        epem = prod_bd(be, dem)
        epep = prod_bd(be, dep)

        ! Calcul de la deformation epsi = tau*BU + (1-tau) * NE en t- et t+
        epsim = 0.d0
        epsip = 0.d0
        epsim(1:2*ndim) = tau(1)*epum(1:2*ndim)+(1-tau(1))*epem(1:2*ndim)
        epsip(1:2*ndim) = tau(1)*epup(1:2*ndim)+(1-tau(1))*epep(1:2*ndim)

        ! -----------------------!
        !   LOI DE COMPORTEMENT  !
        ! -----------------------!

! ----- Integrator
        siefcm = 0.d0
        siefcm(1:2*ndim) = siefm(1:2*ndim, kpg)*vrac2(1:2*ndim)
        siefcp = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    instm, instp, &
                    compor, carcri, multComp, &
                    6, epsim, epsip-epsim, &
                    6, siefcm, &
                    vim(1, kpg), &
                    siefcp, vip(1, kpg), &
                    36, dsdepsi, &
                    cod(kpg))
        if (cod(kpg) .eq. 1) goto 999

        ! ----------------------------------------!
        !   FORCES INTERIEURES ET CONTRAINTES EF  !
        ! ----------------------------------------!
        if (lSigm) then
            siefup = siefcp(1:2*ndim)
            siefep = epup-epsip(1:2*ndim)
        end if

        if (lVect) then
            call add_fint(fint, xu, poids*prod_sb(siefup, bu))
            call add_fint(fint, xe, poids*prod_sb(siefep, be))
        end if

        if (lSigm) then
            siefp(1:2*ndim, kpg) = siefup(1:2*ndim)/vrac2(1:2*ndim)
            siefp(2*ndim+1:4*ndim, kpg) = siefep
        end if

        ! -----------------------!
        !    MATRICE TANGENTE    !
        ! -----------------------!

        if (lMatr) then

            ! Blocs de la matrice tangente
            kefuu = tau(1)*dsdepsi(1:2*ndim, 1:2*ndim)
            kefue = (1-tau(1))*dsdepsi(1:2*ndim, 1:2*ndim)
            kefeu = (1-tau(1))*id(1:2*ndim, 1:2*ndim)
            kefee = (tau(1)-1)*id(1:2*ndim, 1:2*ndim)

            ! Assemblage des blocs de la matrice EF
            call add_matr(matr, xu, xu, poids*prod_bkb(bu, kefuu, bu))
            call add_matr(matr, xu, xe, poids*prod_bkb(bu, kefue, be))
            call add_matr(matr, xe, xu, poids*prod_bkb(be, kefeu, bu))
            call add_matr(matr, xe, xe, poids*prod_bkb(be, kefee, be))

        end if

    end do gauss

! - SYNTHESE DES CODES RETOURS
999 continue
    if (lSigm) then
        call codere(cod, npg, codret)
    end if

end subroutine
