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
! aslint: disable=W1306,W1504
!
subroutine dilele(ds_dil, option, typmod, &
                  BEHInteg, compor, carcri, &
                  ndim, nnos, nnom, npg, nddl, dimdef, &
                  iw, vff, &
                  vffb, idff, idffb, geomi, &
                  lgpg, instam, instap, &
                  ddlm, ddld, siefm, vim, &
                  siefp, vip, fint, matr, &
                  lMatr, lVect, lSigm, codret)
!
    use dil_type
    use bloc_fe_module, only: prod_bd, prod_sb, prod_bkb, add_fint, add_matr
    use Behaviour_type
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/dfdmip.h"
#include "asterfort/dil2gr.h"
#include "asterfort/dilpen.h"
#include "asterfort/nmbeps.h"
#include "asterfort/nmcomp.h"
!
    type(dil_modelisation), intent(in) :: ds_dil
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: typmod(2)
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    integer(kind=8), intent(in) :: ndim, nnos, nnom, npg, nddl, lgpg, dimdef
    integer(kind=8), intent(in) :: iw, idff, idffb
    real(kind=8), intent(in) :: instam, instap
    real(kind=8), intent(in) :: geomi(ndim, nnos+nnom)
    real(kind=8), intent(in) :: vff(nnos+nnom, npg), vffb(nnos, npg)
    real(kind=8), intent(in) :: ddlm(nddl), ddld(nddl)
    real(kind=8), intent(in) :: siefm(dimdef*npg), vim(lgpg*npg)
    real(kind=8), intent(inout) :: siefp(dimdef*npg), vip(lgpg*npg)
    real(kind=8), intent(inout) :: fint(nddl), matr(nddl, nddl)
    aster_logical, intent(in) :: lVect, lMatr, lSigm
    integer(kind=8), intent(inout) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  CALCUL  DES OPTIONS RIGI_MECA_*, RAPH_MECA ET FULL_MECA_*
!           EN PETITES DEFORMATIONS D_PLAN_DIL_U
!           POUR SECOND GRADIENT DE DILATATION : XXXX_DIL
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION  : OPTION DE CALCUL
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNOS    : NOMBRE DE NOEUDS SOMMETS
! IN  NNOM    : NOMBRE DE NOEUDS MILIEUX
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  DIMDEF  : DIMENSION DEF ET CONT GENERALISEES
! IN  NDDL    : NOMBRE DE DEGRES DE LIBERTE D'UN ELEMENT
! IN  DIMDEF  : NOMBRE DE DEFORMATIONS GENERALISEES PAR PG
! IN  IW      : PTR. POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME QUADRATIQUES
! IN  VFFB    : VALEUR  DES FONCTIONS DE FORME LINEAIRES
! IN  IDFF    : PTR. DERIVEE DES FONCTIONS DE FORME QUADRATIQUES ELEMENT DE REF.
! IN  IDFFB   : PTR. DERIVEE DES FONCTIONS DE FORME LINEAIRES ELEMENT DE REF.
! IN  GEOMI   : COORDONNEES DES NOEUDS (CONFIGURATION INITIALE)
! IN  COMPOR  : COMPORTEMENT
! IN  MATE    :
! IN  LGPG    : DIMENSION DU VECTEUR DES VAR. INTERNES POUR 1 PT GAUSS
! IN  CRIT    :
! IN  INSTAM  :
! IN  INSTAP  :
! IN  DDLM    : DDL AU PAS T-
! IN  DDLD    : INCREMENT DE DDL ENTRE T- ET T+
! IN  SIEFM   : CONTRAINTES GENERALISEES EN T-
! IN  VIM     : VARIABLES INTERNES EN T-
! OUT SIEFP   : CONTRAINTES GENERALIEES (RAPH_MECA ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA_*)
! OUT FINT    : FORCES INTERIEURES (RAPH_MECA ET FULL_MECA_*)
! OUT MATR   : MATR. DE RIGIDITE NON SYM. (RIGI_MECA_* ET FULL_MECA_*)
! OUT CODRET  : CODE RETOUR DE L'INTEGRATION DE LA LDC
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: axi = ASTER_FALSE
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8), parameter :: vrac2(6) = [1.d0, 1.d0, 1.d0, rac2, rac2, rac2]
    real(kind=8), parameter :: tiers = 1.d0/3.d0
    real(kind=8), dimension(6), parameter :: kron = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), dimension(6), parameter :: projhyd = (/tiers, tiers, tiers, 0.d0, 0.d0, 0.d0/)
    real(kind=8), dimension(6, 6), parameter:: projdev = reshape( &
                                               (/2*tiers, -tiers, -tiers, 0.d0, 0.d0, 0.d0, &
                                                 -tiers, 2*tiers, -tiers, 0.d0, 0.d0, 0.d0, &
                                                 -tiers, -tiers, 2*tiers, 0.d0, 0.d0, 0.d0, &
                                                 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
                                                 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
                                                 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/6, 6/))
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: kpg, n, i
    integer(kind=8) :: xu(ndim, nnos+nnom), xg(1, nnos), xp(1, nnos)
    integer(kind=8) :: cod(npg)
    integer(kind=8) :: nnu, nng, nnp, ndu, ndg, ndp, neu, neg, nep
    real(kind=8)  :: rpena
    real(kind=8)  :: dum(ndim, nnos+nnom), dup(ndim, nnos+nnom)
    real(kind=8)  :: dgm(1, nnos), dpm(1, nnos)
    real(kind=8)  :: dgp(1, nnos), dpp(1, nnos)
    real(kind=8)  :: r, dff(nnos+nnom, ndim), dffb(nnos, ndim), poids
    real(kind=8)  :: butmp(2*ndim, ndim, nnos+nnom)
    real(kind=8)  :: bu(2*ndim+1, ndim, nnos+nnom)
    real(kind=8)  :: bg(1+ndim, 1, nnos), bp(1, 1, nnos)
    real(kind=8)  :: epefum(2*ndim+1), epefgm(1+ndim), epefpm(1)
    real(kind=8)  :: epefup(2*ndim+1), epefgp(1+ndim), epefpp(1)
    real(kind=8)  :: siefup(2*ndim+1), siefgp(1+ndim), siefpp(1)
    real(kind=8)  :: epl1gm(6), epl1gp(6)
    real(kind=8)  :: sigm1g(6), sigp1g(6)
    real(kind=8)  :: eplcp(ndim), silcp(ndim), deps(2*ndim)
    real(kind=8)  :: dsde1g(6, 6), dsde2g(ndim, ndim)
    real(kind=8)  :: kefuu(2*ndim+1, 2*ndim+1), kefug(2*ndim+1, 1+ndim), kefup(2*ndim+1, 1)
    real(kind=8)  :: kefgu(1+ndim, 2*ndim+1), kefgg(1+ndim, 1+ndim), kefgp(1+ndim, 1)
    real(kind=8)  :: kefpu(1, 2*ndim+1), kefpg(1, 1+ndim), kefpp(1, 1)
    real(kind=8)  :: dev(2*ndim, 2*ndim), hyd(2*ndim), kr(2*ndim)
!
! --------------------------------------------------------------------------------------------------
!
    !Nombre de noeuds
    nnu = nnos+nnom
    nng = nnos
    nnp = nnos
    !Nombre de ddls
    ndu = ndim
    ndg = 1
    ndp = 1
    !Nombre de defs généralisées
    neu = 1+2*ndim
    neg = 1+ndim
    nep = 1

    kr = kron(1:2*ndim)
    hyd = projhyd(1:2*ndim)
    dev = projdev(1:2*ndim, 1:2*ndim)

    if (lVect) fint = 0
    if (lMatr) matr = 0
    cod = 0

! - Initialisation pour comportement second gradient
    call dilpen(BEHInteg%materPara%jvMaterCode, rpena)

! - tableaux de reference bloc (depl,gonf,pres) -> numero du ddl
    forall (i=1:ndu, n=1:nng) xu(i, n) = (n-1)*(ndu+ndg+ndp)+i
    forall (i=1:ndp, n=1:nnp) xp(i, n) = (n-1)*(ndu+ndg+ndp)+ndu+i
    forall (i=1:ndg, n=1:nng) xg(i, n) = (n-1)*(ndu+ndg+ndp)+ndu+ndp+i
    forall (i=1:ndu, n=nng+1:nnu) xu(i, n) = (ndu+ndg+ndp)*nng+(n-1-nng)*ndu+i

! - Decompactage des ddls en t- et t+
    forall (i=1:ndu, n=1:nnu) dum(i, n) = ddlm(xu(i, n))
    forall (i=1:ndu, n=1:nnu) dup(i, n) = ddlm(xu(i, n))+ddld(xu(i, n))
    forall (i=1:ndg, n=1:nng) dgm(i, n) = ddlm(xg(i, n))
    forall (i=1:ndg, n=1:nng) dgp(i, n) = ddlm(xg(i, n))+ddld(xg(i, n))
    forall (i=1:ndp, n=1:nnp) dpm(i, n) = ddlm(xp(i, n))
    forall (i=1:ndp, n=1:nnp) dpp(i, n) = ddlm(xp(i, n))+ddld(xp(i, n))

    gauss: do kpg = 1, npg

        ! -----------------------!
        !  ELEMENTS CINEMATIQUES !
        ! -----------------------!

        ! Calcul des derivees des fonctions de forme P1
        call dfdmip(ndim, nnos, axi, geomi, kpg, iw, vffb(1, kpg), idffb, r, poids, dffb)

        ! Calcul des derivees des fonctions de forme P2
        call dfdmip(ndim, nnu, axi, geomi, kpg, iw, vff(1, kpg), idff, r, poids, dff)

        ! Assemblage de la matrice BU
        call nmbeps(axi, r, vff(:, kpg), dff, butmp)
        bu = 0.d0
        bu(1:2*ndim, :, :) = butmp
        bu(2*ndim+1, :, :) = butmp(1, :, :)+butmp(2, :, :)+butmp(3, :, :)

        ! Assemblage des matrices BG et BP
        bg = 0.d0
        bp = 0.d0
        bg(1, 1, :) = vffb(:, kpg)
        bg(2:neg, 1, :) = transpose(dffb)
        bp(1, 1, :) = vffb(:, kpg)

        ! Calcul des deformations generalisees aux points de Gauss
        ! Instant -
        epefum = prod_bd(bu, dum)
        epefgm = prod_bd(bg, dgm)
        epefpm = prod_bd(bp, dpm)
        ! Instant +
        epefup = prod_bd(bu, dup)
        epefgp = prod_bd(bg, dgp)
        epefpp = prod_bd(bp, dpp)
        ! Increment deps
        deps = epefup(1:2*ndim)-epefum(1:2*ndim)

        ! -------------------------------------------------------!
        !   LOI DE COMPORTEMENT PREMIER GRADIENT                 !
        ! -------------------------------------------------------!

        !Calcul des déformations pour ldc
        epl1gm = 0.d0
        epl1gp = 0.d0
        if (ds_dil%inco) then
            epl1gm(1:2*ndim) = matmul(dev, epefum(1:2*ndim))+epefgm(1)*hyd
            epl1gp(1:2*ndim) = matmul(dev, epefup(1:2*ndim))+epefgp(1)*hyd
        else
            epl1gm(1:2*ndim) = epefum(1:2*ndim)
            epl1gp(1:2*ndim) = epefup(1:2*ndim)
        end if

        !Format sigma pour ldc
        sigm1g = 0.d0
        sigm1g(1:2*ndim) = siefm(1+(kpg-1)*dimdef:2*ndim+(kpg-1)*dimdef)*vrac2(1:2*ndim)
        if (ds_dil%inco) then
            sigm1g(1:2*ndim) = sigm1g(1:2*ndim)+(siefm(dimdef*(kpg-1)+2*ndim+2) &
                                                 -siefm(dimdef*(kpg-1)+2*ndim+1)+epefpm(1) &
                                                 +rpena*(epefum(2*ndim+1)-epefgm(1)))*kr
        end if

        sigp1g = 0.d0
        dsde1g = 0.d0

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- Integrator
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    instam, instap, &
                    compor, carcri, multComp, &
                    6, epl1gm, epl1gp-epl1gm, 6, sigm1g, &
                    vim(1+lgpg*(kpg-1):lgpg*kpg), &
                    sigp1g, vip(1+lgpg*(kpg-1):lgpg*kpg), &
                    36, dsde1g, &
                    cod(kpg))

        ! -------------------------------------------------------!
        !   LOI DE COMPORTEMENT SECOND GRADIENT DE DILATATION    !
        ! -------------------------------------------------------!

        ! Preparation des deformations generalisees de ldc second gradient
        eplcp = epefgp(2:neg)

        ! Comportement second gradient de dilatation
        call dil2gr(BEHInteg%materPara%jvMaterCode, ndim, ndim, eplcp, silcp, dsde2g)

        ! ----------------------------------------!
        !   FORCES INTERIEURES ET CONTRAINTES EF  !
        ! ----------------------------------------!

        if (lSigm) then
            !Format sigma sortie ldc
            if (ds_dil%inco) then
                siefup(1:2*ndim) = matmul(dev, sigp1g(1:2*ndim))
            else
                siefup(1:2*ndim) = sigp1g(1:2*ndim)
            end if

            ! Contraintes generalisees EF par bloc
            siefup(2*ndim+1) = epefpp(1)+rpena*(epefup(2*ndim+1)-epefgp(1))
            if (ds_dil%inco) then
                siefgp(1) = -siefup(2*ndim+1)+dot_product(sigp1g(1:2*ndim), hyd)
            else
                siefgp(1) = -siefup(2*ndim+1)
            end if
            siefgp(2:neg) = silcp
            siefpp(1) = epefup(2*ndim+1)-epefgp(1)

        end if

        if (lVect) then
            ! Forces interieures au point de Gauss kpg
            call add_fint(fint, xu, poids*prod_sb(siefup, bu))
            call add_fint(fint, xg, poids*prod_sb(siefgp, bg))
            call add_fint(fint, xp, poids*prod_sb(siefpp, bp))
        end if

        if (lSigm) then
            ! Stockage des contraintes generalisees
            if (ds_dil%inco) then
                siefp(dimdef*(kpg-1)+1:dimdef*(kpg-1)+2*ndim) = siefup(1:2*ndim)/vrac2(1:2*ndim) &
                                                                +siefup(2*ndim+1)*kr
            else
                siefp(dimdef*(kpg-1)+1:dimdef*(kpg-1)+2*ndim) = siefup(1:2*ndim)/vrac2(1:2*ndim)
            end if
            siefp(dimdef*(kpg-1)+2*ndim+1) = siefup(2*ndim+1)
            siefp(dimdef*(kpg-1)+2*ndim+2:dimdef*(kpg-1)+2*ndim+2+ndim) = siefgp
            siefp(dimdef*kpg) = siefpp(1)
        end if

        ! -----------------------!
        !    MATRICE TANGENTE    !
        ! -----------------------!

        if (lMatr) then

            ! Construction des blocs de la matrice tangente EF
            kefuu = 0.d0
            if (ds_dil%inco) then
                kefuu(1:2*ndim, 1:2*ndim) = matmul(matmul(dev, dsde1g(1:2*ndim, 1:2*ndim)), dev)
            else
                kefuu(1:2*ndim, 1:2*ndim) = dsde1g(1:2*ndim, 1:2*ndim)
            end if
            kefuu(2*ndim+1, 2*ndim+1) = rpena

            kefug = 0.d0
            if (ds_dil%inco) then
                kefug(1:2*ndim, 1) = matmul(matmul(dev, dsde1g(1:2*ndim, 1:2*ndim)), hyd)
            end if
            kefug(2*ndim+1, 1) = -rpena

            kefgu = 0.d0
            if (ds_dil%inco) then
                kefgu(1, 1:2*ndim) = matmul(matmul(hyd, dsde1g(1:2*ndim, 1:2*ndim)), dev)
            end if
            kefgu(1, 2*ndim+1) = -rpena

            kefgg = 0.d0
            if (ds_dil%inco) then
                kefgg(1, 1) = rpena+dot_product(matmul(hyd, dsde1g(1:2*ndim, 1:2*ndim)), hyd)
            else
                kefgg(1, 1) = rpena
            end if
            kefgg(2:neg, 2:neg) = dsde2g

            kefup = 0.d0
            kefup(2*ndim+1, 1) = 1.0d0

            kefgp = 0.d0
            kefgp(1, 1) = -1.0d0

            kefpu = 0.d0
            kefpu(1, 2*ndim+1) = 1.0d0

            kefpg = 0.d0
            kefpg(1, 1) = -1.0d0

            kefpp = 0.d0

            ! Assemblage des blocs de la matrice EF
            call add_matr(matr, xu, xu, poids*prod_bkb(bu, kefuu, bu))
            call add_matr(matr, xu, xg, poids*prod_bkb(bu, kefug, bg))
            call add_matr(matr, xg, xu, poids*prod_bkb(bg, kefgu, bu))
            call add_matr(matr, xg, xg, poids*prod_bkb(bg, kefgg, bg))
            call add_matr(matr, xu, xp, poids*prod_bkb(bu, kefup, bp))
            call add_matr(matr, xg, xp, poids*prod_bkb(bg, kefgp, bp))
            call add_matr(matr, xp, xu, poids*prod_bkb(bp, kefpu, bu))
            call add_matr(matr, xp, xg, poids*prod_bkb(bp, kefpg, bg))
            call add_matr(matr, xp, xp, poids*prod_bkb(bp, kefpp, bp))

        end if

    end do gauss

    if (lSigm) then
        call codere(cod, npg, codret)
    end if

end subroutine
