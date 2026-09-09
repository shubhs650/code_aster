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
subroutine nmfihm(BEHInteg, &
                  ndim, nddl, nno1, nno2, npg, &
                  lgpg, ipg, wref, vff1, vff2, &
                  idf2, dffr2, option, geom, &
                  ddlm, ddld, iu, ip, sigm, &
                  sigp, vect, matr, vim, vip, &
                  tm, tp, carcri, compor, typmod, &
                  lVect, lMatr, lSigm, codret)
!
    use Behaviour_type
    use Behaviour_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/ejcine.h"
#include "asterfort/gedisc.h"
#include "asterfort/nmcomp.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8) :: ndim, npg, ipg, idf2, lgpg, nno1, nno2, nddl, iu(3, 16)
    integer(kind=8) :: ip(8)
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg), dffr2(ndim-1, nno2, npg)
    real(kind=8) :: wref(npg), geom(ndim, nno2), ddlm(nddl), ddld(nddl), tm, tp
    real(kind=8) :: sigm(2*ndim-1, npg), sigp(2*ndim-1, npg)
    real(kind=8) :: vect(nddl), matr(nddl*nddl)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg)
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    character(len=8), intent(in) :: typmod(2)
    aster_logical, intent(in) :: lSigm, lVect, lMatr
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!  OPTIONS DE MECANIQUE NON LINEAIRE POUR JOINT ET JOINT_HYME 2D ET 3D
!
!    OPTIONS DE CALCUL
!       * RAPH_MECA      : DDL = DDL- + DDDL  ->   SIGP , FINT
!       * FULL_MECA      : DDL = DDL- + DDDL  ->   SIGP , FINT , KTAN
!       * RIGI_MECA_TANG : DDL = DDL-       ->                  KTAN
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM   DIMENSION DE L'ESPACE
! IN  NDDL   NOMBRE DE DEGRES DE LIBERTE TOTAL
! IN  NNO1   NOMBRE DE NOEUDS DE LA FACE POUR LES DEPLACEMENTS
! IN  NNO2   NOMBRE DE NOEUDS DE LA FACE POUR LES PRESSIONS ET LA GEOM
! IN  NPG    NOMBRE DE POINTS DE GAUSS
! IN  LGPG   NOMBRE DE VARIABLES INTERNES
! IN  WREF   POIDS DE REFERENCE DES POINTS DE GAUSS
! IN  VFF1   VALEUR DES FONCTIONS DE FORME (DE LA FACE) POUR U
! IN  VFF2   VALEUR DES FONCTIONS DE FORME (DE LA FACE) POUR P ET X
! IN  DFFR2  DERIVEE DES FONCTIONS DE FORME DE REFERENCE DE P ET X EN G
! IN  MATE   MATERIAU CODE
! IN  OPTION OPTION DE CALCUL
! IN  GEOM   COORDONNEES NOEUDS DE PRESSION (2D:SEG2,3D:TRIA3 OU QUAD4)
! IN  DDLM   VALEURS DES DEGRES DE LIBERTE A L'INSTANT MOINS
! IN  DDLD   VALEURS DES INCREMEMNT DES DEGRES DE LIBERTE
! IN  IU     DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT
! IN  IP     DECALAGE D'INDICE POUR ACCEDER AUX DDL DE PRESSION
! IN  SIGM   CONTRAINTES -         (RAPH_MECA   ET FULL_MECA_*)
! IN  VIM    VARIABLES INTERNES AU DEBUT DU PAS DE TEMPS
! IN  TM     INSTANT -
! IN  TP     INSTANT +
! IN  CRIT   VALEURS DE L'UTILISATEUR POUR LES CRITERES DE CONVERGENCE
! IN  COMPOR NOM DE LA LOI DE COMPORTEMENT
! IN  TYPMOD TYPE DE LA MODELISATION
! OUT SIGP    : CONTRAINTES +         (RAPH_MECA   ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA   ET FULL_MECA_*)
! OUT MATR    : MATRICE DE RIGIDITE   (RIGI_MECA_* ET FULL_MECA_*)
! OUT VECT    : FORCES INTERIEURES    (RAPH_MECA   ET FULL_MECA_*)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: ksp = 1
    aster_logical, parameter :: axi = ASTER_FALSE
    aster_logical :: ifhyme
    integer(kind=8) :: i, j, kk, m, n, os, p, q, kpg, cod(npg)
    real(kind=8) :: dsidep(6, 6), b(2*ndim-1, ndim+1, 2*nno1+nno2)
    real(kind=8) :: sigmo(6), sigma(6), epsm(6), deps(6), wg
    real(kind=8) :: coopg(ndim+1, npg), rot(ndim*ndim)
    real(kind=8) :: presgm, presgd, temp
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    cod = 0

! IFHYME = TRUE  : CALCUL COUPLE HYDRO-MECA
! IFHYME = FALSE : CALCUL MECA SANS HYDRO ET ELIMINATION DES DDL DE PRES
! (FINT_P=0, KTAN_PP=IDENTITE, KTAN_UP=0)
!
    if (typmod(2) .eq. 'EJ_HYME') then
        ifhyme = ASTER_TRUE
    elseif (typmod(2) .eq. 'ELEMJOIN') then
        ifhyme = ASTER_FALSE
    else
        ASSERT(ASTER_FALSE)
    end if

! - CALCUL DES COORDONNEES DES POINTS DE GAUSS
    call gedisc(ndim, nno2, npg, vff2, geom, coopg)

! - Loop on Gauss points
    do kpg = 1, npg
!
!       CALCUL DE LA MATRICE CINEMATIQUE
        call ejcine(ndim, axi, nno1, nno2, vff1(1, kpg), &
                    vff2(1, kpg), wref(kpg), dffr2(1, 1, kpg), geom, wg, &
                    kpg, ipg, idf2, rot, b)
!
!       CALCUL DES DEFORMATIONS (SAUTS ET GRADIENTS DE PRESSION)
        epsm = 0.d0
        deps = 0.d0
!
        do i = 1, ndim
            do j = 1, ndim
                do n = 1, 2*nno1
                    epsm(i) = epsm(i)+b(i, j, n)*ddlm(iu(j, n))
                    deps(i) = deps(i)+b(i, j, n)*ddld(iu(j, n))
                end do
            end do
        end do
!
        do i = ndim+1, 2*ndim-1
            do n = 1, nno2
                epsm(i) = epsm(i)+b(i, ndim+1, 2*nno1+n)*ddlm(ip(n))
                deps(i) = deps(i)+b(i, ndim+1, 2*nno1+n)*ddld(ip(n))
            end do
        end do
!
!       CALCUL DE LA PRESSION AU POINT DE GAUSS
        presgm = 0.d0
        presgd = 0.d0
        do n = 1, nno2
            presgm = presgm+ddlm(ip(n))*vff2(n, kpg)
            presgd = presgd+ddld(ip(n))*vff2(n, kpg)
        end do
!
!       STOCKAGE DE LA PRESSION DE FLUIDE AU PG
!       POUR LA VI DE POST-TRAITEMENT DANS LA LDC
        epsm(2*ndim) = presgm
        deps(2*ndim) = presgd
!
!       COOROT : COORDONNEES DU PG + MATRICE DE ROTATION
!       (MATRICE UTILE POUR LES VI DE POST-TRAITEMENT DANS LA LDC)
        do j = 1, ndim
            BEHInteg%behavESVA%behavESVAGeom%coorElga(kpg, j) = coopg(j, kpg)
        end do
        do j = 1, ndim*ndim
            BEHInteg%behavESVA%behavESVAOther%rotpg(j) = rot(j)
        end do
!
!       CONTRAINTES -
        sigmo = 0.d0
        do n = 1, 2*ndim-1
            sigmo(n) = sigm(n, kpg)
        end do

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- Integrator
        sigma = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    tm, tp, &
                    compor, carcri, multComp, &
                    6, epsm, deps, &
                    6, sigmo, &
                    vim(1, kpg), &
                    sigma, vip(1, kpg), &
                    36, dsidep, cod(kpg))
        if (cod(kpg) .eq. 1) then
            goto 999
        end if

! ----- Stresses
        if (lSigm) then
            do n = 1, 2*ndim-1
                sigp(n, kpg) = sigma(n)
            end do
        end if

! ----- Vector
        if (lVect) then
! --------- Vector (DOF: U)
            do n = 1, 2*nno1
                do i = 1, ndim
                    kk = iu(i, n)
                    temp = 0.d0
                    do j = 1, ndim
                        temp = temp+b(j, i, n)*sigp(j, kpg)
                    end do
                    vect(kk) = vect(kk)+wg*temp
                end do
            end do
! --------- Vector (DOF: P)
            do n = 1, nno2
                kk = ip(n)
                temp = 0.d0
                do i = ndim+1, 2*ndim-1
                    temp = temp+b(i, ndim+1, 2*nno1+n)*sigp(i, kpg)
                end do
                if (ifhyme) then
                    vect(kk) = vect(kk)+wg*temp
                else
                    vect(kk) = 0.d0
                end if
            end do
        end if
! ----- Matrix
        if (lMatr) then
! --------- Matrix [U(I,N),U(J,M)]
            do n = 1, 2*nno1
                do i = 1, ndim
                    os = (iu(i, n)-1)*nddl
                    do m = 1, 2*nno1
                        do j = 1, ndim
                            kk = os+iu(j, m)
                            temp = 0.d0
                            do p = 1, ndim
                                do q = 1, ndim
                                    temp = temp+b(p, i, n)*dsidep(p, q)*b(q, j, m)
                                end do
                            end do
                            matr(kk) = matr(kk)+wg*temp
                        end do
                    end do
                end do
            end do
! --------- Matrix [K:P(N),P(M)]
            do n = 1, nno2
                os = (ip(n)-1)*nddl
                do m = 1, nno2
                    kk = os+ip(m)
                    temp = 0.d0
                    do p = ndim+1, 2*ndim-1
                        do q = ndim+1, 2*ndim-1
                            temp = temp+b(p, ndim+1, 2*nno1+n)*dsidep(p, q)*b(q, ndim+1, 2*nno1+m)
                        end do
                    end do
                    if (ifhyme) then
                        matr(kk) = matr(kk)+wg*temp
                    else
                        if (n .eq. m) then
                            matr(kk) = 1.d0
                        else
                            matr(kk) = 0.d0
                        end if
                    end if
                end do
            end do
! --------- Matrix [P(N),U(J,M)]
            do n = 1, nno2
                os = (ip(n)-1)*nddl
                do m = 1, 2*nno1
                    do j = 1, ndim
                        kk = os+iu(j, m)
                        temp = 0.d0
                        do p = ndim+1, 2*ndim-1
                            do q = 1, ndim
!                               A ANNULE AFIN DE PASSER VERS LA MINIMISATION ALTERNEE
                                temp = temp+b(p, ndim+1, 2*nno1+n)*dsidep(p, q)*b(q, j, m)
                            end do
                        end do
                        if (ifhyme) then
                            matr(kk) = matr(kk)+wg*temp
                        else
                            matr(kk) = 0.d0
                        end if
                    end do
                end do
            end do
! --------- Matrix [U(I,N),P(M)]
            do n = 1, 2*nno1
                do i = 1, ndim
                    os = (iu(i, n)-1)*nddl
                    do m = 1, nno2
                        kk = os+ip(m)
                        temp = -b(1, i, n)*vff2(m, kpg)
                        if (ifhyme) then
                            matr(kk) = matr(kk)+wg*temp
                        else
                            matr(kk) = 0.d0
                        end if
                    end do
                end do
            end do
        end if
    end do

! - Final flag
999 continue
    if (lSigm) then
        call codere(cod, npg, codret)
    end if
!
end subroutine
