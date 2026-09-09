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
! aslint: disable=W0413
!
subroutine nmtahe(fami, kpg, ksp, ndim, imate, &
                  carcri, instam, instap, epsm, &
                  deps, sigm, nvi, vim, option, sigp, &
                  vip, dsidep, iret)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmtaac.h"
#include "asterfort/nmtacr.h"
#include "asterfort/nmtadp.h"
#include "asterfort/nmtael.h"
#include "asterfort/nmtama.h"
#include "asterfort/nmtari.h"
#include "asterfort/nmtarl.h"
#include "asterfort/nmtasp.h"
#include "asterfort/nmtaxi.h"
!
    integer(kind=8) :: kpg, ksp, ndim, imate, nvi
    character(len=*) :: fami
    character(len=16) :: option
    real(kind=8) :: carcri(CARCRI_SIZE), instam, instap
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: sigm(6), vim(nvi), sigp(6), vip(nvi), dsidep(6, 6)
! ----------------------------------------------------------------------
!     REALISE LA LOI DE TAHERI POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES DEFORMATIONS
!
!
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  EPSM    : DEFORMATIONS A L'INSTANT DU CALCUL PRECEDENT
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE TANGENTE
!
! OUT IRET    : CODE RETOUR DE L'INTEGRATION DE LA LOI DE TAHERI
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ABSENCE DE CONVERGENCE
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX YY ZZ XY XZ YZ
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ndimsi, niter, k, iret
    integer(kind=8) :: ind
!
    real(kind=8) :: rac2
    real(kind=8) :: matm(3), mat(14)
    real(kind=8) :: sigel(6), epm(6)
    real(kind=8) :: dp, sp, xi
    real(kind=8) :: f, g, fdp, fds, gdp, gds, fdx, gdx, dpmax, sig(6)
    real(kind=8) :: tang(6, 6)
    real(kind=8) :: det, dirdp, dirsp, dirxi, ener, min, rho, rhomax, interi
!
    parameter(rhomax=2.d0, interi=0.99999d0)
!
!
!
! - INITIALISATION
    iret = 0
!   DIMENSION DES TENSEURS ET MISE AUX NORMES
    ndimsi = ndim*2
    rac2 = sqrt(2.d0)
    ASSERT(nvi .eq. 9)
    do k = 4, ndimsi
        vim(2+k) = vim(2+k)*rac2
    end do
!
!    LECTURE DES CARACTERISTIQUES
    call nmtama(fami, kpg, ksp, imate, instam, &
                instap, matm, mat)
!
!
!    CALCUL DES CONTRAINTES ELASTIQUES
    call nmtael(fami, kpg, ksp, imate, ndimsi, &
                matm, mat, sigm, epsm, deps, &
                epm, sigel, sigp)
!
!
! - CALCUL DES CONTRAINTES REELLES ET DES VARIABLES INTERNES
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
!
!      PREDICTION ELASTIQUE
        dp = 0.d0
        xi = 1.d0
        call nmtasp(ndimsi, carcri, mat, sigel, vim, &
                    epm, dp, sp, xi, f, &
                    iret)
        if (iret .eq. 1) goto 999
!
!      CHARGE
        if (f .gt. 0.d0) then
!
!        CALCUL DE DP : EQUATION SCALAIRE F=0 AVEC  0 < DP < DPMAX
            call nmtadp(ndimsi, carcri, mat, sigel, vim, &
                        epm, dp, sp, xi, g, &
                        iret)
!
!        PLASTICITE CLASSIQUE (G<0)
            if (g .le. 0.d0) then
                call nmtaac(2, ndimsi, mat, sigel, vim, &
                            epm, dp, sp, xi, sigp, &
                            vip)
!
!        PLASTICITE A DEUX SURFACES (G>0) -> NEWTON
            else
!
!          ITERATIONS DE NEWTON
                sp = vim(2)
                xi = 1.d0
                do niter = 1, int(carcri(1))
!
!            DIRECTION DE DESCENTE
                    call nmtacr(2, ndimsi, mat, sigel, vim, &
                                epm, dp, sp, xi, f, &
                                g, fds, gds, fdp, gdp, &
                                fdx, gdx, dpmax, sig, tang)
                    det = fdp*gds-fds*gdp
                    dirdp = (g*fds-f*gds)/det
                    dirsp = (f*gdp-g*fdp)/det
                    dirxi = 0.d0
!
!            CORRECTION DE LA DIRECTION POUR RESTER DS P>P- ET S>SP>SP-
                    if (dp+rhomax*dirdp .lt. 0.d0) dirdp = (-dp)/rhomax
                    if (sp+rhomax*dirsp .lt. vim(2)) dirsp = (vim(2)-sp)/rhomax
                    if (sp+rhomax*dirsp .gt. mat(11)) dirsp = (mat(11)-sp)/rhomax
!
!            RECHERCHE LINEAIRE
                    ener = (f**2+g**2)/2.d0
                    min = (f*fdp+g*gdp)*dirdp+(f*fds+g*gds)*dirsp
                    rho = rhomax*interi
                    call nmtarl(2, ndimsi, mat, sigel, vim, &
                                epm, dp, sp, xi, dirdp, &
                                dirsp, dirxi, min, rho, ener)
!
!            ACTUALISATION
                    dp = dp+rho*dirdp
                    sp = sp+rho*dirsp
!
                    if (ener/mat(4)**2 .lt. carcri(3)**2) goto 210
                end do
                iret = 1
                goto 999
210             continue
!
                call nmtaac(3, ndimsi, mat, sigel, vim, &
                            epm, dp, sp, xi, sigp, &
                            vip)
            end if
!
!
!      DECHARGE
        else
!
!        EXAMEN DE LA SOLUTION ELASTIQUE (XI=0)
            dp = 0.d0
            xi = 0.d0
            if (vim(9) .ne. 0.d0) then
                call nmtasp(ndimsi, carcri, mat, sigel, vim, &
                            epm, dp, sp, xi, f, &
                            iret)
                ASSERT(iret .eq. 0)
            end if
!
!        DECHARGE CLASSIQUE
            if (vim(9) .eq. 0.d0 .or. f .le. 0.d0) then
                call nmtaac(0, ndimsi, mat, sigel, vim, &
                            epm, dp, sp, xi, sigp, &
                            vip)
!
!
!
!        PSEUDO DECHARGE  -> EQUATION SCALAIRE F=0 AVEC  0<T<1
            else
!
!          CALCUL DE XI : EQUATION SCALAIRE F=0 AVEC  0 < XI < 1
                call nmtaxi(ndimsi, carcri, mat, sigel, vim, &
                            epm, dp, sp, xi, g, &
                            iret)
!
!          ITERATIONS DE NEWTON
                if (g .gt. 0.d0) then
!
                    dp = 0.d0
                    sp = vim(2)
                    do niter = 1, int(carcri(1))
!
!              DIRECTION DE DESCENTE
                        call nmtacr(3, ndimsi, mat, sigel, vim, &
                                    epm, dp, sp, xi, f, &
                                    g, fds, gds, fdp, gdp, &
                                    fdx, gdx, dpmax, sig, tang)
                        det = fdx*gds-fds*gdx
                        dirxi = (g*fds-f*gds)/det
                        dirsp = (f*gdx-g*fdx)/det
                        dirdp = 0.d0
!
!              CORRECTION DIRECTION POUR RESTER DANS 0<XI<1 ET SP-<SP<S
                        if (xi+rhomax*dirxi .lt. 0.d0) dirxi = (-xi)/rhomax
                        if (xi+rhomax*dirxi .gt. 1.d0) dirxi = (1.d0-xi)/rhomax
                        if (sp+rhomax*dirsp .lt. vim(2)) dirsp = (vim(2)-sp)/rhomax
                        if (sp+rhomax*dirsp .gt. mat(11)) dirsp = (mat(11)-sp)/rhomax
!
!              RECHERCHE LINEAIRE
                        ener = (f**2+g**2)/2.d0
                        min = (f*fdx+g*gdx)*dirxi+(f*fds+g*gds)*dirsp
                        rho = rhomax*interi
                        call nmtarl(3, ndimsi, mat, sigel, vim, &
                                    epm, dp, sp, xi, dirdp, &
                                    dirsp, dirxi, min, rho, ener)
!
!              ACTUALISATION
                        xi = xi+rho*dirxi
                        sp = sp+rho*dirsp
!
                        if (ener/mat(4)**2 .lt. carcri(3)**2) goto 310
                    end do
                    iret = 1
                    goto 999
310                 continue
                end if
!
                call nmtaac(1, ndimsi, mat, sigel, vim, &
                            epm, dp, sp, xi, sigp, &
                            vip)
!
!
            end if
!
        end if
!
    end if
!
!
! -- RIGIDITE TANGENTE (FULL_MECA)
!
    if (option(1:9) .eq. 'FULL_MECA') then
        ind = int(vip(9)+0.5d0)
        call nmtari(ind, ndimsi, mat, sigel, vim, &
                    epm, dp, sp, xi, dsidep)
    end if
!
!
! -- RIGIDITE TANGENTE (RIGI_MECA_TANG) -> MATRICE ELASTIQUE
!
    if (option(1:14) .eq. 'RIGI_MECA_TANG') then
        ind = 0
        dp = 0.d0
        sp = vim(2)
        xi = 1.d0
        call nmtari(ind, ndimsi, mat, sigel, vim, &
                    epm, dp, sp, xi, dsidep)
    end if
!
!
! REMISE AUX NORMES
    if (option(1:14) .eq. 'RIGI_MECA_TANG') then
        do k = 4, ndimsi
            vim(2+k) = vim(2+k)/rac2
        end do
    else
        do k = 4, ndimsi
            vim(2+k) = vim(2+k)/rac2
            vip(2+k) = vip(2+k)/rac2
        end do
    end if
!
!
999 continue
! FIN ------------------------------------------------------------------
end subroutine
