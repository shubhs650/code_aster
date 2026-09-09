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
subroutine lcplas(BEHInteg, &
                  fami, kpg, ksp, relaComp, toler, &
                  itmax, mod, imat, nmat, materd, &
                  materf, nr, nvi, timed, timef, &
                  deps, epsd, sigd, vind, sigf, &
                  vinf, comp, nbcomm, cpmono, pgl, &
                  nfs, nsg, toutms, hsr, &
                  codret, theta, vp, seuil, &
                  devg, devgii, drdy, carcri)
!
    use Behaviour_type
    implicit none
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
!     INTEGRATION IMPLICITE DES COMPORTEMENTS. CALCUL DE SIGF,VINF,DSDE
!     ----------------------------------------------------------------
!     ARGUMENTS
!
!     IN FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!        KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!        LOI    :  MODELE DE COMPORTEMENT
!        TOLER  :  TOLERANCE DE CONVERGENCE LOCALE
!        ITMAX  :  NOMBRE MAXI D'ITERATIONS LOCALES
!        MOD    :  TYPE DE MODELISATION
!        IMAT   :  ADRESSE DU MATERIAU CODE
!        NMAT   :  DIMENSION MATER
!        MATERD :  COEFFICIENTS MATERIAU A T
!        MATERF :  COEFFICIENTS MATERIAU A T+DT
!        NR     :  NB EQUATION DU SYSTEME R(DY)
!        NVI    :  NB VARIABLES INTERNES
!        TIMED  :  INSTANT  T
!        TIMEF  :  INSTANT T+DT
!        DEPS   :  INCREMENT DE DEFORMATION
!        EPSD   :  DEFORMATION A T
!        SIGD   :  CONTRAINTE A T
!        VIND   :  VARIABLES INTERNES A T
!        COMP   :  COMPOR - LOI ET TYPE DE DEFORMATION
!        NBCOMM :  INCIDES DES COEF MATERIAU monocristal
!        CPMONO :  NOM DES COMPORTEMENTS monocristal
!        PGL    :  MATRICE DE PASSAGE
!        TOUTMS :  TENSEURS D'ORIENTATION monocristal
!        HSR    :  MATRICE D'INTERACTION monocristal
!        THETA  :  PARAMETRE DE LA THETA-METHODE
!        VP     :  VALEURS PROPRES DU DEVIATEUR ELASTIQUE(HOEK-BROWN)
!        VECP   :  VECTEURS PROPRES DU DEVIATEUR ELASTIQUE(HOEK-BROWN)
!        CRIT   : CRITERES DE CONVERGENCE LOCAUX
!
!    OUT SIGF   :  CONTRAINTE A T+DT
!        VINF   :  VARIABLES INTERNES A T+DT
!        CODRET :  CODE RETOUR. 0=OK, 1=ECHEC
!        DRDY   :  MATRICE JACOBIENNE
!       ----------------------------------------------------------------
#include "asterfort/Behaviour_type.h"
#include "asterfort/lchobr.h"
#include "asterfort/lcpllg.h"
#include "asterfort/lcplnl.h"
#include "asterfort/lcrous.h"
    integer(kind=8) :: itmax, codret, irtet, kpg, ksp
    integer(kind=8) :: imat, nmat, nvi, nr
!
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: timed, timef, deltat
    real(kind=8) :: toler, theta
    real(kind=8) :: epsd(6), deps(6)
    real(kind=8) :: sigd(6), sigf(6)
    real(kind=8) :: vind(*), vinf(*)
    real(kind=8) :: materf(nmat, 2), materd(nmat, 2)
    real(kind=8) :: seuil, devg(*), devgii
    real(kind=8) :: vp(3), drdy(nr, nr)
!
    character(len=8) :: mod
    character(len=16) :: relaComp
!
    integer(kind=8) :: nbcomm(nmat, 3), nfs, nsg, cutLevel
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg)
    character(len=16) :: comp(*)
    character(len=24) :: cpmono(5*nmat+1)
    character(len=*) :: fami
!       ----------------------------------------------------------------
!
    codret = 0
    deltat = timef-timed
    cutLevel = BEHInteg%behavPara%cutLevel
!
!       ----------------------------------------------------------------
!       CAS PARTICULIERS
!       ----------------------------------------------------------------
!
    if (relaComp(1:8) .eq. 'ROUSS_PR' .or. relaComp(1:10) .eq. 'ROUSS_VISC') then
        call lcrous(fami, kpg, ksp, toler, itmax, &
                    imat, nmat, materd, materf, nvi, &
                    deps, sigd, vind, theta, relaComp, &
                    deltat, sigf, vinf, irtet)
        if (irtet .gt. 0) goto 1
!
    elseif ((relaComp(1:10) .eq. 'HOEK_BROWN') .or. (relaComp(1:14) .eq. 'HOEK_BROWN_EFF')) then
        call lchobr(toler, itmax, nmat, materf, &
                    vind, &
                    vp, cutLevel, sigf, &
                    vinf, irtet)
        if (irtet .gt. 0) goto 1
!
    else if (relaComp(1:6) .eq. 'LAIGLE') then
        call lcpllg(toler, itmax, mod, nmat, materf, &
                    nr, nvi, deps, sigd, vind, &
                    seuil, cutLevel, sigf, vinf, devg, &
                    devgii, irtet)
        if (irtet .gt. 0) goto 1
!
!       ----------------------------------------------------------------
!       CAS GENERAL : RESOLUTION PAR NEWTON
!       ----------------------------------------------------------------
    else
        call lcplnl(BEHInteg, &
                    fami, kpg, ksp, relaComp, toler, &
                    itmax, mod, imat, nmat, materd, &
                    materf, nr, nvi, timed, timef, &
                    deps, epsd, sigd, vind, comp, &
                    nbcomm, cpmono, pgl, nfs, nsg, &
                    toutms, hsr, sigf, vinf, &
                    irtet, drdy, carcri)
        if (irtet .eq. 1) then
            goto 1
        else if (irtet .eq. 2) then
            goto 2
        end if
!
    end if
!
!     CONVERGENCE OK
!
    codret = 0
    goto 999
!
1   continue
!
!     PB INTEGRATION ou ITMAX ATTEINT : redecoupage local puis global
    codret = 1
    goto 999
!
2   continue
!     ITMAX ATTEINT : redecoupage du pas de temps global
    codret = 2
!
999 continue
end subroutine
