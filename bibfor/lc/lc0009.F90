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
! aslint: disable=W0104
!
subroutine lc0009(fami, kpg, ksp, ndim, imate, &
                  instam, instap, epsm, &
                  deps, sigm, vim, option, &
                  sigp, vip, typmod, &
                  nvi, dsidep, codret)
!
    implicit none
!
#include "asterfort/lclbr1.h"
#include "asterfort/Behaviour_type.h"
!
    integer(kind=8) :: kpg, ksp, ndim, imate
    integer(kind=8) :: codret
    character(len=*) :: fami
    character(len=8) :: typmod(*)
    character(len=16) :: option
    integer(kind=8) :: nvi
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(*), deps(*)
    real(kind=8) :: sigm(*), vim(nvi), sigp(*), vip(nvi), dsidep(*)

!       IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!       IN      KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!       IN      NDIM    DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!       IN      TYPMOD  : TYPE DE MODELISATION
!               IMATE    ADRESSE DU MATERIAU CODE
!               COMPOR    COMPORTEMENT DE L ELEMENT
!                     COMPOR(1) = RELATION DE COMPORTEMENT (CHABOCHE...)
!                     COMPOR(2) = NB DE VARIABLES INTERNES
!                     COMPOR(3) = TYPE DE DEFORMATION (PETIT,JAUMANN...)
!               CRIT    CRITERES  LOCAUX
!                       CRIT(1) = NOMBRE D ITERATIONS MAXI A CONVERGENCE
!                                 (ITER_INTE_MAXI == ITECREL)
!                       CRIT(2) = TYPE DE JACOBIEN A T+DT
!                                 (TYPE_MATR_COMP == MACOMP)
!                                 0 = EN VITESSE     > SYMETRIQUE
!                                 1 = EN INCREMENTAL > NON-SYMETRIQUE
!                       CRIT(3) = VALEUR DE LA TOLERANCE DE CONVERGENCE
!                                 (RESI_INTE == RESCREL)
!                       CRIT(5) = NOMBRE D'INCREMENTS POUR LE
!                                 REDECOUPAGE LOCAL DU PAS DE TEMPS
!                                 (ITER_INTE_PAS == ITEDEC)
!                                 0 = PAS DE REDECOUPAGE
!                                 N = NOMBRE DE PALIERS
!               INSTAM   INSTANT T
!               INSTAP   INSTANT T+DT
!               DEPS   INCREMENT DE DEFORMATION TOTALE
!               SIGM    CONTRAINTE A T
!               VIM    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!    ATTENTION  VIM    VARIABLES INTERNES A T MODIFIEES SI REDECOUPAGE
!               OPTION     OPTION DE CALCUL A FAIRE
!                             'RIGI_MECA_TANG'> DSIDEP(T)
!                             'FULL_MECA'     > DSIDEP(T+DT) , SIG(T+DT)
!                             'RAPH_MECA'     > SIG(T+DT)
!       OUT     SIGP    CONTRAINTE A T+DT
!               VIP    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSIDEP    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
! ======================================================================

    call lclbr1(fami, kpg, ksp, ndim, typmod, &
                imate, epsm, deps, &
                option, sigp, dsidep)
end subroutine
